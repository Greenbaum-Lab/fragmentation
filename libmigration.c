// libmigration.c – extended for asymmetric yet conservative migration matrices
// Compatible with Wilkinson‑Herbots (1998) structured‑coalescent equations.
// The solver now works with any migration matrix M that satisfies
//              Σ_j M_ij  =  Σ_i M_ij   for every deme i.
// ---------------------------------------------------------------------------
// Build:
//     gcc -O3 -fPIC -shared libmigration.c -lgsl -lgslcblas -o libmigration.so
// Optional debug:
//     gcc -O3 -fPIC -shared -DDEBUG_CONS libmigration.c -lgsl -lgslcblas -o libmigration.so
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>

/*-----------------------------------------------------------------------------*/
#define CONS_TOL 1e-6      /* maximum allowed |row‑sum − col‑sum| */
/*-----------------------------------------------------------------------------*/
/*  Helper: map unordered pair (i,j) with 0≤i≤j<n  to compact index 0..n(n+1)/2‑1 */
static inline int idx_T(int n, int i, int j)
{
    if (i > j) { int tmp = i; i = j; j = tmp; }
    return i * n - (i * (i - 1)) / 2 + (j - i);
}

/*-----------------------------------------------------------------------------*/
/*  Check conservativeness: row sums equal column sums to within tol            */
/*-----------------------------------------------------------------------------*/
static int is_conservative(const double *M, int n, double tol)
{
    double maxdiff = 0.0;
    for (int i = 0; i < n; ++i) {
        double row = 0.0, col = 0.0;
        for (int j = 0; j < n; ++j) {
            row += M[i * n + j];
            col += M[j * n + i];
        }
        double diff = fabs(row - col);
        if (diff > maxdiff) maxdiff = diff;
        if (diff > tol) {
#ifdef DEBUG_CONS
            printf("max diff %.3e exceeds tol %.3e\n", diff, tol);
#endif
            return 0;   /* not conservative */
        }
    }
#ifdef DEBUG_CONS
    printf("max diff %.3e (tol %.3e)\n", maxdiff, tol);
#endif
    return 1;           /* conservative */
}

/*=============================================================================*/
/*  Build coefficient matrix A (dim×dim) for WH Eq. 9 and RHS vector b         */
/*=============================================================================*/
static double *build_coefficient_matrix(const double *M, int n)
{
    const int dim = n + n * (n - 1) / 2;               /* unknowns = equations */
    double *A = calloc((size_t)dim * dim, sizeof(double));
    if (!A) return NULL;

    /* ---------------- 1.  Diagonal equations  ---------------------------- */
    for (int i = 0; i < n; ++i) {
        double Mi = 0.0;                                 /* Σ_k M_ik */
        for (int k = 0; k < n; ++k) Mi += M[i * n + k];

        A[i * dim + idx_T(n, i, i)] = 1.0 + Mi;          /* (1+Mi)·T_ii */

        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            int col = idx_T(n, i < k ? i : k, i < k ? k : i);
            A[i * dim + col] = -M[i * n + k];            /* -M_ik·T_ik */
        }
    }

    /* ---------------- 2.  Off‑diagonal equations  ------------------------ */
    int eq = n;  /* current row in A */
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j, ++eq) {
            double Mi = 0.0, Mj = 0.0;
            for (int k = 0; k < n; ++k) {
                Mi += M[i * n + k];
                Mj += M[j * n + k];
            }
            A[eq * dim + idx_T(n, i, j)] = Mi + Mj;      /* (Mi+Mj)·T_ij */

            /* - M_ik · T_kj */
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                int col = idx_T(n, (k < j ? k : j), (k < j ? j : k));
                A[eq * dim + col] -= M[i * n + k];
            }
            /* - M_jk · T_ik */
            for (int k = 0; k < n; ++k) {
                if (k == j) continue;
                int col = idx_T(n, (k < i ? k : i), (k < i ? i : k));
                A[eq * dim + col] -= M[j * n + k];
            }
        }
    }
    return A;
}

static double *build_rhs(int dim)
{
    double *b = malloc(dim * sizeof(double));
    if (!b) return NULL;
    for (int i = 0; i < dim; ++i) b[i] = 2.0;            /* WH Eq. 9 RHS = 1 */
    return b;
}

/*=============================================================================*/
/*  Public API – called from Python via ctypes                                   */
/*=============================================================================*/

double *coefficient_matrix_from_migration(double *M, int n)
{
    if (!is_conservative(M, n, CONS_TOL)) {
        fprintf(stderr, "Error: migration matrix is not conservative.\n");
        return NULL;
    }
    return build_coefficient_matrix(M, n);
}

/*  Solve A·x = b and return full symmetric T‑matrix (n×n)                      */
double *coalescence_from_migration(double *M, int n)
{
    if (!is_conservative(M, n, CONS_TOL)) {
        fprintf(stderr, "Error: migration matrix is not conservative.\n");
        return NULL;
    }

    const int dim = n + n * (n - 1) / 2;
    double *A = build_coefficient_matrix(M, n);
    double *b = build_rhs(dim);
    if (!A || !b) {
        fprintf(stderr, "Memory allocation failure.\n");
        return NULL;
    }

    gsl_matrix_view  Am = gsl_matrix_view_array(A, dim, dim);
    gsl_vector_view  bv = gsl_vector_view_array(b, dim);
    gsl_vector      *x  = gsl_vector_alloc(dim);
    gsl_permutation *p  = gsl_permutation_alloc(dim);
    if (!x || !p) {
        fprintf(stderr, "GSL allocation failure.\n");
        return NULL;
    }

    int signum;
    if (gsl_linalg_LU_decomp(&Am.matrix, p, &signum)) {
        fprintf(stderr, "GSL: LU decomposition failed (singular matrix).\n");
        return NULL;
    }
    gsl_linalg_LU_solve(&Am.matrix, p, &bv.vector, x);

    /* --- reconstruct symmetric T matrix --- */
    double *T = malloc((size_t)n * n * sizeof(double));
    if (!T) return NULL;

    int cur = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j, ++cur) {
            double val = gsl_vector_get(x, cur);
            T[i * n + j] = val;
            T[j * n + i] = val;   /* symmetry */
        }
    }

    /* tidy up */
    gsl_vector_free(x);
    gsl_permutation_free(p);
    free(A);
    free(b);

    return T;
}

