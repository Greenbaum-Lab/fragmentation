import numpy as np
import math


class Coalescence:
    def __init__(self, matrix: np.ndarray) -> None:
        """
        Initialize a coalescence times matrix object
        :param matrix: input Coalescence time matrix
        """
        self.matrix = matrix
        self.shape = matrix.shape[0]

    def produce_fst(self) -> np.ndarray:
        """
        produces and returns the corresponding Fst matrix
        :return: The corresponding Fst matrix
        """
        F_mat = np.zeros((self.shape, self.shape))
        for i in range(self.shape):
            for j in range(i + 1, self.shape):
                t_S = (self.matrix[i, i] + self.matrix[j, j]) / 2
                t_T = (self.matrix[i, j] + t_S) / 2
                if np.isinf(t_T):
                    F_i_j = 1
                else:
                    F_i_j = (t_T - t_S) / t_T
                F_mat[i, j], F_mat[j, i] = F_i_j, F_i_j
        return F_mat


class Migration:
    def __init__(self, matrix: np.ndarray) -> None:
        """
        Initialize a migration matrix object
        :param matrix: input migration matrix
        """
        self.matrix = matrix
        self.shape = matrix.shape[0]

    def produce_coalescence(self) -> np.ndarray:
        """
        produces and returns the corresponding coalescence matrix
        :return: The corresponding coalescence matrix
        """
        A = self.produce_coefficient_matrix()
        b = self.produce_solution_vector()
        x = np.linalg.solve(A, b)
        T_mat = np.zeros((self.shape, self.shape))
        cur_ind = 0
        for i in range(self.shape):
            for j in range(i, self.shape):
                T_mat[i, j] = x[cur_ind]
                T_mat[j, i] = x[cur_ind]
                cur_ind += 1
        return T_mat

    def calculate_first_coefficients(self, j: int, i: int, same_pop: int, lower_bound: int, upper_bound: int,
                                     p_list: list, counter: list) -> float:
        """
        calculates the coefficients for the first n equations
        :param j: column of coefficient matrx
        :param i: row of coefficient matrix
        :param same_pop: The column corresponding to T(i,i)
        :param lower_bound: j values in range [lower_bound, upper_bound] correspond to coefficient -M(i,i+counter)
        :param upper_bound: j values in range [lower_bound, upper_bound] correspond to coefficient -M(i,i+counter)
        :param p_list: all values that are smaller than i
        :param counter: counts the number of time j was in the interval [lower_bound,upper_bound]
        :return: The coefficient in place [i,j], for i in [n-1]
        """
        n = self.matrix.shape[0]
        if j == same_pop:
            return 1 + np.sum(self.matrix[i, :])
        if lower_bound <= j <= upper_bound:
            counter[0] += 1
            return -1 * self.matrix[i, i + counter[0] - 1]
        for p in p_list:
            if j == (i - p) + np.sum([n - k for k in range(p)]):
                return -1 * self.matrix[i, p]
        return 0

    def calculate_last_coefficients(self, j, cur_pop, other_pop) -> float:
        """
        calculates the coefficients for the last (n choose 2) rows in the coefficient matrx
        :param j: the column in the coefficient matrix
        :param cur_pop: the index of the population that corresponds to the current value
        :param other_pop: the index of the other population that corresponds to the current value
        :return: The coefficient in the coefficient matrix according to certain conditions deduced from
        Wilkinson-Herbots' equations.
        """
        n = self.matrix.shape[0]
        if j == np.sum([n - k for k in range(other_pop)]) + (cur_pop - other_pop):
            return float(np.sum(self.matrix[[cur_pop, other_pop], :]))
        for p in range(n):
            for t in [other_pop, cur_pop]:
                if t == other_pop:
                    not_t = cur_pop
                else:
                    not_t = other_pop
                if p != not_t:
                    min_t_p = min(t, p)
                    max_t_p = max(t, p)
                    if j == np.sum([n - k for k in range(min_t_p)]) + max_t_p - min_t_p:
                        return -1 * self.matrix[not_t, p]
        return 0

    def produce_coefficient_matrix(self) -> np.ndarray:
        """
        produce and return the coefficient matrix used to calculate the T matrix(coalescence).
        :return: Coefficient matrix corresponding to object's migration matrix
        """
        n = self.shape
        n_last_equations = comb(n, 2)
        n_first_equations = n
        mat_size = n_first_equations + n_last_equations
        coefficient_mat = np.zeros((mat_size, mat_size))
        for i in range(n_first_equations):
            same_population = int(np.sum([n - k for k in range(i)]))
            lower_bound = same_population + 1
            upper_bound = np.sum([n - k for k in range(i + 1)]) - 1
            smaller_ind_lst = [p for p in range(i)]
            counter = [1]
            for j in range(mat_size):
                coefficient_mat[i, j] = self.calculate_first_coefficients(j, i, same_population, lower_bound,
                                                                          upper_bound, smaller_ind_lst, counter)
        cur_population = 1
        other_population = 0
        for i in range(n_last_equations):
            if other_population == cur_population:
                other_population = 0
                cur_population += 1
            for j in range(mat_size):
                coefficient_mat[n + i, j] = self.calculate_last_coefficients(j, cur_population, other_population)
            other_population += 1

        return coefficient_mat

    def produce_solution_vector(self):
        """
        produce the solution vector(b), according to Wilkinson-Herbot's equations
        :return: solution vector b
        """
        n = self.shape
        n_first = np.repeat(1, n)
        n_last = np.repeat(2, comb(n, 2))
        return np.hstack((n_first, n_last))


def comb(n: int, k: int) -> int:
    """
    calculate and return n Choose k
    :param n: number of objects
    :param k: number of selected objects
    :return: n Choose k
    """
    return int(math.factorial(n) / (math.factorial(k) * math.factorial(n - k)))


def m_to_f(m: np.ndarray) -> np.ndarray:
    """
    Receives a migration matrix(a squared, positive matrix with zeroes on the diagonal), and returns its
    corresponding Fst matrix according to Wilkinson-Herbot's equations. Note that matrices that represent migration
    networks where one or more of the subpopulations are isolated will cause an error, since their coalescence times
    are interpreted as infinite, and Fst as 0. In this case, use the migration matrix that doesn't include the isolated
    subpopulations.
    :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
    :return: Corresponding Fst matrix according to Wilkinson-Herbot's equations. If there is no solution, an error will
    occur.
    """
    M = Migration(m)
    t = M.produce_coalescence()
    T = Coalescence(t)
    return T.produce_fst()


def m_to_t(m: np.ndarray) -> np.ndarray:
    """
       Receives a migration matrix (a squared, positive matrix with zeroes on the diagonal), and returns its
       corresponding Coalescent times (T) matrix according to Wilkinson-Herbot's equations.
       Note that matrices that represent migration networks where one or more of the subpopulations are isolated
       will cause an error, since their coalescence times are interpreted as infinite, and Fst as 1.
       In this case, use the migration matrix that doesn't include the isolated subpopulations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return: Corresponding T matrix according to Wilkinson-Herbot's equations. If there is no solution,
       an error will occur.
       """
    M = Migration(m)
    return M.produce_coalescence()


