import ctypes

import numpy as np
import math
from collections import deque


lib = ctypes.cdll.LoadLibrary('./libmigration.so')

lib.coefficient_matrix_from_migration.restype = ctypes.POINTER(ctypes.c_double)
lib.coefficient_matrix_from_migration.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]

def coefficient_matrix_from_migration_wrapper(migration_matrix):
    n = migration_matrix.shape[0]
    mat_size = n + (n * (n - 1)) // 2  # size of the coefficient matrix
    migration_matrix_c = migration_matrix.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    result_c = lib.coefficient_matrix_from_migration(migration_matrix_c, n)
    result = np.ctypeslib.as_array(result_c, shape=(mat_size*mat_size,)).reshape((mat_size, mat_size))

    return result

lib.coalescence_from_migration.restype = ctypes.POINTER(ctypes.c_double)
lib.coalescence_from_migration.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]




# This is a file for other users to copy to their projects.

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

    # def produce_coalescence(self) -> np.ndarray:
    #     """
    #     produces and returns the corresponding coalescence matrix
    #     :return: The corresponding coalescence matrix
    #     """
    #     A = self.produce_coefficient_matrix()
    #     b = self.produce_solution_vector()
    #     x = np.linalg.solve(A, b)
    #     T_mat = np.zeros((self.shape, self.shape))
    #     cur_ind = 0
    #     for i in range(self.shape):
    #         for j in range(i, self.shape):
    #             T_mat[i, j] = x[cur_ind]
    #             T_mat[j, i] = x[cur_ind]
    #             cur_ind += 1
    #     return T_mat

    def produce_coalescence(self) -> np.ndarray:
        n = self.matrix.shape[0]
        migration_matrix_c = self.matrix.flatten().ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        result_c = lib.coalescence_from_migration(migration_matrix_c, n)
        T_mat = np.ctypeslib.as_array(result_c, shape=(n * n,)).reshape((n, n))
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


def find_fst(m: np.ndarray) -> np.ndarray:
    """
    Receives a migration matrix with one connected component(a squared, positive matrix with zeroes on the diagonal),
    and returns it's corresponding Fst matrix according to Wilkinson-Herbot's equations.
    :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
    :return: Corresponding Fst matrix according to Wilkinson-Herbot's equations. If there is no solution, an error will
    occur.
    """
    if m.shape[0] == 1:
        return np.array([[0]])
    M = Migration(m)
    t = M.produce_coalescence()
    T = Coalescence(t)
    return T.produce_fst()


def find_coalescence(m: np.ndarray) -> np.ndarray:
    """
       Receives a migration matrix with one connected component
       (a squared, positive matrix with zeroes on the diagonal), and returns it's corresponding Coalescent times
       (T) matrix according to Wilkinson-Herbot's equations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return: Corresponding T matrix according to Wilkinson-Herbot's equations. If there is no solution,
       an error will occur.
       """
    if m.shape[0] == 1:
        return np.array([[1]])
    M = Migration(m)
    return M.produce_coalescence()


def find_components(matrix: np.ndarray) -> dict:
    """
    Find connected components in a connected graph represented by adjacency matrix.
    :param matrix: adjacency matrix representing a directed graph
    :return:something
    """
    components = 1
    n = matrix.shape[0]
    queue = deque()
    visited = set()
    not_visited = set([i for i in range(1, n)])
    visited.add(0)
    comp_dict = {components: [0]}
    queue.append(0)
    while len(not_visited) != 0:
        while len(queue) != 0:
            cur_vertex = queue.popleft()
            for i in range(n):
                if i not in visited and (matrix[cur_vertex, i] != 0 or matrix[i, cur_vertex] != 0):
                    queue.append(i)
                    visited.add(i)
                    not_visited.remove(i)
                    comp_dict[components].append(i)
        for vertex in not_visited:
            components += 1
            queue.append(vertex)
            visited.add(vertex)
            not_visited.remove(vertex)
            comp_dict[components] = [vertex]
            break
    return comp_dict


def split_migration_matrix(migration_matrix: np.ndarray, connected_components: list) -> list:
    """
    Splits a migration matrix to sub-matrices according to it's connected components.
    :param migration_matrix: A valid migration matrix.
    :param connected_components: list of lists, where each list represents a connected component's vertices
    (populations).
    :return: A list of sub-matrices, where each sun-matrix is the migration matrix of a connected component. Note that
    in order to interpret which populations are described in each sub matrix the connected components list is needed.
    """
    sub_matrices = []
    for component in connected_components:
        sub_matrix = migration_matrix[np.ix_(component, component)]
        sub_matrices.append(sub_matrix)

    return sub_matrices


def split_migration(migration_matrix: np.ndarray) -> tuple:
    """
    Finds a migration matrix connected components, and splits the matrix to it's connected components.
    :param migration_matrix: A valid migration matrix.
    :return: A tuple (sub_matrices, components). Sub matrices is a list of numpy arrays, where each array is a
    component's migration matrix. components is a list of lists, where each list represents a component vertices
    (populations). The order of the components corresponds to the order of the sub-matrices.
    """
    components = list(find_components(migration_matrix).values())
    sub_matrices = split_migration_matrix(migration_matrix, components)
    return sub_matrices, components


def reassemble_matrix(sub_matrices: list, connected_components: list, which: str) -> np.ndarray:
    num_nodes = sum(len(component) for component in connected_components)
    if which == "fst":
        adjacency_matrix = np.ones((num_nodes, num_nodes), dtype=float)
    else:
        adjacency_matrix = np.full((num_nodes, num_nodes), np.inf)

    for component, sub_matrix in zip(connected_components, sub_matrices):
        indices = np.array(component)
        adjacency_matrix[np.ix_(indices, indices)] = sub_matrix

    return adjacency_matrix


def m_to_f(m: np.ndarray) -> np.ndarray:
    """
    Receives a migration matrix(a squared, positive matrix with zeroes on the diagonal) with any number
    of connected components, and returns it's corresponding Fst matrix according to Wilkinson-Herbot's equations and
    Slatkin equations.
    :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
    :return: Corresponding Fst matrix according to Wilkinson-Herbot's equations. If there is no solution, an error will
    occur.
    """
    split = split_migration(m)
    sub_matrices, components = split[0], split[1]
    f_matrices = []
    for matrix in sub_matrices:
        f_matrices.append(find_fst(matrix))
    return reassemble_matrix(f_matrices, components, "fst")


def m_to_t(m: np.ndarray) -> np.ndarray:
    """
       Receives a migration matrix(a squared, positive matrix with zeroes on the diagonal) with any number
       of connected components, and returns it's corresponding Coalescent times (T) matrix according to
       Wilkinson-Herbot's equations.
       :param m: Migration matrix- squared, positive, with zeroes on the diagonal.
       :return: Corresponding T matrix according to Wilkinson-Herbot's equations. If there is no solution,
       an error will occur.
       """
    split = split_migration(m)
    sub_matrices, components = split[0], split[1]
    t_matrices = []
    for matrix in sub_matrices:
        t_matrices.append(find_coalescence(matrix))
    return reassemble_matrix(t_matrices, components, "coalescence")
