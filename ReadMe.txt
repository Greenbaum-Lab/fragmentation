
Networks fragmentation analysis

Code for creating multiple population networks and evaluate the genetic diversity and genetic differentiation of populations

1. create a well-connected network: ER, RGG
2. create fragmentation processes (random, correlated, distance-dependent)
    plot first and last network
    stop when network is no longer connected
3. create a list of migration networks
4. calculate the pairwise Fst for all populations for each step along the fragmentation
    excluding the diagonal
5. calculate the coalescence time for all populations for each step along the fragmentation
    heterozygosity is the diagonal
6. calculate betweenes and clustering average for each step
7. transform list of lists to dataframe of results
8. plot distributions of fst and heterozygosity
9. plot average, median across steps
10. plot average, median across centrality measures

run multiple iteration to get distributions
continue with broken network and track giant component until all network is broken
normalize heterozygosity so that all rows have equal sum






Transformations
*Written by: Eyal Haluts
*Email: eyal.haluts@mail.huji.ac.il

*Before running the code make sure the packages math and numpy are installed in your enviroment.

How to run the code:
- Move the file Transformation.py to your project.
- Import the function m_to_f from Transformation.py.
- import the function m_to_t from Transformation.py.
- Usuage of the functions is detailed in it's docstring.

Code example:

import numpy as np
from Transformation import m_to_f

M_1 = np.array([[0, 2, 0, 1], [0, 0, 1, 2], [2, 1, 0, 0], [1, 0, 2, 0]])
F_1 = m_to_f(M_1)
print(F_1, "\n")
M_2 = np.array([[0, 1.87, 1.48, 0.74], [0.65, 0, 1.74, 0.17], [1.73, 0, 0, 1.95], [1.7, 0.68, 0.46, 0]])
F_2 = m_to_f(M_2)
print(np.round(F_2, decimals=2))

........................................................................................................................................................................................
Prints:

[[0.         0.11111111 0.11111111 0.11111111]
 [0.11111111 0.         0.11111111 0.11111111]
 [0.11111111 0.11111111 0.         0.11111111]
 [0.11111111 0.11111111 0.11111111 0.        ]]

[[0.   0.1  0.09 0.1 ]
 [0.1  0.   0.12 0.13]
 [0.09 0.12 0.   0.1 ]
 [0.1  0.13 0.1  0.  ]]
