Networks fragmentation analysis

Code for creating multiple population networks and evaluate the genetic diversity and genetic differntiation of populations

1. create a well connected network (ER, RGG)
2. create a fragmentation process (random, correlated, distance dependent)
3. calculate pairwise Fst for all population for each step along the fragmentation 
4. calculate coalesence times for all populations for each step along the fragmentation (heterozygosity is the diagonal)
5. transform list of lists to dataframe of results
6. plot distributions
7. plot average, median






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
