# 

**Fragmentation** is a Python-based toolkit developed for analyzing the genetic effects of reduced connectivity in population networks. The repository provides streamlined scripts for processing network and genetic data, computing a variety of statistics, and generating visualizations, with optional C acceleration for computationally intensive tasks.

## Key Components

- **`Transformation.py`**  
  Implements mathematical transformations between migration matrices, coalescence times, and Fst statistics, including routines for component detection and conservative migration matrix generation, using both Python and optional C acceleration.

- **`centrality.py`**  
  Computes node centralities (degree, betweenness) for single or multiple networks, aiding network fragmentation analysis.

- **`centrality_corr.py`**  
  Analyzes correlations between node centrality measures and heterozygosity, with tools for merging, filtering, and plotting these relationships across fragmentation experiments.

- **`distance_matrices.py`**  
  Calculates shortest path, Euclidean, and random walk distance matrices for network nodes using NetworkX and parallel processing.

- **`distributions.py`**  
  Extracts, filters, and visualizes distributions of genetic diversity measures (e.g., heterozygosity, Fst) at fixed fragmentation intervals.

- **`early_warning.py`**  
  Computes early warning indicators for genetic collapse or transitions (e.g., standard deviation, skewness, kurtosis, return rate) and provides plotting utilities.

- **`funcs.py`**  
  Defines data structures and utility functions for loading, processing, and summarizing fragmentation simulation results, including node/step annotation and component analysis.

- **`funcs_initial_data.py`**  
  Generates and normalizes initial networks, runs fragmentation replicates, and summarizes genetic statistics for various network models and fragmentation types.

- **`giant_comp.py`**  
  Analyzes the relationship between the size of the largest network component and genetic diversity, including binning and plotting tools.

- **`libmigration.c` / `libmigration.so`**  
  Supplies C routines and a compiled library for efficient migration, coalescence, and Fst matrix computations, used as an optional backend for high-performance operations.


## License
This repository is provided for academic and research purposes. Please see the repository or contact the authors for licensing details.

## Contact
For questions, suggestions, or contributions, please contact ohad.peled@mail.huji.ac.il
If you need further customization or more detailed usage examples, just ask!
