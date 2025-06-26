# Fragmentation

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





Fragmentation
Fragmentation is a Python toolkit for analyzing the genetic effects of reduced connectivity in population networks. It includes a suite of scripts to process network data, compute genetic and network statistics, and visualize results, with optional C acceleration for heavy computations.

Table of Contents
Overview
Repository Structure
Module Descriptions
Getting Started
License
Contact
Overview
This repository enables simulation and analysis of population network fragmentation and its genetic consequences, such as changes in heterozygosity and Fst. It provides tools for:

Running fragmentation experiments with various models
Calculating and visualizing genetic statistics
Analyzing network structure and centrality
Correlating genetic and network metrics

Module Descriptions
Transformation.py: Implements mathematical transformations between migration matrices, coalescence times, and Fst statistics, including routines for connected component detection and conservative migration matrix generation, optionally using C for speed.
centrality.py: Computes node centralities (degree and betweenness) for single or multiple networks, facilitating the analysis of network fragmentation effects.
centrality_corr.py: Analyzes correlations between node centralities and heterozygosity, providing tools to merge, filter, and plot these relationships across network fragmentation experiments.
distance_matrices.py: Calculates shortest path, Euclidean, and random walk distance matrices for network nodes, leveraging parallel processing for efficiency.
distributions.py: Extracts, filters, and visualizes the distributions of genetic diversity measures (e.g., heterozygosity, Fst) at fixed fragmentation intervals.
early_warning.py: Computes early warning indicators (e.g., standard deviation, skewness, kurtosis, return rate) of genetic collapse or transitions in network fragmentation, with plotting and export utilities.
funcs.py: Defines data structures (like FragmentationResult) and utility functions for loading, processing, and summarizing simulation runs, including node annotation and component analysis.
funcs_initial_data.py: Generates and normalizes initial networks, runs fragmentation replicates, and summarizes genetic statistics; supports various network models and fragmentation types.
giant_comp.py: Evaluates how the fraction of nodes in the largest network component (giant component) relates to genetic diversity, with binning and plotting functions.
libmigration.c / libmigration.so: Provides C routines for efficient matrix computations underlying migration, coalescence, and Fst calculations, used for performance-critical operations.

License
This repository is provided for academic and research purposes. Please see the repository or contact the authors for licensing details.

Contact
For questions, suggestions, or contributions, please contact ohad.peled@mail.huji.ac.il

If you need further customization or more detailed usage examples, just ask!
