# Fragmentation

**Fragmentation** is a Python-based toolkit developed by Ohad Peled for analyzing the genetic effects of reduced connectivty in population networks. The repository provides streamlined scripts for simulating gene flow in (mostly) Random Geometric Graphs (RGG) and generating statistical summaries and visualizations.

## Key Components

- **`fragmentation.py`**  
  Processes input data (e.g., BAM or CSV files), calculates fragment-size distributions, and produces summary statistics. It includes options for filtering fragments by length and computing key metrics.

- **`plot_fragments.py`**  
  Generates visualizations such as histograms and density plots from the summary data, enabling a quick assessment of the quality and characteristics of sequencing libraries.

- **`utils.py`**  
  Contains helper functions for data parsing, filtering, and statistical calculations used by the analysis and plotting scripts.

