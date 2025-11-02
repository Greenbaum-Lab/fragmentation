import pandas as pd

from centrality_corr import compute_het_central_correlation, filter_correlations, plot_correlation, merge_centrality_het
from funcs import load_data
from centrality import compute_centrality_types

###scripts
##### compute centrality for all fragmentation types
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)

# centrality_df = compute_centrality_types(data, fragmentation_types)

#### merge centrality with heterozygosity data
# frag_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(frag_types)
# centrality_df = pd.read_csv('centrality_sig03.csv')
# merged_df = merge_centrality_het(centrality_df, data, frag_types)
# #
# #### compute correlation between centrality and heterozygosity
# centrality_df = pd.read_csv('centrality_het_sig03.csv')
# corr_df = compute_het_central_correlation(
#     df=centrality_df,
#     centrality='degree',
# )
#
### plot correlation between centrality and heterozygosity
corr_df = pd.read_csv('het_degree_correlation_sig03.csv')
filtered_corr_df = filter_correlations(corr_df, min_replicates=5)
print(filtered_corr_df)
plot_correlation(
    corr_df=filtered_corr_df,
    output_path='corr_degree_sig03.svg'
)

