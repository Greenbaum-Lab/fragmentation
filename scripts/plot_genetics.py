import os

from data_manipulation.manp_genetics import filter_intervals, select_random_nodes, extract_nodes, process_variance
from funcs import  load_data
from viz_funcs.viz_genetics import plot_genetics, plot_distribution, plot_het_nodes, plot_variance

######################
###################### plot data
print(os.getcwd())

fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
fragmentation_types = ['rand']
data = load_data(fragmentation_types)
plot_genetics(data, measure='het')

#############plot distributions
#### one frag type each time
# fragmentation_types = ['rand']
# data = load_data(fragmentation_types)
# df = filter_intervals(data['rand'], measure='fst', interval_pct=25)
# plot_distribution(df,measure='fst', frag_type='rand')
#
# ###################### plot individual nodes
# ###############
# fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt', 'wrst']
# fragmentation_types = ['rand']
# data = load_data(fragmentation_types)
# df = data.get('rand').het_dist
# x=select_random_nodes(df, per_replica=1)
# df_selected = extract_nodes(df, x)
# plot_het_nodes(df_selected, n_nodes=10)
#
