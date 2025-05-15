
from funcs import FragmentationResult, percent_step, load_data, assign_node_numbers
from viz_genetics import plot_genetics, plot_distribution, plot_het_nodes, plot_variance


######################
###################### plot data
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
data = load_data(fragmentation_types)
plot_genetics(data, measure='het')

#############plot distributions
#### one frag type each time
fragmentation_types = ['rand']
data = load_data(fragmentation_types)
df = filter_intervals(data['rand'], measure='fst', interval_pct=25)
plot_distribution(df,measure='fst', frag_type='rand')

###################### plot individual nodes
###############
fragmentation_types = ['rand', 'cor', 'intr', 'dist', 'reg', 'div', 'opt', 'wrst']
fragmentation_types = ['rand']
data = load_data(fragmentation_types)
df = data.get('rand').het_dist
x=select_random_nodes(df, per_replica=1)
df_selected = extract_nodes(df, x)
plot_het_nodes(df_selected, n_nodes=10)

############ calculate and plot variance across nodes in the network
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
data = load_data(fragmentation_types)
#for single frag type use fragmentation_types[x]
var = process_variance(data, fragmentation_types)
plot_variance(var)