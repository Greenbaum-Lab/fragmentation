from funcs import load_data
from viz_funcs.viz_net import plot_het_component

#### plot heterozygosity-giant compoent
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
data = load_data(fragmentation_types)
plot_het_component(data, frag_types=fragmentation_types, n_bins=20, output='./figs/het_component.svg')



# ############ calculate and plot variance across nodes in the network
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
# data = load_data(fragmentation_types)
# #for single frag type use fragmentation_types[x]
# var = process_variance(data, fragmentation_types)
# plot_variance(var)