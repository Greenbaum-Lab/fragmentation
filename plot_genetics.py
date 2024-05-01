from statistics import mean
import pickle
import pandas as pd
from joypy import joyplot
from matplotlib import pyplot as plt

from funcs_analysis import load_data, plot_data, filter_intervals, plot_all_distributions, plot_distribution

from funcs_initial_data import make_networks, make_replicates_new

# pd.set_option('display.max_rows',None)


########## run pipeline
print("here i start!")
n = 50  # no. of nodes
n_rep = 100  # no. of replicates
net = "RGG"
ignore = False
#
# # # create list off nets
nets = make_networks(n_nets=n_rep, n_nodes=n, net_type=net)

# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='rand', ignore=ignore)
print("1")
cor = make_replicates_new(nets=nets, frag_type='cor', ignore=ignore)
print("2")
int = make_replicates_new(nets=nets, frag_type='int', ignore=ignore)
print("3")
reg = make_replicates_new(nets=nets, frag_type='reg', ignore=ignore)
print("4")
div = make_replicates_new(nets=nets, frag_type='div', ignore=ignore)
print("5")
dist = make_replicates_new(nets=nets, frag_type='dist', ignore=ignore)
print("6")
opt = make_replicates_new(nets=nets, frag_type='opt', ignore=ignore)
print("7")
opt2 = make_replicates_new(nets=nets, frag_type='opt2', ignore=ignore)
print("8")
wrst = make_replicates_new(nets=nets, frag_type='wrst', ignore=ignore)
print("9")



# save files as tuple
pickle_filename = f'{net}, rand_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)

pickle_filename = f'{net}, cor_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(cor, file)

pickle_filename = f'{net}, int_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(int, file)

pickle_filename = f'{net}, reg_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(reg, file)

pickle_filename = f'{net}, div_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(div, file)

pickle_filename = f'{net}, dist_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(dist, file)

pickle_filename = f'{net}, opt_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(opt, file)

pickle_filename = f'{net}, opt2_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(opt2, file)

pickle_filename = f'{net}, wrst_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(wrst, file)
########## finish pipeline



fragmentation_types = ['rand', 'cor', 'int', 'dist', 'reg', 'div', 'opt', 'opt2', 'wrst']
net = 'RGG'
ignore = False
data = load_data(fragmentation_types, net, ignore)


# Plot fst and het along fragmentation
plot_data(data, 5, 'Pairwise Fst',measure='fst', save=False)
plot_data(data, 3, 'Heterozygosity',measure='heterozygosity', save=False)


##############plot distributions

# rand = filter_intervals(rand[2])
# plot_distribution(rand,'het','rand')
# plot_all_distributions(data)