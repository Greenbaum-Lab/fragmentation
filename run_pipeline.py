import pickle
import pandas as pd
from funcs_initial_data import make_networks, make_replicates_new, make_rgg

######### run pipeline
print("here i start!")
n_nodes = 50  # no. of nodes
n_rep = 100
n_edges = 250
net = "RGG"
ignore = False

# # # create list off nets
nets =  make_rgg(n_nets=n_rep, n_nodes=n_nodes, target_edges=n_edges)

# run the pipeline for all fragmentation types
rand = make_replicates_new(nets=nets, frag_type='rand', ignore=ignore)
print("1")
pickle_filename = f'{net}, rand_ignore_{ignore}.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)
del rand

cor = make_replicates_new(nets=nets, frag_type='cor', ignore=ignore)
print("2")
with open(pickle_filename, 'wb') as file:
    pickle.dump(cor, file)
del cor

intr = make_replicates_new(nets=nets, frag_type='intr', ignore=ignore)
print("3")
with open(pickle_filename, 'wb') as file:
    pickle.dump(intr, file)
del intr

reg = make_replicates_new(nets=nets, frag_type='reg', ignore=ignore)
print("4")
with open(pickle_filename, 'wb') as file:
    pickle.dump(reg, file)
del reg

div = make_replicates_new(nets=nets, frag_type='div', ignore=ignore)
print("5")
with open(pickle_filename, 'wb') as file:
    pickle.dump(div, file)
del div

dist = make_replicates_new(nets=nets, frag_type='dist', ignore=ignore)
print("6")
with open(pickle_filename, 'wb') as file:
    pickle.dump(dist, file)
del dist

opt = make_replicates_new(nets=nets, frag_type='opt', ignore=ignore)
print("7")
with open(pickle_filename, 'wb') as file:
    pickle.dump(opt, file)
del opt

# opt2 = make_replicates_new(nets=nets, frag_type='opt2', ignore=ignore)
# print("8")
# wrst = make_replicates_new(nets=nets, frag_type='wrst', ignore=ignore)
# print("9")

print("finish")
########## finish pipeline