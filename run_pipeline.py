import pickle
from funcs_initial_data import run_replicates
from networks_generator import make_rgg

######### run pipeline
print("here i start!")
n_nodes = 50  # no. of nodes
n_rep = 100
n_edges = 250
net = "RGG"

# # # create list off nets
nets =  make_rgg(n_graphs=n_rep, n_nodes=n_nodes, target_edges=n_edges,asymmetric=True)
print("nets created")


# run the pipeline for all fragmentation types
rand = run_replicates(nets=nets, frag_key='rand',n_workers=20)
print("1")
pickle_filename = f'{net}, rand_asymmetric_TEST.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)
del rand

cor = run_replicates(nets=nets, frag_key='cor')
print("2")
with open(f'{net}, cor_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(cor, file)
del cor

intr = run_replicates(nets=nets, frag_key='intr')
print("3")
with open(f'{net}, intr_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(intr, file)
del intr

reg = run_replicates(nets=nets, frag_key='reg')
print("4")
with open(f'{net}, reg_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(reg, file)
del reg

div = run_replicates(nets=nets, frag_key='div')
print("5")
with open(f'{net}, div_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(div, file)
del div
#
dist = run_replicates(nets=nets, frag_key='dist')
print("6")
with open(f'{net}, dist_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(dist, file)
del dist
#
opt = run_replicates(nets=nets, frag_key='opt')
print("7")
with open(f'{net}, opt_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(opt, file)
del opt
#
wrst = run_replicates(nets=nets, frag_key='wrst')
print("8")
with open(f'{net}, wrst_asymmetric_TEST.pickle', 'wb') as file:
    pickle.dump(wrst, file)
del wrst
print("finish")
########## finish pipeline