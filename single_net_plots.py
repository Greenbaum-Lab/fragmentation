



#############single net plots
from funcs3 import make_fragmentation

n = 50  # no. of nodes
p = 0.4  # probability to connect nodes
# net = nx.erdos_renyi_graph(n=n, p=p)  # create network
# net = nx.erdos_renyi_graph(n=n,p=0.8)
net = nx.random_geometric_graph(n, p)


x = make_fragmentation(net=net, frag_type='rand', ignore_isolated=False)
print(x[4][['fst', 'step']])

plt.plot(x[5]['step'], x[5]['avg'], label="avg rand", color=color_palette(0))
plt.show()


plt.figure()
joyplot(
    data=x[4][['fst', 'step']],
    by='step', ylim=0, overlap=0.5,
    colormap=plt.cm.autumn, fade=True,
    figsize=(12, 8)
)
plt.title('pairwise Fst along distance fragmentation', fontsize=16, pad=-20, y=1.02, verticalalignment='bottom')
plt.show()


