














def calculate_mantel(net, fst_matrix, dist_type, perms):
    """
    Perform a Mantel test for two distance matrices.

    :param net: NetworkX graph
    :param fst_matrix: Genetic distance matrix
    :param dist_type: Type of distance ('euclidean' or 'path')
    :param perms: Number of permutations for the Mantel test
    :return: Weighted mean of the correlation and p-value by component size
    """
    # calculate the distance matrix based on network connectivity
    if dist_type == 'euclidean':
        distance_matrix = get_euclidean_matrix(net)
    if dist_type == 'path':
        distance_matrix = get_shortest_path_matrix(net)
    if dist_type == 'random':
        distance_matrix = get_random_walk_matrix(net)
    # if the network is connected, calculate the mantel test directly
    if nx.is_connected(net):
        r, p, _ = test(X=distance_matrix, Y=fst_matrix, perms=perms, method='pearson', ignore_nans=True)
        return r, p

    r_values, p_values, weights = [], [], []
    # in case the network is not connected, calculate the mantel test for each component
    for comp in find_connected_components(net):
        comp_dist_matrix = distance_matrix[np.ix_(comp, comp)]
        comp_fst_matrix = fst_matrix[np.ix_(comp, comp)]
        r, p, _ = test(X=comp_dist_matrix, Y=comp_fst_matrix, perms=perms, method='pearson', ignore_nans=True)
        r_values.append(r)
        p_values.append(p)
        if r is np.nan:
            weights.append(0)
        else:
            weights.append(len(comp))
    # Check if r_values or p_values is empty (happens when the network has its last component >3)
    if not r_values or not p_values:
        return None
    # symmetric 3*3 matrices get na in shortest path distance, so drop them in all distance matrices
    if len(r_values) == 1:
        return r_values[0], p_values[0]
    # Mask NaN values that result from symmetric matrices with no variation
    masked_r_values = np.ma.masked_array(r_values, np.isnan(r_values))
    masked_p_values = np.ma.masked_array(p_values, np.isnan(p_values))

    # calculate the weighted mean of the correlation and p-value
    weighted_r = np.ma.average(masked_r_values, weights=weights)
    weighted_p = np.ma.average(masked_p_values, weights=weights)

    return weighted_r, weighted_p























########### analysis of fst-distance

def calculate_mantel_for_process(data, perms, dist_type, replica):
    """"
    calculate mantel correlation and p value for each step along fragmentation.
    :param data: raw data of fragmentation type
    """
    results = []
    networks = access_networks(data)[replica]
    fst_matrices = access_fst_matrices(data)[replica]

    for step, (net, fst) in enumerate(zip(networks, fst_matrices)):
        result = calculate_mantel(net=net, perms=perms, fst_matrix=fst, dist_type=dist_type)
        if result is None:
            break
        r, p = result
        results.append({'step': step, 'r_val': r, 'p_val': p, 'replica': replica})
    cor_data = pd.DataFrame(results)
    return cor_data


def calculate_mantel_for_process(data, perms, dist_type, replica, num_workers=None):
    results = []
    networks = access_networks(data)[replica]
    fst_matrices = access_fst_matrices(data)[replica]

    def process_step(step, net, fst):
        result = calculate_mantel(net=net, perms=perms, fst_matrix=fst, dist_type=dist_type)
        if result is None:
            return None
        r, p = result
        return {'step': step, 'r_val': r, 'p_val': p, 'replica': replica}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_step, step, net, fst) for step, (net, fst) in enumerate(zip(networks, fst_matrices))]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)

    cor_data = pd.DataFrame(results)
    return cor_data


def calculate_mantel_replicas(data,perms, dist_type):
    """
    calculate mantel correlation and p value across fragmentation for all replicas.
    """
    results = []
    networks = access_networks(data)

    for replica in range(len(networks)):
        print(replica)
        cor_data = calculate_mantel_for_process(data, perms, dist_type=dist_type, replica=replica)
        results.append(cor_data)
    cor_data = pd.concat(results)

    return cor_data


def calculate_mantel_all(data, perms, dist_type='euclidean'):
    """"
    calculate mantel correlation and p value across fragmentation for all fragmnetation types.
    """
    results = []
    for frag_type in data.keys():
        print(frag_type)
        cor_data = calculate_mantel_replicas(data[frag_type], perms, dist_type)
        cor_data['fragmentation_type'] = frag_type
        results.append(cor_data)
        cor_data.to_csv(f'./csv/corl_fst_{frag_type}_{dist_type}.csv', index=False)

    cor_data = pd.concat(results)

    # write data as csv
    cor_data.to_csv(f'./TEST-corl_fst_{dist_type}_all.csv', index=False)
    return cor_data



### calculate mantel for all for euclidean
# fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)
# perms = 999
# cor_data = calculate_mantel_all(data, perms, dist_type='path')

def plot_cor_fst(df):
    """
    plot mantel correlation for a single type or replica.
    mini test.
    """
    # df = normalize_steps(df)
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='r_val', data=df)
    plt.xlabel('Fragmentation (%)')
    plt.ylabel('Correlation')
    plt.show()

def plot_mantel_all(df):
    """
    plot mantel correlation for all fragmentation types.
    """
    df['step'] = df['step'] / df['step'].max() * 100  # Normalize steps to percentage
    df = df[df['p_val'] < 0.05]  # Filter rows with pval < 0.05
    # Filter steps with at least 5 unique replicas

    df = (
        df
        .groupby(['fragmentation_type', 'step'])
        .filter(lambda g: g['replica'].nunique() >= 5)
    )
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='r_val', hue='fragmentation_type', data=df, errorbar='sd', legend=False)
    plt.xlabel('% fragmentation', fontsize=28)
    plt.ylabel('Correlation (r)', fontsize=28)
    plt.tick_params(axis='both', labelsize=25)
    plt.ylim(-0.05, 1.1)
    plt.savefig(f'./figs/cor_fst_euclidean_pval.svg', format="svg")
    plt.show()



#### plot single correlation fst-distance
# fragmentation_types = ['wrst']
# data = load_data(fragmentation_types)
# data = data[fragmentation_types[0]]
#
# steps = [0, 75, 150]
# fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
#
# for i, step in enumerate(steps):
#     net = data[1][5][step]
#     fst = data[7][5][step]
#     distance_matrix = get_random_walk_matrix(net)
#     r, p = calculate_mantel(net=net, fst_matrix=fst, dist_type='random', perms=999)
#     print(r, p)
#
#     flat_matrix1 = distance_matrix.flatten()
#     flat_matrix2 = fst.flatten()
#     flat_matrix1 = flat_matrix1[flat_matrix1 != 0]
#     flat_matrix2 = flat_matrix2[flat_matrix2 != 0]
#     df = pd.DataFrame({'distance': flat_matrix1, 'fst': flat_matrix2})
#     df = df.dropna()
#     # filter inf
#     df = df[~df['distance'].isin([np.inf, -np.inf])]
#     print(df)
#
#     sns.regplot(x='distance', y='fst', data=df, fit_reg=True, order=1, ax=axes[i])
#     axes[i].set_xlabel('Distance', fontsize=30)
#     axes[i].set_ylabel(r'Pairwise $F_{ST}$' if i == 0 else '', fontsize=30)
#     axes[i].tick_params(axis='both', labelsize=25)
#     axes[i].set_ylim(0, 0.5)
#     axes[i].text(0.05, 1.2, f'r={r:.2f}\np={p:.2e}', fontsize=20, transform=axes[i].transAxes)
#
# plt.tight_layout()
# plt.savefig(f'./figs/random_fst_steps.svg', format="svg")
# plt.show()

# df = pd.read_csv('./TEST-corl_fst_euclidean_all.csv')
# print(df)
# plot_mantel_all(df)

