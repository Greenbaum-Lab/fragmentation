
####### plotting of individual nodes heterozygosity for single fragmentation type #######

def select_random_nodes(
    df: pd.DataFrame,
    per_replica: int = 1,
    nodes_per_step: int = 50
) -> Dict[int, np.ndarray]:
    """
    For each replica, choose `per_replica` random node indices.

    :param df: DataFrame containing the heterozygosity data.
    :param per_replica: Number of random nodes to select per replica.
    :param nodes_per_step: Number of nodes in each step.
    :return: A dictionary with replica ids as keys and arrays of selected node indices as values.
    """
    # Ensure `node_number` is assigned
    df = assign_node_numbers(df, nodes_per_step)

    selections = {}
    for rep, sub in df.groupby("replica"):
        n_nodes = sub["node_number"].nunique()
        picks = np.random.choice(n_nodes, min(per_replica, n_nodes), replace=False)
        selections[int(rep)] = picks
    return selections


def extract_nodes(
    df: pd.DataFrame,
    selections: Dict[int, np.ndarray],
    nodes_per_step: int = 50
) -> pd.DataFrame:
    """
    Extract the heterozygosity data of selected nodes for each replica and step.

    :param df: DataFrame with 'node_number', 'step', 'replica', and 'het' values.
    :param selections: A dictionary with replicas as keys and lists of selected node indices as values.
    :param nodes_per_step: Number of nodes per step in the data (should match original assignment).
    :return: DataFrame containing only selected nodes, including a 'node_replica_id'.
    """
    out = []
    for rep, nodes in selections.items():
        sub = df[df["replica"] == rep]
        for node in nodes:
            node_df = sub[sub["node_number"] == node].copy()
            node_df["id"] = f"n{node}_r{rep}"
            out.append(node_df)

    return pd.concat(out, ignore_index=True).drop(columns=['replica', 'node_number'])





def plot_het_nodes(
    df: pd.DataFrame,
    n_nodes: int = 10,
) -> None:
    """
    Plot the heterozygosity for selected nodes across steps.

    :param df: DataFrame with 'step', 'node_replica_id', and 'het' values.
    :param n_nodes: Number of nodes to plot (choose top `n_nodes` nodes based on their node_replica_id).
    :param measure: The column to plot ('het' or 'fst').
    :param title: The plot's title.
    """
    node_ids = df['id'].unique()

    # Sample n_nodes
    selected_nodes = np.random.choice(node_ids, n_nodes, replace=False)

    df = percent_step(df, step_col='step', pct_col='step_pct')
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6,4))

    for node_id in selected_nodes:
        # Filter data for each node_replica_id
        node_data = df[df['id'] == node_id]

        # Plot the line for the node's data
        ax.plot(node_data['step_pct'], node_data['het'],color='grey', alpha=0.5)

    # Customize plot
    ax.set_xlabel('Time', fontsize=16)
    ax.set_ylabel("Heterozygosity", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.show()

