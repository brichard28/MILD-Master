def mne_cluster_based_adjacency_nirs(info, pos, threshold_mm):
    """

    :param info:
    :param pos:
    :param threshold_mm:
    :return:
    """

    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse import csr_matrix
    import mne
    from matplotlib import pyplot as plt

    adjacency, ch_names = mne.channels.find_ch_adjacency(info, ch_type=None)
    adjacency_array = adjacency.toarray()  # Convert to a dense array if it is sparse

    nirs_positions_adjacency = pos.copy()[:, 0:2]
    nirs_positions_adjacency[:, 1] = nirs_positions_adjacency[:, 1] + 0.01

    # Threshold distance in millimeters (20 mm)
    distances = squareform(pdist(pos))*1000
    threshold_distance = threshold_mm
    adjacency_array[distances > threshold_distance] = 0

    adjacency_sparse = csr_matrix(adjacency)

    plot_adjacency = True

    if plot_adjacency is True:
        # plot the adjacency matrix
        plt.figure(figsize=(10, 10))
        plt.imshow(adjacency_array, cmap='gray', origin='upper')
        plt.title('Adjacency Matrix')
        plt.xlabel('Channels')
        plt.ylabel('Channels')
        plt.colorbar(label='Adjacency (1 = connected, 0 = not connected)')
        plt.show()

        # Plot the sensors and their thresholded adjacency as lines
        fig, ax = plt.subplots()
        mne.viz.plot_sensors(info, axes=ax, show=False)

        for ch_idx in range(adjacency_array.shape[0]):
            for neighbor_idx in range(adjacency_array.shape[1]):
                if adjacency_array[ch_idx, neighbor_idx]:
                    ax.plot([nirs_positions_adjacency[ch_idx, 0], nirs_positions_adjacency[neighbor_idx, 0]],
                            [nirs_positions_adjacency[ch_idx, 1], nirs_positions_adjacency[neighbor_idx, 1]], 'k-', alpha=0.5)

        plt.show()

    return adjacency_sparse
