def mne_cluster_based_permutation_test_cpu(adjacency_sparse, X):
    """
    Perform cluster-based permutation test using MNE, running only on CPU.

    :param adjacency_sparse: Sparse adjacency matrix for clustering
    :param X: Data array [(n_samples, n_tests),...]
    :return: t-values, clusters, cluster p-values, H0 distribution, significant channels
    """
    import mne
    import numpy as np
    from ftest_rel_no_p_custom import ftest_rel_no_p_custom

    # Use None for threshold to allow permutation-based significance testing
    # Perform the permutation cluster test
    F, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_test(
        X,
        n_permutations=1000,
        adjacency=adjacency_sparse,
        stat_fun= ftest_rel_no_p_custom,
        n_jobs=1  # Ensures CPU usage
    )

    # Identify significant clusters
    sig_clusters = np.where(cluster_p_values < 0.05)[0]
    clusters_temp = np.array([cluster[0] for cluster in clusters], dtype=object)

    is_valid = bool(clusters_temp[sig_clusters].size) and any(
        subarray.size for subarray in clusters_temp[sig_clusters]
    )

    significant_channels = (
        np.concatenate(clusters_temp[sig_clusters]) if is_valid else None
    )

    return F, clusters, cluster_p_values, H0, significant_channels
