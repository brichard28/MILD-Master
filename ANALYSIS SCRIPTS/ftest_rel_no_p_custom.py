def ftest_rel_no_p_custom(*args):
    """Perform one-sample t-test, converted into a F-stat for cluster testing.

    This is a modified version of :func:`scipy.stats.ttest_1samp` that avoids
    a (relatively) time-consuming p-value calculation, and can adjust
    for implausibly small variance values :footcite:`RidgwayEtAl2012`.

    Parameters
    ----------
    X : list of arrays
    sigma : float
        The variance estimate will be given by ``var + sigma * max(var)`` or
        ``var + sigma``, depending on "method". By default this is 0 (no
        adjustment). See Notes for details.
    method : str
        If 'relative', the minimum variance estimate will be sigma * max(var),
        if 'absolute' the minimum variance estimate will be sigma.

    Returns
    -------
    t : array
        T-values, potentially adjusted using the hat method.

    Notes
    -----
    To use the "hat" adjustment method :footcite:`RidgwayEtAl2012`, a value
    of ``sigma=1e-3`` may be a reasonable choice.

    References
    ----------
    footbibliography::
    """
    import numpy as np
    from functools import reduce

    X = [a for a in args]

    condition1 = X[0]  # Array of shape (n_subjects, n_channels) for condition 1
    condition2 = X[1]  # Array of shape (n_subjects, n_channels) for condition 2

    # Calculate the difference between the two conditions for each subject
    diff = condition1 - condition2

    # Calculate the variance of the differences (for each channel)
    var = np.nanvar(diff, axis=0, ddof=1)  # Using ddof=1 for unbiased variance
    sigma = 0  # no correction to start, but we can try 1e-3

    # Apply regularization if sigma > 0
    limit = sigma * np.max(var)
    var += limit

    # Compute the mean of the differences (for each channel)
    mean_diff = np.nanmean(diff, axis=0)
    
    # compute the standard error
    sem = np.sqrt(var) / np.sqrt(diff.shape[0])

    t_stat = mean_diff / sem
    # before, the above denonminator (sem) was np.sqrt(var / diff.shape[0]), which wouldn't quite be the standard error

    f_stat = np.square(t_stat)

    # Calculate t-statistics: mean difference divided by standard error
    return f_stat