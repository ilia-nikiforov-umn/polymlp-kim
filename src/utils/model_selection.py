"""Functions for model selection in converter of polynomial MLP repository."""

from collections import defaultdict
import numpy as np
from sklearn.cluster import SpectralClustering


def _get_representatives(time: np.ndarray, rmse_e: np.ndarray):
    """Return representatives using a clustering."""

    data = np.log(time).reshape((-1, 1))
    if data.shape[0] <= 5:
        return np.arange(data.shape[0])
        
    clustering = SpectralClustering(
        n_clusters=5,
        assign_labels='discretize',
        random_state=0,
    ).fit(data)
    groups = defaultdict(list)
    for i, lab in enumerate(clustering.labels_):
        groups[lab].append(i)

    reps = []
    threshold = 0.1
    for g in groups.values():
        if len(g) == 1:
            reps.append(g[0])
        else:
            rmse = rmse_e[np.array(g)]
            diff_rmse = rmse[1:] - rmse[:-1]
            if np.any(np.abs(diff_rmse) > threshold):
                t = np.log10(time[np.array(g)])
                diff_t = t[1:] - t[:-1]
                grads = diff_rmse / diff_t
                idx = np.argmin(grads) + 1
            else:
                idx = 0
            reps.append(g[idx])
    return np.array(reps)


def _screen_with_rmse(data: np.ndarray, e_thresholds: float = 5.0):
    """Eliminate MLPs with large RMSEs."""
    rmse_e = data[:, 2].astype(float)
    match = (rmse_e < e_thresholds)
    if np.count_nonzero(match) < 3:
        match = (rmse_e < e_thresholds * 2.0)
        if np.count_nonzero(match) == 0:
            raise RuntimeError("No optimal polymlp.")
    return data[match]


def choose_distributed_mlps(
    data: np.ndarray,
    e_thresholds: float = 5.0,
    time_scale: float = 1.0,
):
    """Find mlps on convex hull and select better MLPs."""
    data = _screen_with_rmse(data, e_thresholds=e_thresholds)
    times = data[:, 0].astype(float)
    rmse_e = data[:, 2].astype(float)
    rmse_f = data[:, 3].astype(float)
    score = rmse_e + rmse_f * 200
 
    reps = _get_representatives(times, score)
    data = data[reps]

    data[:, 0] = data[:, 0].astype(float) / time_scale
    data[:, 1] = data[:, 1].astype(float) / time_scale
    optimal_ids = data[:, 4]
    return optimal_ids, data
