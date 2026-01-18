"""Functions for model selection in converter of polynomial MLP repository."""

from collections import defaultdict
import numpy as np
from sklearn.cluster import KMeans


def _get_representatives(time: np.ndarray, score: np.ndarray, n_clusters: int = 5):
    """Return representatives using a clustering."""

    data = score.reshape((-1, 1))
    if data.shape[0] <= n_clusters:
        return np.arange(data.shape[0])

    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    labels = kmeans.fit_predict(data)
    groups = defaultdict(list)
    for i, lab in enumerate(labels):
        groups[lab].append(i)

    reps = []
    threshold = 0.5
    for g in groups.values():
        if len(g) == 1:
            reps.append(g[0])
        else:
            labels = np.array(g)
            y = score[labels]
            t = np.log10(time[labels])
            dy = y[-1] - y[0]
            dt = t[-1] - t[0]
            if (dy / dt) > -threshold:
                reps.append(g[0])
            else:
                dy = y[1:] - y[:-1]
                dt = t[1:] - t[:-1]
                grads = dy / dt
                idx = np.argmin(grads) + 1
                reps.append(g[idx])
    reps = np.sort(reps)
    return reps


def _screen_with_rmse(data: np.ndarray, e_thresholds: float = 5.0):
    """Eliminate MLPs with large RMSEs."""
    rmse_e = data[:, 2].astype(float)
    match = (rmse_e < e_thresholds)
    if np.count_nonzero(match) < 3:
        match = (rmse_e < e_thresholds * 1.6)
        if np.count_nonzero(match) == 0:
            match = (rmse_e < e_thresholds * 2.0)
            if np.count_nonzero(match) == 0:
                raise RuntimeError("No optimal polymlp.")
    return data[match]


def choose_distributed_mlps(
    data: np.ndarray,
    e_thresholds: float = 5.0,
    time_scale: float = 1.0,
    n_clusters: int = 5,
):
    """Find mlps on convex hull and select better MLPs."""
    data = _screen_with_rmse(data, e_thresholds=e_thresholds)
    times = data[:, 0].astype(float)
    rmse_e = data[:, 2].astype(float)
    rmse_f = data[:, 3].astype(float)
    score = rmse_e + rmse_f * 200
 
    reps = _get_representatives(times, score, n_clusters=n_clusters)
    data = data[reps]

    data[:, 0] = data[:, 0].astype(float) / time_scale
    data[:, 1] = data[:, 1].astype(float) / time_scale
    optimal_ids = data[:, 4]
    return optimal_ids, data
