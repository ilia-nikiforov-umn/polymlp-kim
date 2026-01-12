"""Converter of polynomial MLP repository to KIM models."""

import glob
import numpy as np

from model_converter import convert_polymlp_to_kim_model

def eval_score(
    time: np.ndarray,
    rmse_e: np.ndarray,
    rmse_f: np.ndarray,
    weight_t: float = 2.0,
):
    """Define score."""
    return weight_t * (np.log10(time) + 6) + (rmse_e + rmse_f * 300)


def choose_distributed_mlps(
    summary_file: str,
    e_thresholds: float = 5.0,
    weights_time: tuple = (1, 10, 100),
):
    """Find mlps on convex hull."""
    data = np.loadtxt(summary_file, dtype=str)
    rmse_e = data[:, 2].astype(float)
    match = (rmse_e < e_thresholds)
    if np.count_nonzero(match) == 0:
        match = (rmse_e < e_thresholds * 2.0)
        if np.count_nonzero(match) == 0:
            match = (rmse_e < e_thresholds * 3.0)
            if np.count_nonzero(match) == 0:
                raise RuntimeError("No optimal polymlp.")

    data = data[match]
    time = data[:, 0].astype(float) / 2840
    rmse_e = data[:, 2].astype(float)
    rmse_f = data[:, 3].astype(float)
    ids = data[:, 4]

    optimal_ids = []
    for weight_t in weights_time:
        scores = eval_score(time, rmse_e, rmse_f, weight_t=weight_t)
        idmin = np.argmin(scores)
        optimal_ids.append(ids[idmin])
    return np.unique(optimal_ids)


def get_polymlp_attrs(polymlp_file: str):
    """Return attributes of polymlp from directory name."""
    split = polymlp_file.split("/")
    repository_id = split[-4]
    polymlp_id = split[-2]
    polymlp_year = repository_id.split("-")[-3]
    return repository_id, polymlp_id, polymlp_year


if __name__ == "__main__":
    target = 2
    if target == 1:
        path_repository = "/home/seko/mlip/repository/1-unary-2024/"
        mlps = sorted(glob.glob(path_repository + "*/polymlps/polymlp-*/polymlp.lammps*"))
    
        for i, polymlp_file in enumerate(mlps):
            print(polymlp_file, flush=True)
            repository_id, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            id1 = polymlp_id.replace("polymlp-", "")
            project_id = str(target) + str(i).zfill(5) + id1 
            convert_polymlp_to_kim_model(
                polymlp_file, 
                repository_id,
                polymlp_id,
                polymlp_year,
                project_id,
            )
    
    elif target == 2 or target == 3:
        if target == 2:
            path_repository = "/home/seko/mlip/repository/2-binary-alloy/"
            repositories = sorted(glob.glob(path_repository + "*"))
            mlps = []
            for rep in repositories:
                for polymlp_id in choose_distributed_mlps(rep + "/summary/pareto"):
                    pot = rep + "/pot/" + polymlp_id + "/mlp.lammps"
                    mlps.append(pot)
        elif target == 3:
            path_repository = "/home/seko/mlip/repository/3-ternary-alloy/"
            mlps = sorted(glob.glob(path_repository + "*/pot/*/*.lammps"))
    
        for i, polymlp_file in enumerate(mlps):
            print(polymlp_file, flush=True)
            repository_id, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            id1 = polymlp_id.split("-")[-1]
            id2 = 1 if "pair" in polymlp_id else 2
            project_id = str(target) + str(i).zfill(5) + str(id2) + id1.zfill(4)
            convert_polymlp_to_kim_model(
                polymlp_file, 
                repository_id,
                polymlp_id,
                polymlp_year,
                project_id,
            )
