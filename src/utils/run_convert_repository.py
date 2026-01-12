"""Converter of polynomial MLP repository to KIM models."""

import yaml
import glob
import numpy as np

from model_converter import convert_polymlp_to_kim_model
from scipy.spatial import ConvexHull


def find_convex_mlps(
    d_array: np.ndarray,
    use_force: bool = False,
    use_logscale_time: bool = False,
):
    """Find convex polymlps from data."""
    d_target = d_array[:, [0, 3]] if use_force else d_array[:, [0, 2]]
    if use_logscale_time:
        time = np.log10(d_target[:, 0].astype(float))
        error = d_target[:, 1].astype(float)
        d_target = np.vstack([time, error]).T

    ch = ConvexHull(d_target)
    v_convex = np.unique(ch.simplices)

    data_values = d_array[v_convex, :4].astype(float)
    data_str = d_array[v_convex, 4:]

    v_convex_l = []
    for i1, d1 in enumerate(data_values):
        lower_convex = True
        time1, rmse_e1, rmse_f1 = d1[0], d1[2], d1[3]
        for i2, d2 in enumerate(data_values):
            if i1 == i2:
                continue
            time2, rmse_e2, rmse_f2 = d2[0], d2[2], d2[3]
            if time2 < time1 and rmse_e2 < rmse_e1 and rmse_f2 < rmse_f1:
                lower_convex = False
                break
        if lower_convex:
            v_convex_l.append(i1)
    v_convex_l = np.array(v_convex_l)
    d_convex = np.hstack([data_values[v_convex_l], data_str[v_convex_l]])
    d_convex = d_convex[d_convex[:, 2].astype(float) < 30]
    return d_convex


def eval_score(
    time: np.ndarray,
    rmse_e: np.ndarray,
    rmse_f: np.ndarray,
    weight_t: float = 2.0,
):
    """Define score."""
    return weight_t * (np.log10(time) + 6) + (rmse_e + rmse_f * 300)


def parse_summary_yaml(summary_yaml: str):
    """Parse summary.yaml."""
    yamldata = yaml.safe_load(open(summary_yaml))
    data = [
        [d["cost_single"], d["cost_openmp"], d["rmse_energy"], d["rmse_force"], d["id"]]
        for d in yamldata["polymlps"]
    ]
    return np.array(data)


def choose_distributed_mlps(
    data: np.ndarray,
    e_thresholds: float = 5.0,
    weights_time: tuple = (1, 10, 100),
    n_atom: int = 1,
):
    """Find mlps on convex hull."""
    data = find_convex_mlps(data)

    rmse_e = data[:, 2].astype(float)
    match = (rmse_e < e_thresholds)
    if np.count_nonzero(match) < 2:
        match = (rmse_e < e_thresholds * 2.0)
        if np.count_nonzero(match) < 2:
            match = (rmse_e < e_thresholds * 3.0)
            if np.count_nonzero(match) == 0:
                raise RuntimeError("No optimal polymlp.")

    data = data[match]
    time = data[:, 0].astype(float) / n_atom
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
    target = 1
    if target == 1:
        path_repository = "/home/seko/mlip/repository/1-unary-2024/"
        path_summary = "/home/seko/mlip/1-unary/"
        repositories = sorted(glob.glob(path_repository + "*"))
        mlps = []
        for rep in repositories:
            system = rep.split("/")[-1].split("-")[0]
            summary = path_summary + system \
                    + "/dataset-2024/5-opt/polymlp_summary_convex.yaml"
            data = parse_summary_yaml(summary)
            for polymlp_id in choose_distributed_mlps(
                data, weights_time=(1, 3, 4, 5, 10),
            ):
                pot = glob.glob(rep + "/polymlps/" + polymlp_id + "/polymlp.lammps*")[0]
                mlps.append(pot)
 
        for i, polymlp_file in enumerate(mlps):
            print(polymlp_file, flush=True)
            repository_id, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            id1 = polymlp_id.replace("polymlp-", "")
            if "manual" in id1:
                id1 = id1.replace("manual", "")
                id1 = "1" + id1.zfill(4)
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
                data = np.loadtxt(rep + "/summary/pareto", dtype=str)
                for polymlp_id in choose_distributed_mlps(data, n_atom=2840):
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
