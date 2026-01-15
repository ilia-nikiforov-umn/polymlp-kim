"""Converter of polynomial MLP repository to KIM models."""

import yaml
import glob
import numpy as np

from model_converter import convert_polymlp_to_kim_model
from scipy.spatial import ConvexHull


def parse_summary_yaml(summary_yaml: str):
    """Parse summary.yaml."""
    yamldata = yaml.safe_load(open(summary_yaml))
    data = [
        [d["cost_single"], d["cost_openmp"], d["rmse_energy"], d["rmse_force"], d["id"]]
        for d in yamldata["polymlps"]
    ]
    return np.array(data)


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
    return d_convex


def get_representatives_from_times(time: np.ndarray):
    """Return representatives using a clustering of computational times."""
    logtime = np.log2(time)
    logtime_int = logtime[1:] - logtime[:-1]
    groups = []
    acc = 0.0
    for i, interval in enumerate(logtime_int):
        if i == 0:
            ids = [0]
        acc += interval
        if acc <= 1.0:
            ids.append(i+1)
        else:
            groups.append(ids)
            acc = 0.0
            ids = [i + 1]
    if len(ids) > 0:
        groups.append(ids)
    reps = np.array([g[0] for g in groups])
    return reps 


def choose_distributed_mlps(
    data: np.ndarray,
    e_thresholds: float = 5.0,
    time_scale: float = 1.0,
):
    """Find mlps on convex hull."""
    score = data[:, 2].astype(float) + data[:, 3].astype(float) * 200
    data = np.hstack([data[:, 0:2], score.reshape((-1, 1)), data[:, 2:]])
    data = find_convex_mlps(data)

    rmse_e = data[:, 3].astype(float)
    match = (rmse_e < e_thresholds)
    if np.count_nonzero(match) < 3:
        match = (rmse_e < e_thresholds * 2.0)
        if np.count_nonzero(match) < 3:
            match = (rmse_e < e_thresholds * 3.0)
            if np.count_nonzero(match) == 0:
                raise RuntimeError("No optimal polymlp.")

    data = data[match]
    reps = get_representatives_from_times(data[:, 0].astype(float))
    data = data[reps]

    data[:, 0] = data[:, 0].astype(float) / time_scale
    data[:, 1] = data[:, 1].astype(float) / time_scale
    optimal_ids = data[:, 5]
    return optimal_ids, data[:, np.array([0, 1, 3, 4, 5])]


def get_polymlp_attrs(polymlp_file: str):
    """Return attributes of polymlp from directory name."""
    split = polymlp_file.split("/")
    repository_id = split[-4]
    polymlp_id = split[-2]
    polymlp_year = repository_id.split("-")[-3]
    return repository_id, polymlp_id, polymlp_year


if __name__ == "__main__":
    target = 3
    if target == 1:
        path_repository = "/home/seko/mlip/repository/1-unary-2024/"
        path_summary = "/home/seko/mlip/1-unary/"
        repositories = sorted(glob.glob(path_repository + "*"))

        seq = 0
        for rep in repositories:
            repository_id = rep.split("/")[-1]
            system = repository_id.split("-")[0]
            summary = path_summary + system \
                    + "/dataset-2024/5-opt/polymlp_summary_convex.yaml"
            data = parse_summary_yaml(summary)
            polymlp_ids, data = choose_distributed_mlps(data)
            mlps = [rep + "/polymlps/" + i + "/polymlp.lammps" for i in polymlp_ids]

            for i, (polymlp_file, d) in enumerate(zip(mlps, data)):
                print(polymlp_file, flush=True)
                time1, rmse_e, rmse_f = float(d[0]), float(d[2]), float(d[3])
                _, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
                id1 = polymlp_id.replace("polymlp-", "")
                if "manual" in id1:
                    id1 = id1.replace("manual", "")
                    id1 = "1" + id1.zfill(4)
                project_id = str(target) + str(seq).zfill(5) + id1 
                performance_level = i + 1

                convert_polymlp_to_kim_model(
                    polymlp_file, 
                    polymlp_id,
                    polymlp_year,
                    performance_level,
                    project_id,
                    repository_id,
                    time1, 
                    rmse_e, 
                    rmse_f,
                )
                seq += 1
    
    elif target == 2:
        path_repository = "/home/seko/mlip/repository/2-binary-alloy/"
        repositories = sorted(glob.glob(path_repository + "*"))

        seq = 0
        for rep in repositories:
            repository_id = rep.split("/")[-1]
            data = np.loadtxt(rep + "/summary/pareto", dtype=str)
            polymlp_ids, data = choose_distributed_mlps(data, time_scale=2.84)
            print(data)
            mlps = [rep + "/pot/" + i + "/mlp.lammps" for i in polymlp_ids]
            for i, (polymlp_file, d) in enumerate(zip(mlps, data)):
                print(polymlp_file, flush=True)
                time1, rmse_e, rmse_f = float(d[0]), float(d[2]), float(d[3])
                _, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
                id1 = polymlp_id.split("-")[-1]
                id2 = 1 if "pair" in polymlp_id else 2
                project_id = str(target) + str(i).zfill(5) + str(id2) + id1.zfill(4)
                performance_level = i + 1

                convert_polymlp_to_kim_model(
                    polymlp_file, 
                    polymlp_id,
                    polymlp_year,
                    performance_level,
                    project_id,
                    repository_id,
                    time1, 
                    rmse_e, 
                    rmse_f,
                )
                seq += 1

    elif target == 3:
        path_repository = "/home/seko/mlip/repository/3-ternary-alloy/"
        mlps = [path_repository + "Cu-Ag-Au-2024-05-09/pot/gtinv-182/polymlp.lammps"]
        seq = 0
        for i, polymlp_file in enumerate(mlps):
            print(polymlp_file, flush=True)
            repository_id, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            id1 = polymlp_id.split("-")[-1]
            id2 = 1 if "pair" in polymlp_id else 2
            project_id = str(target) + str(seq).zfill(5) + str(id2) + id1.zfill(4)
            performance_level = 1
            time1 = 0.880
            rmse_e = 1.44
            rmse_f = 0.0186
            convert_polymlp_to_kim_model(
               polymlp_file, 
               polymlp_id,
               polymlp_year,
               performance_level,
               project_id,
               repository_id,
               time1, 
               rmse_e, 
               rmse_f,
            )
            seq += 1
            

        # mlps = sorted(glob.glob(path_repository + "*/pot/*/*.lammps"))
        # for i, polymlp_file in enumerate(mlps):
        #     print(polymlp_file, flush=True)
        #     repository_id, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
        #     id1 = polymlp_id.split("-")[-1]
        #     id2 = 1 if "pair" in polymlp_id else 2
        #     project_id = str(target) + str(i).zfill(5) + str(id2) + id1.zfill(4)
        #     convert_polymlp_to_kim_model(
        #         polymlp_file, 
        #         repository_id,
        #         polymlp_id,
        #         polymlp_year,
        #         project_id,
        #     )
