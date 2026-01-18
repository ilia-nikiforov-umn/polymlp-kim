"""Converter of polynomial MLP repository to KIM models."""

import numpy as np
import glob

from model_converter import convert_polymlp_to_kim_model
from model_selection import choose_distributed_mlps
from utils import parse_summary_yaml, get_polymlp_attrs


target = 1

if target == 1:
    path_summary = "/home/seko/mlip/1-unary/"
    path_repository = "/home/seko/mlip/repository/1-unary-hybrid-2024/"
    repositories = sorted(glob.glob(path_repository + "*"))

    seq = 0
    for rep in repositories:
        repository_id = rep.split("/")[-1]
        system = repository_id.split("-")[0]
        summary_location = "/dataset-2024/7-opt-hybrid/polymlp_summary_convex.yaml"
        summary = path_summary + system + summary_location
        data = parse_summary_yaml(summary)
        polymlp_ids, data = choose_distributed_mlps(data, n_clusters=4)
        is_hybrid = np.char.find(polymlp_ids, "hybrid") > 0
        if np.count_nonzero(is_hybrid) > 2:
            polymlp_ids, data = choose_distributed_mlps(data, n_clusters=3)
            is_hybrid = np.char.find(polymlp_ids, "hybrid") > 0
            if np.count_nonzero(is_hybrid) > 2:
                polymlp_ids, data = choose_distributed_mlps(data, n_clusters=2)
                is_hybrid = np.char.find(polymlp_ids, "hybrid") > 0
        polymlp_ids = polymlp_ids[is_hybrid]
        data = data[is_hybrid]

        mlps = []
        for i in polymlp_ids:
            files = glob.glob(rep + "/polymlps/" + i + "/polymlp.lammps*")
            mlps.append(files)

        for i, (polymlp_files, d) in enumerate(zip(mlps, data)):
            print(polymlp_files, flush=True)
            time1, rmse_e, rmse_f = float(d[0]), float(d[2]), float(d[3])
            _, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_files[0])
            performance_level = i + 1
            print(polymlp_id, polymlp_year)

            id1 = polymlp_id.replace("polymlp-", "").replace("hybrid-", "")
            id1 = id1.split("-")
            id1 = "".join([id1[0][-3:], id1[1][-2:]]).zfill(5)
            project_id = "1" + str(target) + str(seq).zfill(5) + id1

            convert_polymlp_to_kim_model(
                polymlp_files,
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
