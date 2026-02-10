"""Converter of polynomial MLP repository to KIM models."""

import numpy as np
import glob

from model_converter import convert_polymlp_to_kim_model
from model_selection import choose_distributed_mlps
from utils import parse_summary_yaml, get_polymlp_attrs


target = 1

if target == 1:
    path_summary = "/home/seko/mlip/1-unary/"
    path_repository = "/home/seko/mlip/repository/1-unary-2024/"
    repositories = sorted(glob.glob(path_repository + "*"))

    seq = 0
    for rep in repositories:
        repository_id = rep.split("/")[-1]
        system = repository_id.split("-")[0]
        summary_location = "/dataset-2024/5-opt/polymlp_summary_convex.yaml"
        summary = path_summary + system + summary_location
        data = parse_summary_yaml(summary)
        polymlp_ids, data = choose_distributed_mlps(data)
        mlps = [rep + "/polymlps/" + i + "/polymlp.lammps" for i in polymlp_ids]

        for i, (polymlp_file, d) in enumerate(zip(mlps, data)):
            print(polymlp_file, flush=True)
            time1, rmse_e, rmse_f = float(d[0]), float(d[2]), float(d[3])
            _, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            performance_level = i + 1

            if "manual" in polymlp_id:
                id1 = polymlp_id.replace("manual", "")
                id1 = "1" + id1.zfill(4)
            else:
                id1 = polymlp_id.replace("polymlp-", "").zfill(5)
            project_id = str(target) + str(seq).zfill(5) + id1

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
        mlps = [rep + "/pot/" + i + "/mlp.lammps" for i in polymlp_ids]

        for i, (polymlp_file, d) in enumerate(zip(mlps, data)):
            print(polymlp_file, flush=True)
            time1, rmse_e, rmse_f = float(d[0]), float(d[2]), float(d[3])
            _, polymlp_id, polymlp_year = get_polymlp_attrs(polymlp_file)
            performance_level = i + 1

            id1 = polymlp_id.split("-")[-1]
            id2 = 1 if "pair" in polymlp_id else 2
            id1 = str(id2) + id1.zfill(4)
            project_id = str(target) + str(seq).zfill(5) + id1

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
        performance_level = 1

        id1 = polymlp_id.split("-")[-1]
        id2 = 1 if "pair" in polymlp_id else 2
        project_id = str(target) + str(seq).zfill(5) + str(id2) + id1.zfill(4)

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
