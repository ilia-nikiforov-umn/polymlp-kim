"""Converter of polynomial MLP repository to KIM models."""

import glob
import numpy as np

from model_converter import convert_polymlp_to_kim_model


target = 2
if target == 1:
    path_repository = "/home/seko/mlip/repository/1-unary-2024/"
    mlps = sorted(glob.glob(path_repository + "*/polymlps/polymlp-*/polymlp.lammps*"))

    for i, polymlp_file in enumerate(mlps):
        print(polymlp_file, flush=True)
        split = polymlp_file.split("/")
        repository_id = split[-4]
        polymlp_id = split[-2]
        polymlp_year = repository_id.split("-")[-3]

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
    elif target == 3:
        path_repository = "/home/seko/mlip/repository/3-ternary-alloy/"

    mlps = sorted(glob.glob(path_repository + "*/pot/*/*.lammps"))

    for i, polymlp_file in enumerate(mlps):
        print(polymlp_file, flush=True)
        split = polymlp_file.split("/")
        repository_id = split[-4]
        polymlp_id = split[-2]
        polymlp_year = repository_id.split("-")[-3]

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
