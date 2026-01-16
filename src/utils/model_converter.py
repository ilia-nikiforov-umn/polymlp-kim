"""Converter of polynomial MLP to KIM model."""

import sys
import os
import shutil

import numpy as np
import tarfile

from pypolymlp.core.io_polymlp import convert_to_yaml, load_mlp
from pypolymlp.utils.kim_utils import generate_kim_files


def convert_polymlp_to_kim_model(
    polymlp_file: str,
    polymlp_id: str,
    polymlp_year: int,
    performance_level: int,
    project_id: int,
    repository_id: str,
    time_single_core: float,
    rmse_energy: float,
    rmse_force: float,
    project_version: int = 0,
    model_driver: str = "Polymlp__MD_000000123456_000",
):
    """Convert polymlp to KIM model."""

    tmp_path = "./Polymlp__MO_tmp/"
    os.makedirs(tmp_path, exist_ok=True)

    polymlp_yaml = tmp_path + "polymlp.yaml"
    if ".lammps.tar.gz" in polymlp_file:
        tarpath = "/".join(polymlp_file.split("/")[:-1])
        with tarfile.open(polymlp_file) as tar:
            tar.extractall(path=tarpath)
        convert_to_yaml(polymlp_file.replace(".tar.gz",""), yaml=polymlp_yaml)
    elif ".lammps" in polymlp_file:
        convert_to_yaml(polymlp_file, yaml=polymlp_yaml)
    else:
        shutil.copy(polymlp_file, polymlp_yaml)

    params, _ = load_mlp(polymlp_yaml)
    elements = params.elements

    description = (
        "Polynomial machine learning potential (MLP) for", 
        "-".join(elements),
        "system. This potential is",
        polymlp_id, 
        "potential in",
        repository_id,
        "taken from Polynomial MLP Repository.",
        "The RMS errors for energy and forces are",
        str(np.round(rmse_energy, 3)),
        "meV/atom and",
        str(np.round(rmse_force, 4)),
        "eV/angstrom, respectively.",
        "The estimated computational cost required for a single core calculation is",
        str(np.round(time_single_core, 2)),
        "ms/atom/step.",
    )
    description = " ".join(description)

    project = generate_kim_files(
        tmp_path,
        elements,
        polymlp_year=polymlp_year,
        performance_level=performance_level,
        project_id=project_id,
        project_version=project_version,
        description=description,
        model_driver=model_driver,
    )
    shutil.move(tmp_path, project)
