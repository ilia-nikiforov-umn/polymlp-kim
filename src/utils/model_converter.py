"""Converter of polynomial MLP to KIM model."""

import sys
import os
import shutil

import numpy as np

from pypolymlp.core.io_polymlp import convert_to_yaml, load_mlp


def convert_polymlp_to_kim_model(
    polymlp_file: str,
    repository_id: str,
    polymlp_id: str,
    polymlp_year: int,
    project_id: int,
    project_version: int = 0,
    model_driver: str = "Polymlp__MD_000000123456_000",
):
    """Convert polymlp to KIM model."""

    tmp_path = "./Polymlp__MO_tmp/"
    os.makedirs(tmp_path, exist_ok=True)

    polymlp_yaml = tmp_path + "polymlp.yaml"
    if ".lammps" in polymlp_file:
        convert_to_yaml(polymlp_file, yaml=polymlp_yaml)
    else:
        shutil.copy(polymlp_file, polymlp_yaml)

    params, _ = load_mlp(polymlp_yaml)
    elements = params.elements
    system = "-".join(elements)
    project = "Polymlp_Seko_" + str(polymlp_year) + "_" \
            + "".join(elements) + "__MO_" \
            + str(project_id).zfill(12) + "_" + str(project_version).zfill(3)

    with open(tmp_path + "polymlp.settings", "w") as f:
        print(" ".join(elements), file=f)
        print(file=f)

    np.set_printoptions(formatter={'str_kind': lambda x: f'"{x}"'})
    with open(tmp_path + "CMakeLists.txt", "w") as f:
        print("#", file=f)
        print("# Author: Atsuto Seko", file=f)
        print("#", file=f)
        print(file=f)
        print("cmake_minimum_required(VERSION 3.10)", file=f)
        print(file=f)
        print("list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})", file=f)
        print("find_package(KIM-API-ITEMS 2.2 REQUIRED CONFIG)", file=f)
        print(file=f)
        print("kim_api_items_setup_before_project(ITEM_TYPE \"portableModel\")", file=f)
        print("project(" + project + ")", file=f)
        print("kim_api_items_setup_after_project(ITEM_TYPE \"portableModel\")", file=f)
        print(file=f)
        print("add_kim_api_model_library(", file=f)
        print("  NAME            ${PROJECT_NAME}", file=f)
        print("  DRIVER_NAME     \"" + model_driver + "\"", file=f)
        print("  PARAMETER_FILES \"polymlp.settings\"", file=f)
        print("                  \"polymlp.yaml\"", file=f)
        print("  )", file=f)

    with open(tmp_path + "kimspec.edn", "w") as f:
        print("{\"domain\" \"openkim.org\"", file=f)
        print(" \"extended-id\" \"" + project + "\"", file=f)
        print(" \"kim-api-version\" \"2.2\"", file=f)
        print(" \"model-driver\" \"" + model_driver + "\"", file=f)
        print(" \"potential-type\" \"polymlp\"", file=f)
        print(" \"publication-year\" \"" + str(polymlp_year) + "\"", file=f)
        print(" \"species\"",  np.array(elements), file=f)
        print(
            " \"title\" \"Polynomial machine learning potential for",
            system, 
            "developed by Seko (" + str(polymlp_year) + ")",
            "v" + str(project_version).zfill(3) + "\"", 
            file=f,
        )
        print(" }", file=f)

    shutil.move(tmp_path, project)

    # TODO: change directory name.


if __name__ == "__main__":
    polymlp_file = sys.argv[1]
    split = polymlp_file.split("/")
    repository_id = split[-4]
    polymlp_id = split[-2]
    polymlp_year = repository_id.split("-")[-3]

    project_id = 325

    convert_polymlp_to_kim_model(
        polymlp_file, 
        repository_id,
        polymlp_id,
        polymlp_year,
        project_id,
    )
