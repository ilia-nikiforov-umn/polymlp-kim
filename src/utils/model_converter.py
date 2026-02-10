"""Converter of polynomial MLP to KIM model."""

from typing import Union

import sys
import os
import shutil

import numpy as np
import tarfile

from pypolymlp.core.io_polymlp import convert_to_yaml
from pypolymlp.utils.kim_utils import generate_kim_files, copy_mlps


def convert_polymlp_to_kim_model(
    polymlp_files: Union[str, list[str]],
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
    elements = copy_mlps(polymlp_files, path=tmp_path)

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

    content_origin = (
        "Polymlp Repository (Kyoto University)", 
        "https://cms.mtl.kyoto-u.ac.jp/seko/mlp-repository/index.html"
    )
    content_origin = " ".join(content_origin)

    user_id = "b3113743-4f85-48da-86e1-85acf6bb3388"
    contributor_id = user_id
    developer = [user_id]
    maintainer_id = user_id

    recordkey = ("MO", str(project_id).zfill(12), str(project_version).zfill(3))
    recordkey = "_".join(recordkey)
    citation1 = {
        "article-number": "{011101}",
        "author": "Seko, Atsuto",
        "doi": "10.1063/5.0129045",
        "journal": "{J. Appl. Phys.}",
        "eissn": "{1089-7550}",
        "issn": "{0021-8979}",
        "orcid-numbers": "{Seko, Atsuto/0000-0002-2473-3837}",
        "recordkey": recordkey + "a",
        "recordprimary": "recordprimary",
        "recordtype": "article",
        "unique-id": "{WOS:000908391700010}",
        "title": "{Tutorial: Systematic development of polynomial machine learning potentials for elemental and alloy systems}",
        "volume": "{133}",
        "number": "{1}",
        "year": "{2023}",
        "month": "{Jan}",
    }
    citations = [citation1]

    project = generate_kim_files(
        tmp_path,
        elements,
        polymlp_year=polymlp_year,
        performance_level=performance_level,
        project_id=project_id,
        project_version=project_version,
        description=description,
        model_driver=model_driver,
        content_origin=content_origin,
        contributor_id=contributor_id,
        developer=developer,
        maintainer_id=maintainer_id,
        citations=citations,
    )
    shutil.copy('LICENSE', tmp_path)
    shutil.move(tmp_path, project)
