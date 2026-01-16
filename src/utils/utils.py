"""Utility functions for converter of polynomial MLP repository."""

import yaml
import numpy as np


def parse_summary_yaml(summary_yaml: str):
    """Parse summary.yaml."""
    yamldata = yaml.safe_load(open(summary_yaml))
    data = [
        [d["cost_single"], d["cost_openmp"], d["rmse_energy"], d["rmse_force"], d["id"]]
        for d in yamldata["polymlps"]
    ]
    return np.array(data)


def get_polymlp_attrs(polymlp_file: str):
    """Return attributes of polymlp from directory name.

    When polymlp_file is path/Ag-2024-06-05/polymlps/polymlp-00135/polymlp.lammps,
    three quantities of
        repository_id = "Ag-2024-06-05"
        polymlp_id = "polymlp-00135"
        polymlp_year = "2024"
    will be returned.
    """
    split = polymlp_file.split("/")
    repository_id = split[-4]
    polymlp_id = split[-2]
    polymlp_year = repository_id.split("-")[-3]
    return repository_id, polymlp_id, polymlp_year
