"""Compute the cohesive energy, forces, and pressure of given POSCAR."""

from pypolymlp.api.pypolymlp_utils import PypolymlpUtils

polymlp = PypolymlpUtils()
polymlp.generate_kim_model(
    "polymlp.yaml",
    performance_level=1,
    project_id=1234567,
    project_version=1,
    model_driver="Polymlp__MD_000000123456_000",
)
