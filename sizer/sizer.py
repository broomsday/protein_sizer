"""
Functions towards computing accurate oligomer states from SEC data of non-globular designs.
"""


from pathlib import Path

import pandas as pd

import biotite.structure as bts
from biotite.structure.io import load_structure


def compute_radius_of_gyration(structure_file: Path) -> float:
    """
    Return the radius of gyration for the structure in a PDB.
    """
    structure = load_structure(structure_file)

    return bts.gyration_radius(structure)


def compute_partition_coefficient(elution: float, void: float, total: float) -> float:
    """
    Compute the SEC partition coefficient.

    K = (Ve - Vo) / (Vt - Vo)
    Ve = elution volume (peak max height)
    Vt = total volume, computed as pi*h*r^2 (cylindrical volume) of the packed bed
    Vo = void volume, elution volume of something so large as to be completely excluded

    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8410290/
    """
    return (elution - void) / (total - void)


def estimate_radius_of_gyration(
    partition: float, slope: float, intercept: float
) -> float:
    """
    Given linear fit parameters and the partition coefficient, estimate the radius of gyration.
    """
    log_10_rg = (slope * partition) + intercept
    return 10**log_10_rg


def get_closest_oligomer(radius_gyration: float, structure_data: pd.DataFrame) -> int:
    """
    Report the oligomer structure that most closely matches the radius of gyration.
    """
    structure_data["diff"] = abs(structure_data["radius_gyration"] - radius_gyration)
    structure_data = structure_data.sort_values(by="diff")

    return structure_data.iloc[0].oligomer
