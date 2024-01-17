"""
Functions towards computing accurate oligomer states from SEC data of non-globular designs.
"""


from pathlib import Path

import biotite.structure as bts
from biotite.structure.io import load_structure


def compute_radius_of_gyration(structure_file: Path) -> float:
    """
    Return the radius of gyration for the structure in a PDB.
    """
    structure = load_structure(structure_file)

    return bts.gyration_radius(structure)
