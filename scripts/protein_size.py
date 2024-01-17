"""
A tool to compute the expected oligomer state of a protein based on it's structure and SEC elution times.

This is intended to be more accurate for non-globular proteins.
"""


from pathlib import Path

import typer
import pandas as pd

from sizer.sizer import compute_radius_of_gyration


def get_structure_path(pdb_name: str, pdb_top_dir: Path) -> Path:
    """
    Given the name of a PDB file and the PDB dir, guess the PDB Path.
    """
    top_dir = Path(pdb_top_dir)
    for file_path in top_dir.rglob(pdb_name):
        return file_path

    return None


def main(
    standards: Path = typer.Argument(
        ..., help=".csv of stanard elution times and structures"
    ),
    sample_elutions: Path = typer.Argument(..., help=".csv of sample elution times"),
    sample_structures: Path = typer.Argument(
        ..., help=".csv of sample structures and oligomer states"
    ),
    pdb_dir: Path = typer.Argument(..., help="Where all the PDB structures are held"),
    out_dir: Path = typer.Argument(..., help="Where to save output .csvs"),
):
    """
    Compute the expected oligomerization of a protein based on SEC elution time(s).
    """
    # load data
    standard_data = pd.read_csv(standards)
    sample_elution_data = pd.read_csv(sample_elutions)
    sample_structure_data = pd.read_csv(sample_structures)

    # compute radius of gyration
    standard_data["structure_path"] = standard_data["structure"].apply(
        get_structure_path, args=[pdb_dir]
    )
    standard_data["radius_gyration"] = standard_data["structure_path"].apply(
        compute_radius_of_gyration
    )
    standard_data = standard_data.drop(labels=["structure_path"], axis=1)

    sample_structure_data["structure_path"] = sample_structure_data["structure"].apply(
        get_structure_path, args=[pdb_dir]
    )
    sample_structure_data["radius_gyration"] = sample_structure_data[
        "structure_path"
    ].apply(compute_radius_of_gyration)
    sample_structure_data = sample_structure_data.drop(
        labels=["structure_path"], axis=1
    )

    # TODO: analyze
    # slope, intercept = np.polyfit(x, y, 1)
    # correlation_matrix = np.corrcoef(x, y)
    # correlation_xy = correlation_matrix[0,1]
    # r_squared = correlation_xy**2

    # TODO: report
    out_dir.mkdir(exist_ok=True)
    standard_data.to_csv(out_dir / "standards.csv")
    sample_structure_data.to_csv(out_dir / "structures.csv")
    sample_elution_data.to_csv(out_dir / "elutions.csv")

    print(standard_data)
    print(sample_elution_data)
    print(sample_structure_data)


if __name__ == "__main__":
    typer.run(main)
