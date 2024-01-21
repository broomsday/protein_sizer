"""
A tool to compute the expected oligomer state of a protein based on it's structure and SEC elution times.

This is intended to be more accurate for non-globular proteins.
"""


from pathlib import Path

import typer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sizer.sizer import (
    compute_radius_of_gyration,
    compute_partition_coefficient,
    estimate_radius_of_gyration,
    get_closest_oligomer,
)


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
    column_parameters: Path = typer.Argument(
        ..., help=".csv giving the total and void volumes of the column"
    ),
    pdb_dir: Path = typer.Argument(..., help="Where all the PDB structures are held"),
    out_dir: Path = typer.Argument(..., help="Where to save output .csvs"),
    debug: bool = typer.Option(default=False, help="Show standard curve plot"),
):
    """
    Compute the expected oligomerization of a protein based on SEC elution time(s).

    Example:
    python scripts/protein_size.py data/nanopore_2XQS/standards.csv data/nanopore_2XQS/elutions.csv data/nanopore_2XQS/structures.csv data/nanopore_2XQS/pdbs/ data/nanopore_2XQS/output
    """
    # load data
    standard_data = pd.read_csv(standards)
    sample_elution_data = pd.read_csv(sample_elutions)
    sample_structure_data = pd.read_csv(sample_structures)
    column_parameter_data = pd.read_csv(column_parameters)

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

    # pull out needed column parameters
    void_volume = (
        column_parameter_data[column_parameter_data["parameter"] == "void"]
        .iloc[0]
        .elution
    )
    total_volume = (
        column_parameter_data[column_parameter_data["parameter"] == "total"]
        .iloc[0]
        .elution
    )

    # compute the partition coefficient for each elution volume
    standard_data["partition"] = standard_data["elution"].apply(
        compute_partition_coefficient, args=(void_volume, total_volume)
    )
    sample_elution_data["partition"] = sample_elution_data["elution"].apply(
        compute_partition_coefficient, args=(void_volume, total_volume)
    )

    # save the standards and structures rg calculations
    out_dir.mkdir(exist_ok=True)
    standard_data.to_csv(out_dir / "standards.csv")
    sample_structure_data.to_csv(out_dir / "structures.csv")

    # a plot of K vs. Log10(size) (e.g radius of gyration) is expected to be linear
    if debug:
        plt.scatter(
            standard_data["partition"], np.log10(standard_data["radius_gyration"])
        )
        plt.show()
    slope, intercept = np.polyfit(
        standard_data["partition"], np.log10(standard_data["radius_gyration"]), 1
    )
    correlation_matrix = np.corrcoef(
        standard_data["partition"], np.log10(standard_data["radius_gyration"])
    )
    correlation_xy = correlation_matrix[0, 1]
    r_squared = correlation_xy**2

    print(f"Standard curve R^2: {r_squared}\n")

    # use the slope and intercept to compute the radius_gyration for the samples
    sample_elution_data["radius_gyration"] = sample_elution_data["partition"].apply(
        estimate_radius_of_gyration, args=(slope, intercept)
    )

    # use the radius_gyration of the structures to estimate the closest matching oligomer
    sample_elution_data["closest_oligomer"] = sample_elution_data[
        "radius_gyration"
    ].apply(get_closest_oligomer, args=(sample_structure_data,))

    # save the estimated sample data
    sample_elution_data.to_csv(out_dir / "elutions.csv")

    # report the data
    print(standard_data, "\n")
    print(sample_elution_data, "\n")
    print(sample_structure_data, "\n")


if __name__ == "__main__":
    typer.run(main)
