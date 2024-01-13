"""
A tool to compute the expected oligomer state of a protein based on it's structure and SEC elution times.

This is intended to be more accurate for non-globular proteins.
"""


from pathlib import Path

import typer
import pandas as pd


def main(
    standards: Path = typer.Argument(..., help=".csv of stanard elution times and structures"),
    sample_elutions: Path = typer.Argument(..., help=".csv of sample elution times"),
    sample_structures: Path = typer.Argument(..., help=".csv of sample structures and oligomer states"),
):
    """
    Compute the expected oligomerization of a protein based on SEC elution time(s).
    """
    standards_data = pd.read_csv(standards)
    sample_elution_data = pd.read_csv(sample_elutions)
    sample_structure_data = pd.read_csv(sample_structures)

    print(standards_data)
    print(sample_elution_data)
    print(sample_structure_data)



if __name__ == "__main__":
    typer.run(main)