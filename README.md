Size exclusion chromatography (SEC) is a commonly used experimental techinque to estimate the MW of a protein.

However, elution times from a SEC colomn do not, in fact, depend on MW, but roughly upon radius of gyration.
For globular proteins, this distinction is irrelevant, but as the protein shape becomes non-globular, particularly if there are large internal voids, this matters greatly.

This code takes a list of protein elution times along with PDB structures, for both a standard ladder and the sample(s) of interest, and estimates the MW of the samples.
