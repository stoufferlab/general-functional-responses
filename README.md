# General Functional Responses
Currently, this repo contains all code (and most data) for reproducing the analyses associated with two manuscripts:

Stouffer & Novak (2020) *Hidden layers of density dependence in consumer feeding rates.* [https://doi.org/10.1101/2020.08.25.263806](https://doi.org/10.1101/2020.08.25.263806)

Novak & Stouffer (2020) *Systematic bias in studies of consumer functional responses.* bioRxiv. [https://doi.org/10.1101/2020.08.25.263814](https://doi.org/10.1101/2020.08.25.263814)

Please email Daniel Stouffer (daniel.stouffer@canterbury.ac.nz) or Mark Novak (mark.novak@oregonstate.edu) with any questions.

## Organization
The repo is mostly *not* split-up by manuscript because most code pertaining to *"Systematic bias..."* is also used in *"Hidden layers..."*.  Rather, it is primarily split-up by dataset type, distinguishing (1) functional-response datasets entailing variation in the abundance of both a (single) consumer species and its (single) resource species (*"One_Predator_One_Prey"*) from (2) functional-response datasets entailing variation in the abundance of two resources but no variation in the abundance of the consumer (*"One_Predator_Two_Prey"*).

### [code](code/)
##### [Mathematica](code/Mathematica)
Primarily for the calculation of Fisher Information Matrices, derivations of qualitative bias in functional-response parameter MLEs, and derivations of generalized functional-response models.
##### [R](code/R)
Contains a [library of scripts](code/R/lib) used in the model-fitting and subsequent analyses of both *One_Predator_One_Prey* and *One_Predator_Two_Prey* datasets.  Within each of the two latter folders, `fit_all_datasets.R` does as implied; otherwise use `RUNME.R` scripts.

### [data](data)
Contains all functional-response datasets (for which permission for public posting was granted) used in the analyses.  See [data/README.md](data/README.md) and the supplementary data tables in the supplementary materials of our manuscripts for sources and citations.

### [results](results)
##### [Mathematica](results/Mathematica)
Contains results used to generate Fig. 1 of *"Systematic bias..."*.
##### [R](results/R)
Contains all other results, figures and tables.

# Warranty
All code is provided "as is" and without warranty.  If you know of more efficient or elegant ways of doing anything (of which there are likely many), we’d love to learn from you.
