# Tailoring mathematical models to stem-cell derived cardiomyocyte lines can improve predictions of drug-induced changes to their electrophysiology

This page contains the code accompanying the article
_Tailoring mathematical models to stem-cell derived cardiomyocyte lines can_
_improve predictions of drug-induced changes to their electrophysiology_ by
Chon Lok Lei, Ken Wang, Michael Clerx, Ross H. Johnstone, Maria P.
Hortigon-Vinagre, Victor Zamora, Andrew Allan, Godfrey L. Smith, David J.
Gavaghan, Gary R. Mirams and Liudmila Polonchuk.

See: https://doi.org/10.3389/fphys.2017.00986

## Requirements

Python 2.7, including the modules `numpy`, `scipy`, `matplotlib` and `cma`
 (which can be installed using `pip install cma`).

CVODE 2.6 or higher (see https://computation.llnl.gov/projects/sundials or
 http://myokit.org/download).

## How to use

To run all code, and recreate all figures, simply run the bash script
 `run-all`. On windows systems, the Python files called from this script can be
 invoked manually, e.g. `python figure1a.py`.
