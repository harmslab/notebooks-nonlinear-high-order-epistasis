# Detecting high-order epistasis in nonlinear genotype-phenotype maps

This repository contains all materials to reproduce the results and analysis in the following paper: "Detecting high-order epistasis in nonlinear genotype-phenotype maps". See the paper [here]().

## Try out the notebooks!

All results and analysis are reproducible in Jupyter notebooks. Try the the notebooks now, graciously provided by the Binder web service.

## Download and Install

If you'd like to run the notebooks locally, clone this repository and make sure all the necessary dependencies are installed. Here's a list of everything you'll need:

1. Jupyter
2. IPython
3. Numpy
4. Scipy
5. Scikit-learn
6. Matplotlib
7. ipywidgets
8. epistasis
9. seqspace

All packages can be installed using `pip`. 

## Data formats

The following formats are used in to create all data, metadata, and figures for the paper. The point of listing them is the provide documentation for anyone interested in sifting through through the data.

- Experimental data, and final results are saved in JSON format (`.json` files). This makes the data easily accessible to basically any programming language -- an attempt to practice true reproducibility. Also, this format is fairly human-readable.
- Bootstrapping results are saved as *.pickle* files, which are Python serialized pickle format. This contains sets of random samples drawn from normal distributions around the experimental phenotypes. It also includes epistatic coefficients estimated for each sample. 
- Most raw images are saved as `SVG`s. 
- Final images are exported as high-resolution `PNG`s.
- Data, analysis, raw-figures, code, etc. were all done in Jupyter Notebooks saved as `.ipynb` files. Guaranteed to work with Python 3.

## Types of Content

- **Images**: Plots, subplots, panels, etc. for all figures in the paper.
- **Notebooks**: Jupyter Notebooks with reproducible results for all figures in the paper.

## Table of Contents

- [Figure 1](): **Genotypic epistasis can be quantified using Walsh polynomials.** 
- [Figure 2](): **Nonlinearity in phenotype creates spurious high-order epistatic coefficients.**
- [Figure 3](): **Genotypic epistasis and nonlinear scale induce different patterns of nonadditivity.**
- [Figure 4](): **Experimental genotype-phenotype maps exhibit nonlinear phenotypes.** residuals.
- [Figure 5](): **High-order epistasis is present in genotype-phenotype maps.**
- [Figure 6](): **Nonlinear phenotypes distort measured epistatic coefficients.**