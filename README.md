# Detecting high-order epistasis in nonlinear genotype-phenotype maps

This repository contains all materials to reproduce all results and analysis in the following paper: "Detecting high-order epistasis in nonlinear genotype-phenotype maps". See the paper [here]().

# Try out the notebooks!

All results and analysis are reproducible in Jupyter notebooks. Try the the notebooks now, graciously provided by the Binder web service.

# Data formats

The following formats are used in to create all data, metadata, and figures for the paper. The point of listing them is the provide documentation for anyone interested in sifting through through the data.

- Experimental data, and final results are saved in JSON format (`.json` files). This makes the data easily accessible to basically any programming language -- an attempt to practice true reproducibility. Also, this format is fairly human-readable.
- Bootstrapping results are saved as *.pickle* files, which are Python serialized pickle format. This contains sets of random samples drawn from normal distributions around the experimental phenotypes. It also includes epistatic coefficients estimated for each sample. 
- Most raw images are saved as `SVG`s. 
- Final images are exported as high-resolution `PNG`s.
- Data, analysis, raw-figures, code, etc. were all done in Jupyter Notebooks saved as `.ipynb` files. Guaranteed to work with Python 3.

# Types of Content

- **Images**: Plots, subplots, panels, etc. for all figures in the paper.
- **Notebooks**: Jupyter Notebooks with reproducible results for all figures in the paper.
- 

# Table of Contents

- [Figure 1](): **Genotypic epistasis can be quantified using Walsh polynomials.** A) A genotype-phenotype map exhibiting negative epistasis. Axes are genotype at position 1 (g1), genotype at position 2 (g2), and phenotype (P). For genotypic axes, “0” denotes wildtype and “1” denotes a mutant. Phenotype is encoded both on the P-axis and as a spectrum from white to blue. The map exhibits negative epistasis: relative to wildtype, the effect of the mutations together (P11 = 2) is less than the sum of the individual effects of mutations (P10 + P01 = 1 + 2 = 3). B) The map can be decomposed into epistatic coefficients using a Walsh polynomial, which measures the effects of each mutation relative to the geometric center of the genotype-phenotype map (green sphere). The additive coefficients β1 and β2 (red arrows) are the average effect of each mutation in all backgrounds. The epistatic coefficient β12 (orange arrow) the variation not accounted for by β1 and β2. Geometrically, it is the distance between the center of the map and the “fold” given by vector connecting P00 and P11.

- [Figure 2](): **Nonlinearity in phenotype creates spurious high-order epistatic coefficients.** A) Simulated, random, first-order epistatic coefficients. The mutated site is indicated by panel below the bar graph; bar indicates magnitude and sign of the epistatic coefficient. B) A nonlinear map between a linear phenotype and a saturating, nonlinear phenotype. The first-order coefficients in panel A are used to generate a linear phenotype, which is then transformed by the function shown in B. C) Epistatic coefficients extracted from the genotype-phenotype map generated in panels A and B. Bars denote coefficient magnitude and sign. Color denotes the order of the coefficient: first (βi, red), second (βij, orange), third (βijk, green), fourth (βijkl, purple), and fifth (βijklm, blue). Filled squares in the grid below the bars indicate the identity of mutations that contribute to the coefficient.

- [Figure 3](): **Genotypic epistasis and nonlinear scale induce different patterns of nonadditivity.** A) Patterns of nonadditivity for increasing genotypic epistasis and nonlinear scale. Main panel shows grid ranging from no epistasis (bottom left) to high genotypic epistasis and nonlinearity (top right). Insets in sub-panels show added nonlinearity. Going from left to right: K = 0, K = 2, 421 K = 4. Epistatic coefficient plots to right show magnitude of input genotypic epistasis, with colors and annotation as in Fig 2C. B) Plot of Pobs against Pˆadd for the middle sub panel in panel A. Red line is the fit of the power transform to these data. C) Correlation between epistatic coefficients input into the simulation and extracted from the simulation after correction for nonlinearity with the power transform. Each point is an epistatic coefficient, colored by order. The Pearson’s correlation coefficient is shown in the upper-left quadrant. D) Correlation between epistatic coefficients input into the simulation and extracted from the simulation without application of the power transform.

- [Figure 4](): **Experimental genotype-phenotype maps exhibit nonlinear phenotypes.** Plots show observed phenotype Pobs plotted against Pˆadd (Eq. 1) for data sets I through IV. Points are individual genotypes. Error bars are experimental standard deviations in phenotype. Red lines are the fit of the power transform to the data set. Pearson’s coefficient for each fit are shown on each plot. Dashed lines are Padd = Pobs. Bottom panels in each plot show residuals between the observed phenotypes and the red fit line. Points are the individual residuals. Errorbars are the experimental standard deviation of the phenotype. The horizontal histograms show the distribution of residuals across 10 bins. The red lines are the mean of the residuals.

- [Figure 5](): **High-order epistasis is present in genotype-phenotype maps.** A) Panels show epistatic coefficients extracted from data sets I-IV (Table 1, data set label circled above each graph). Bars denote coefficient magnitude and sign; error bars are propagated measurement uncertainty. Color denotes the order of the coefficient: first (βi, red), second (βij, orange), third (βijk, green), fourth (βijkl, purple), and fifth (βijklm, blue). Bars are colored if the coefficient is significantly different than zero (Z-score with p-value < 0.05 after Bonferroni correction for multiple testing). Stars denote relative significance: p < 0.05 (*), p < 0.01 (* *), p < 0.001 (* * *). Filled squares in 446 the grid below the bars indicate the identity of mutations that contribute to the coefficient. The 447 names of the mutations, taken from the original publications, are indicated to the left of the grid squares. B) Sub-panels show fraction of variation accounted for by first through fifth order epistatic coefficients for data sets I-IV (colors as in panel A). Fraction described by each order is proportional to area.

- [Figure 6](): **Nonlinear phenotypes distort measured epistatic coefficients.** Sub-panels show correlation plots between epistatic coefficients extracted without accounting for nonlinearity (x-axis) and accounting for linearity (y-axis) for data sets I-IV. Each point is an epistatic coefficient, colored by order. Error bars are standard deviations from bootstrap replicates of each fitting approach.