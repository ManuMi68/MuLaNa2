# MuLaNa2

Muzzio Lab Neuroscience Analysis -2

This is a set of statistical analysis codes written in R. It is oriented towards the analysis of experimental designs, 
both parametric ANOVA and robust ANOVA, including one-way, two-way, three-way and four-way between and within. 
On the other hand, Bayesian statistics and modelling based analysis are also implemented.

The data are included in the NeuroSpatialData folder.

The code is organised in an html file that contains all the explanations and details: "NeuralSpatial.html".
Alternatively, the code is included in the RStudio file: "NeuralSpatial.Rmd".

The NeuralSpatial.html file is designed as a tutorial with instructions and explanations of the data, analysis commands,
figure generation, descriptive statistics (including Effect Size) to facilitate potential future meta-analyses or power estimations.

The complete code is included for the analysis of all data from the paper 
"Distinct neural mechanisms for heading retrieval and context recognition in the hippocampus during spatial reorientation".
The code can be copied from the various chunks organised by the figures in the paper. 
When the code is executed, the results (files with anovas and graphics) are stored in the "NeuroSpatial" folder,
while they are displayed in the output window.

Once in the R program, install the libraries:
install.packages(c("rtf"), repos = "http://cran.r-project.org").
install.packages(c("parallel","data.table","ggplot2","dplyr", "ARTool", "lubridate", "kableExtra","caret","hrbrthemes", "modeest")).
devtools::install_github("kassambara/ggpubr").
For Robust Analysis:
Source Rand Wilcox Robust Analysis code from https://dornsife.usc.edu/labs/rwilcox/software/.
This program uses the version: source("Rallfun-v40vMM.txt"), 
that corresponds to the version Rallfun-v40.txt with some improvements
Adapted functions have the suffix "_vMM"

Authors.
Manuel-Miguel Ramos-Álvarez (mramos@ujaen.es) - collaborator, Isabel Muzzio Lab -.
See also the list of contributors who participated in this project.

Licence. This project is PRIVATE to the members of Muzzio Lab. All rights reserved 2024.

