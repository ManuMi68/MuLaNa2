# MuLaNa2

Muzzio Lab Neuroscience Analysis -2

This is a set of statistical analysis codes written in R.
It is oriented to the analysis of experimental designs, both parametric ANOVA and Robust ANOVA,
including on-way, two-way, three-way, and four-way between and within.
On the other hand, Bayesian statistics and modeling-based analysis are also implemented.

The data is included in the folder "NeuroSpatialData".

The code is organized in an html file that includes all the explanations and details:
"NeuralSpatial.html".
Alternatively the code is also included in the RStudio file:
"NeuralSpatial.Rmd".

The complete code is included for the analysis of all the data of the paper
"Distinct neural mechanisms for heading retrieval and context recognition in the hippocampus during spatial reorientation".
The code can be copied from the different chunks, organized from the Figures on the paper.
As the code executes, the results (files with anovas and graphics) are stored in the folder "NeuroSpatial",
while they are displayed in the output window.

Once in the R program, install the libraries:
install.packages(c("rtf"), repos = "http://cran.r-project.org").
install.packages(c("parallel","data.table","ggplot2","dplyr", "ARTool", "lubridate", "kableExtra","caret","hrbrthemes", "modeest")).
devtools::install_github("kassambara/ggpubr").
For Robust Analysis:
Source Rand Wilcox Robust Analysis Code from https://dornsife.usc.edu/labs/rwilcox/software/.
This program uses the version: source("Rallfun-v40vMM.txt"), 
that corresponds to the version Rallfun-v40.txt with a few improvements
Adapted functions include the suffix "_vMM"

Authors.
Manuel-Miguel Ramos-Álvarez (mramos@ujaen.es) - collaborator, Muzzio Lab -.
See also the list of contributors who participated in this project.

License.
This project is PRIVATE to Muzzio Lab members. All Rights Reserved 2024.
