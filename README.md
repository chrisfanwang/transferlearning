# transferlearning

This repository contains code and data used in the paper Transfer learning for piecewise-constant mean estimation: \\ Optimality, $\ell_1$- and $\ell_0$-penalisation. Detection by Fan Wang and Yi Yu.

All experiments and analysis are conducted in R (version >= 4.1).

Dependencies
Some R-packages are required to reproduce all results. To install them, one can run the following lines:

list.of.packages <- c("genlasso", "changepoints", "latex2exp", "egg", "tictoc", "ggplot2",
                      "plotrix")
                      
                      
new_packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
for (package in new_packages){
  install.packages(package, repos='http://cran.us.r-project.org', quiet = TRUE)
}

