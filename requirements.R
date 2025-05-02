packages <- c("VGAM", "MASS", "tictoc", "emdbook",
              "LaplacesDemon", "truncnorm", "optparse",
              "dplyr", "ggplot2", "reshape2", "gridExtra",
              "gtable", "writexl", "sns",
              "doParallel", "foreach", "invgamma", "jsonlite")
install.packages(setdiff(packages, rownames(installed.packages())))