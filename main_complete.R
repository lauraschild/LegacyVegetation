#run entire repos

rm(list = ls())

#please adjust the working directory to the repository path
setwd(".")

{####working directory check ####
  #script to check the working directory
  cat("\n\t\tINFO\nYou will need to change the working directory to the repository path, i.e. '../LegacyVegetation'. 
    \nAlternatively you can open the R Project located in the repository first and then open this script. The working directory will then be set automatically.")
  
  if(!grepl(pattern = "LegacyVegetation",getwd())){
    warning("It looks like you didn't adjust the working directory to the repository path yet.")
    warning("\n(Unless you renamed the directory. In that case you may ignore this warning.)")
  }
}

#check packages
source("scripts/packages.R")

source("main_REVEALS.R")

source("main_reconstruct_forests.R")
