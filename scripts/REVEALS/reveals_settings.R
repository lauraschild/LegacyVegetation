#PARAMETERS####
{
  n_ppe_variations <- 2000
  region_steps <- 40 # 40 for final runs
  regionCutoff <- 1e+06 # m
  compare_k_factors <- F
  verbose <- TRUE
  use_fixed_steps_lake_model <- T # modified version by Peter as Theuerkaufs parametrization produced unnecessary many steps for large lakes
  n_steps_lake_model <- 500 # 500 for final runs # tested that this is enough steps between basin center and basin radius such that influx integration over basin area does not deviate more than 1% from Theuerkaufs solution
  #set functions to use
  dwm = "gpm unstable"
  pollenfun = rmultinom_reveals
  ppefun = rnorm_reveals
}