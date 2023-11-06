#script to run the optimization
library(optimParallel)
library(parallel)

source("scripts/supportive/optimization_function.R", echo=FALSE)

#run optim 
results <- optim(par = parameters$PPE[1:m],
                  fn = optimize,
                  gr = NULL,
                 method = "L-BFGS",
                  lower = bounds$lower[1:m],
                  upper = bounds$upper[1:m],
                  control = list(maxit = 10000))
