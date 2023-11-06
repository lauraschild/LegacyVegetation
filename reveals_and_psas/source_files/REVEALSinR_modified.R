checkREVEALS_mod = function (pollen, params, tBasin, dBasin, dwm, n, regionCutoff, 
                          writeresults, verbose){
  if ((dBasin/2) > regionCutoff) 
    stop("Stop: Basin larger than the size of the region (regionCutoff) - that does'nt work!")
  if (!all(colnames(pollen[, -1]) %in% rownames(params), na.rm = FALSE)) 
    stop("STOP: Not all pollen taxa found in the parameter list")
  if (!is.numeric(pollen[, 1])) 
    stop("STOP: some or all values in first column of pollen record (ages/depths) not numeric")
  if (!isTRUE((dwm == "lsm unstable") | (dwm == "gpm neutral") | 
              (dwm == "gpm unstable") | (dwm == "1overd"))) 
    stop("distance weighting method not defined; should be 'LSM unstable', 'gpm neutral', 'gpm unstable' or '1oved' ")
  if (!isTRUE((tBasin == "peatland") | (tBasin == "lake"))) 
    stop("basin type ('tBasin') not defined; should be 'peatland' or 'lake'")
  if (!isTRUE(dBasin == floor(dBasin))) 
    stop("basin size ('dBasin') should be an integer value")
  if (dBasin < 10) 
    stop("basin diameter ('dBasin') should be at least 10 m")
  if (!isTRUE(n == floor(n))) 
    stop("'n' should be an integer value")
  if (n < 1000) 
    message("Warning: for sensible error estimates, 'n' should be at least 1000")
  if (!"fallspeed" %in% tolower(names(params))) 
    stop("STOP: no 'fallspeed' column in the parameters ")
  if (max(params$fallspeed) > 0.15) 
    stop("fallspeed(s) too high (>0.15 m/s), please check")
  if (min(params$fallspeed) < 0.01) 
    stop("fallspeed(s) too low (<0.01 m/s), please check")
  if (!"ppe" %in% tolower(names(params))) 
    stop("STOP: no 'ppe' column in the parameters ")
  if (!"ppe.error" %in% tolower(names(params))) 
    stop("STOP: no 'ppe.error' column in the parameters ")
  return(T)
}


################################################################################

gpm_mod = function (fallSpeed, cz, n, u){
  x = seq(from = 0, to = 100, by = 0.1)
  x = x^3
  y = exp((-4 * fallSpeed * x^(n/2))/(n * u * sqrt(pi) * cz))
  return(loess(y ~ x, span = 0.015, model = FALSE, degree = 2))
}


################################################################################
  
  # see Sugita 2007, Appendix A, K_i for lakes, while x is only is only one 
  # distance here, called in function DispersalFactorK_mod
  #
  # for further insights read Sugita, 1993
  
LakeModel_mod = function (x, dBasin, disModel, stepSize, regionCutoff){
  rBasin = dBasin/2
  r2 = seq(from = rBasin - x, to = 2 * rBasin, by = 2)
  r2Cut = r2
  r2Cut[r2Cut > (rBasin + x)] = rBasin + x
  alpha = acos((r2Cut^2 + x^2 - rBasin^2)/(2 * r2Cut * x))
  propOut = 1 - alpha/pi
  a = diff(propOut)/2
  propOut = propOut[-1] - a
  airborne = predict(disModel, r2)
  influxRing = abs(diff(airborne))
  influxRing = influxRing * propOut
  influx = sum(influxRing) + predict(disModel, 2 * rBasin) - 
    predict(disModel, regionCutoff)
  ringArea = pi * x^2 - pi * (x - stepSize)^2
  return(influx * ringArea)
}


################################################################################

DispersalFactorK_mod = function (vg, tBasin, dBasin, 
                              dwm, regionCutoff, 
                              u = 3, verbose = F){
  rBasin = dBasin/2
  if (dwm == "gpm neutral") {
    disModel = gpm_mod(fallSpeed = vg, cz = 0.12, n = 0.25, 
                    u = u)
  }
  else if (dwm == "gpm unstable") {
    disModel = gpm_mod(fallSpeed = vg, cz = 0.21, n = 0.2, u = u)
  }
  else if (dwm == "1overd") {
    disModel = OneOverD.f()
  }
  else if (dwm == "lsm unstable low") {
    vgr = 200 * round(vg, 2) + 1
    vgr[vgr == 1] = 2
    disModel = lsm_unstable_low[[vgr]]
  }
  else if (dwm == "lsm unstable") {
    vgr = 100 * round(vg, 2) + 1
    disModel = lsm_unstable[[vgr]]
  }
  else stop("no valid distance weighting method selected")
  
  if (tBasin == "peatland") {
    # Dispersal factor K = exp(-b_i * z^gamma) | 
    #   evaluated on R_k (upper) and z_max(lower)
    if (regionCutoff <= 1e+06){ 
      getInflux_mod = function(regionCutoff, 
                            disModel=disModel,
                            rBasin=rBasin,
                            dBasin=dBasin){
        influx = predict(disModel, rBasin) - predict(disModel,regionCutoff) 
      }
    }
    else{
      getInflux_mod = function(regionCutoff, 
                            disModel=disModel,
                            rBasin=rBasin,
                            dBasin=dBasin){
        influx = predict(disModel, rBasin)
      }  
    }
  }
  else if (tBasin == "lake") {
    getInflux_mod = function(regionCutoff, 
                          disModel=disModel,
                          rBasin=rBasin,
                          dBasin=dBasin){
      # changing stepsize within the LakeModel_mod function to rSteps, which it might
      # be meant to, however, does not change the PSA results but 
      # increases computation speed a lot (factor 4 for current settings)
      # only relevant for large lakes
      # if (dBasin >= 20000)
      #   stepSize = 50
      # if (dBasin >= 200000)
      #   stepSize = 150
      
      stepSize = step_size # original version with step_size=2
      if (dBasin >= 2000) 
        stepSize = 5*step_size
      rSteps = seq(from = stepSize, to = rBasin - 1, by = stepSize)
      lakeInflux = do.call(rbind, lapply(rSteps, LakeModel_mod, 
                                          dBasin = dBasin, disModel = disModel, 
                                          stepSize = stepSize, 
                                          regionCutoff = regionCutoff))
      influx = sum(lakeInflux)/(pi * (max(rSteps))^2)
    }
  }
  
  # evaluate pollen influx share-distance relationship
  full_dist = regionCutoff - rBasin
  nr_low = round(region_steps*3/4,digit=0)
  step_size_low =  round(full_dist/(4*(nr_low+1)),digit=0)
  nr_high = round(region_steps*1/4,digit=0)+1
  step_size_high =  round(3*full_dist/(4*nr_high),digit=0)
  dists = c(seq(from=rBasin+step_size_low, to=rBasin+full_dist/4, 
                 by=step_size_low),
             seq(from=rBasin+full_dist/4+step_size_high, to=regionCutoff, 
                 by=step_size_high))
  influxes_over_dists = do.call(rbind, lapply(dists, getInflux_mod, 
                                               disModel = disModel,
                                               rBasin = rBasin,
                                               dBasin = dBasin))
  influxes = data.frame(cbind(dists,influxes_over_dists))
  colnames(influxes) = c("distance","influx")
  return(influxes)
}


################################################################################

rnorm_reveals_mod = function (n, mean = rep(1, n), sds = rep(1, n)){
  x = 0 * mean - 1
  while (any(x < 0)) {
    x = rnorm(n = n, mean = mean, sd = sds)
  }
  return(x)
}
 

################################################################################

REVEALS_mod = function(pollen, n = n, ppes, ppefun_mod = rnorm_reveals_mod, 
                       deposition, pollenfun_mod = rmultinom_reveals_mod, 
                       writeresults = writeresults, verbose = verbose){
  library(ParallelLogger)
    
  if (verbose) {
    #cat(paste("\n\t\t\tAge Slice ",pollen[1]))
    logTrace(paste("\n\t\t\tAge Slice ", pollen[1]))
  }
  coverM = replicate(n = n, {
    ppes[, 3] = ppefun_mod(length(ppes[, "ppe"]), ppes[, "ppe"],
                            ppes[, "ppe.error"])
    pollenPer = pollenfun_mod(n = 1, pollen = pollen[-1])
    disp_ppe = deposition * ppes[, 3]
    s = sum(pollenPer/disp_ppe)
    round(100 * pollenPer/(disp_ppe * s), digits = 3)
  }, simplify = T)
  rownames(coverM) = rownames(pollen[-1])
  if (writeresults) {
    fileName = paste("results_cover", as.character(pollen[1]), 
                      ".csv", sep = "")
    write.csv2(t(coverM), fileName)
  }
  coverMean = apply(coverM, 1, mean)
  coverMedian = apply(coverM, 1, median)
  coverSD = apply(coverM, 1, sd)
  coverq10 = apply(coverM, 1, function(u) quantile(probs = 0.1, 
                                                    x = u))
  coverq90 = apply(coverM, 1, function(u) quantile(probs = 0.9, 
                                                    x = u))
  return(invisible(list(meansim = coverMean, sdsim = coverSD, 
                        mediansim = coverMedian, 
                        q90sim = coverq90, q10sim = coverq10)))
}


REVEALS_mod_vary_PPE_sets = function(pollen, param_grid, deposition, 
                                     pollenfun,
                                     writeresults = writeresults, 
                                     verbose = verbose){
  if (verbose) {
    #cat(paste("\n\t\t\ttime slice:",pollen[1]," BP"))
    logTrace(paste("\n\t\t\tAge Slice ", pollen[1]))
  }
  taxa = colnames(param_grid[,-1])
  coverM = data.frame(t(apply(param_grid, 1, function(u){
    ppes = u[-1]
    # what happens here?
    pollenPer = pollenfun(m = 1, pollen = pollen[-1])
    disp_ppe = deposition * ppes
    s = sum(pollenPer/disp_ppe)
    round(100 * pollenPer/(disp_ppe * s), digits = 3)
  })))
  row.names(coverM) = NULL
  colnames(coverM) = taxa
  coverM["PPE_Var_ID"] = param_grid$PPE_Var_ID
  coverM["Age_BP"] = pollen[1]
  coverM = coverM[,c("Age_BP","PPE_Var_ID", taxa)]

  return(coverM)
}


################################################################################

rmultinom_reveals_mod = function (m = 1, pollen, ...) {
  100 * rmultinom(m, sum(pollen), pollen/sum(pollen))/sum(pollen, na.rm = T)
}


################################################################################

run_reveals_mod = function (pollen, params, 
                            tBasin, dBasin, 
                            sites,
                            dwm = "lsm unstable", dwm_ind=2,
                            n = 1000, regionCutoff = 1e+06, 
                            ppefun_mod = rnorm_reveals_mod, 
                            pollenfun_mod = rmultinom_reveals_mod, 
                            writeresults = FALSE, 
                            verbose = T, site="NA",region="NA",
                            rpp_set="NA", rpp_ind=1, 
                            influx_threshs=c(80), influx_thresh_use=80){
  
  dwm = tolower(dwm)
  tBasin = tolower(tBasin)
  names(params) = tolower(names(params))
  pollen = as.data.frame(pollen)
  params = as.data.frame(params)
  row.names(params) = params[, 1]
  params = params[, -1]
  if (!checkREVEALS_mod(pollen = pollen, params = params, tBasin = tBasin, 
                    dBasin = dBasin, dwm = dwm, n = n, 
                    regionCutoff = regionCutoff, 
                    writeresults = writeresults, 
                    verbose = verbose)) 
    stop("Some parameter is not suitable - check and correct")
  params = params[colnames(pollen[, -1]), ]
  ppes = params[, c("ppe", "ppe.error")]
  vg = params[, "fallspeed"]
  if (verbose) {
    cat(paste("\n\tSite",which(sites==site),"with DatasetID",site))
    cat("\n\tREVEALS:\n\t\t* Calculating deposition factors per species.")
    cat(paste0("\n\t\t* Varying region cutoff to find ",
               "\n\t\t  pollen source radius, where",
               "\n\t\t  on average over taxa ", influx_threshs,"%",
               "\n\t\t  of total pollen influx comes from:"))
  }
  # here, for any fallspeed in vg, compute the influxes for different
  # distances between dBasin and regionCutoff, that afterwards are used
  # to take the 80% from regionCutoff's influx
  influxes = lapply(vg, DispersalFactorK_mod,
                    tBasin = tBasin, dBasin = dBasin,
                    dwm = dwm, regionCutoff = regionCutoff,
                    verbose = verbose)
  # the following is modified. Usually here only deposition is computed
  # the function compute_and_plot_PSA additionally computes Pollen Source Areas
  PSA_results =
    compute_and_plot_PSA(vg=vg, influxes=influxes,
                         influx_threshs=influx_threshs, 
                         params=params,
                         regionCutoff=regionCutoff,
                         dwm=dwm, dwm_ind=dwm_ind,
                         region=region, site=site,
                         rpp_set=rpp_set, rpp_ind=rpp_ind,
                         tBasin=tBasin, dBasin=dBasin)
  deposition = PSA_results[[1]]
  source_radius = PSA_results[[2]]
  source_radius_std = PSA_results[[3]]
  if (verbose){
    cat(paste0("\n\t\t  Taxa Average +/- St. Dev = ",
               "\n\t\t  (", round(source_radius/1000,2), " +/- ", 
               round(source_radius_std/1000,2), ") km"))
    cat("\n\t\t* Reconstruction of vegetation cover:")
  }
  # apply REVEALS seperately to rows (time slices)
  # do this n times with normal distributed PPE variations
  # and average (sd, q10, q90) result
  all = apply(pollen, 1, REVEALS_mod, n = n, ppes = ppes, ppefun_mod = ppefun_mod, 
               deposition = deposition, pollenfun_mod = pollenfun_mod,
               writeresults = writeresults, 
               verbose = verbose)
  results_df = data.frame(dwm, tBasin, 
                           pollen[1], t(sapply(all,
                                               FUN = function(x) x$meansim)),
                           t(sapply(all, FUN = function(x) x$mediansim)), 
                           t(sapply(all, FUN = function(x) x$q90sim)), 
                           t(sapply(all, FUN = function(x) x$q10sim)), 
                           t(sapply(all, FUN = function(x) x$sdsim)))
  colnames(results_df) = c("Distance.weighting", 
                            "Basin.type", 
                            names(pollen)[1], 
                            paste(names(pollen[-1]),
                                  rep(c("mean",
                                        "median",
                                        "q90",
                                        "q10",
                                        "sd"),
                                      each = length(pollen[1, ]) - 1),
                                  sep = "."))
  if (writeresults) {
    nn = names(results_df)
    k = c(1:4, grep(x = nn, pattern = "mediansim"), 
           grep(x = nn, pattern = "q10sim"), grep(x = nn, 
                                                  pattern = "q90sim"))
    write.table(results_df[, k], 
                file = "ResultStatistics_confidence_limits.csv", 
                append = FALSE, row.names = FALSE, col.names = T, 
                sep = ";", dec = ",")
    k = c(1:4, grep(x = nn, pattern = "meansim"), 
           grep(x = nn, pattern = "sdsim"))
    write.table(results_df[, k], file = "ResultStatistics_SD.csv", 
                append = FALSE, row.names = FALSE, col.names = T, 
                sep = ";", dec = ",")
  }
  results =list("reveals_results"=results_df, 
                 "source_radius"=source_radius,
                 "source_radius_std"=source_radius_std) 
  invisible(results)
}



run_reveals_mod_vary_PPE_sets = function(pollen, 
                                         params, param_grid, 
                                         tBasin, dBasin, 
                                         dwm = "lsm unstable", 
                                         n = 1000,
                                         regionCutoff = 1e+05, 
                                         pollenfun = rmultinom_reveals_mod,
                                         writeresults = F, 
                                         verbose = T){
  
  dwm = tolower(dwm)
  tBasin = tolower(tBasin)
  names(params) = tolower(names(params))
  pollen = as.data.frame(pollen)
  params = as.data.frame(params)
  row.names(params) = params[, 1]
  params = params[, -1]
  if (!checkREVEALS(pollen = pollen, params = params, tBasin = tBasin, 
                    dBasin = dBasin, dwm = dwm, n = n, 
                    regionCutoff = regionCutoff, 
                    writeresults = writeresults, verbose = verbose)) 
    stop("Some parameter is not suitable - check and correct")
  params = params[colnames(pollen[, -1]), ]
  vg = params[, "fallspeed"]
  if (verbose) {
    message("REVEALS:\tCalculating deposition factor for each species ...")
  }
  deposition = do.call(rbind, lapply(vg, DispersalFactorK, 
                                     tBasin = tBasin, dBasin = dBasin,
                                     dwm = dwm, regionCutoff = regionCutoff))
  results_df = data.frame(do.call(rbind, apply(pollen, 1, REVEALS_mod_vary_PPE_sets, param_grid = param_grid, 
                       deposition = deposition, pollenfun = pollenfun, 
                       writeresults = writeresults, verbose = verbose)))
  results_df["Distance.weighting"] = dwm
  results_df["Basin.type"] = tBasin
  results_df = results_df[,c("Distance.weighting","Basin.type",
                             "Age_BP","PPE_Var_ID",
                             names(pollen[-1]))]

  invisible(results_df)
}



################################################################################

# modified plotREVEALS to plot reconstructed cover
plotREVEALS_mod_compare_settings = function (d, type = "single", 
                                              ytitle = "Age cal. BP", 
                                              reverse = T,
                                              setting_titles=setting_titles){
  if (!requireNamespace("forcats", quietly = T)) {
    stop("Package forcats needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  age = d[[1]][[3]]
  d.mean = lapply(d,function(x){
    x[grep("mean", names(d[[1]]), value = T)]})
  for (ii in c(1:length(d.mean))){
    names(d.mean[[ii]]) = sub(".mean", "", names(d.mean[[ii]]))
  }
  
  yr = per = low = up = taxa = NULL
  df.means = lapply(d.mean,function(x){
    data.frame(yr = rep(age, ncol(x)), per = as.vector(as.matrix(x)), 
               taxa = as.factor(rep(colnames(x), each = nrow(x))))
    })
  df.mean = df.means[[1]]
  cols = colnames(df.mean); cols[2]=setting_titles[[1]]; colnames(df.mean)=cols
  df.mean[unlist(setting_titles[[2]])] = as.vector(as.matrix(df.means[[2]][,"per"]))
  df.mean[unlist(setting_titles[[3]])] = as.vector(as.matrix(df.means[[3]][,"per"]))
  df.mean[unlist(setting_titles[[4]])] = as.vector(as.matrix(df.means[[4]][,"per"]))
  
  theme_new = ggplot2::theme(#panel.grid.major = element_blank(), 
                              panel.grid.minor.x = element_blank(), 
                              panel.background = element_blank(), 
                              axis.line = element_line(colour = "black"), 
                              strip.text.x = element_text(size = 10, 
                                                          angle = 90, 
                                                          vjust = 0.5, 
                                                          hjust = 0.2), 
                              strip.background = element_rect(color=NA,
                                                              fill="grey98"),
                              #element_blank(), 
                              strip.text.y = element_text(angle = 0), 
                              legend.position = "right", 
                              panel.border = element_blank(), 
                              axis.text.x = element_text(angle = 0, 
                                                         hjust = 0.5))
  p = ggplot2::ggplot(df.mean) +  
    geom_line(aes(yr,NH_NE),
              col = "darkblue", size = 0.5)+{if (reverse)
                scale_x_reverse()
              } +
    geom_line(aes(yr,CF_NE),
              col = "firebrick4", size = 0.5)+{if (reverse)
                scale_x_reverse()
              } +
    geom_line(aes(yr,NH_US),
              col = "steelblue3", size = 0.5)+{if (reverse)
                scale_x_reverse()
              } +
    geom_line(aes(yr,CF_US),
              col = "coral2", size = 0.5)+{if (reverse)
                scale_x_reverse()
              } +
    scale_y_continuous(breaks=seq(0, 100, 10)) + 
    aes(ymin = 11) + 
    xlab(ytitle) + 
    ylab("%") + 
    coord_flip() + 
    theme_new + 
    facet_grid(~forcats::fct_inorder(df.mean$taxa), 
               scales = "free", 
               space = "free")
  return(p)
}

