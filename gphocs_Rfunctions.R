##########################################################
#####           Ne estimation - G-PhoCS analysis     #####
#####              Leishmania Braziliensis           #####
#####            Genomic population structure        #####
#####                    South America               #####
##########################################################
##### GOAL: FUNCTIONS for XXX.R

### Calculate ESS:
##################
calc_ess <- function(mig.mod = c('No', 'AM2ATL', 'ATL2AM','BiDir'), sub_sample = c('1_Subset1','2_Subset2','3_Subset3')){
  print(" Function to Calculate the ESS (effective Sample size) for each G-PhoCS run for a particular sample subset.
  IMPORTANT: empty log files should be removed from the log-file directory - otherwise the function will give an error message

  INPUT:
  mig.mod     - Migration model that was used OPTIONS: No, AM2ATL, ATL2AM, BiDir
  sub_sample  - certain sub-sample of individuals that was used in the G-PhoCS run
  ")
  ## Required package:
  library(tracerer)

  ## Set paths for the respective model that was used
  if(mig.mod == 'No') {
    main_path <- './3_GPhoCS_output/1_NoMigration/'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate ESS:
    ess.all <- list()
    for (i in 1:length(logs)){
      tmp <- parse_beast_tracelog_file(paste0(main_path, sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-9]
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "Data.ld.ln" ,"Full.ld.ln")
      tmp$Sample <- as.numeric(tmp$Sample)
      tmp.ess <- calc_esses(tmp, sample_interval = 1000)
      rownames(tmp.ess) <- c(paste0('chr_',log.files[[i]]))
      ess.all[[i]] <- tmp.ess
    }

  } else if (mig.mod == 'AM2ATL') {
    main_path <- './3_GPhoCS_output/2_AM2ATL/'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate ESS:
    ess.all <- list()
    for (i in 1:length(logs)){
      tmp <- parse_beast_tracelog_file(paste0(main_path, sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_AM2ATL","Data.ld.ln" ,"Full.ld.ln")
      tmp$Sample <- as.numeric(tmp$Sample)
      tmp.ess <- calc_esses(tmp, sample_interval = 1000)
      rownames(tmp.ess) <- c(paste0('chr_',log.files[[i]]))
      ess.all[[i]] <- tmp.ess
    }

  } else if (mig.mod == 'ATL2AM'){
    main_path <- './3_GPhoCS_output/3_ATL2AM'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate ESS:
    ess.all <- list()
    for (i in 1:length(logs)){
      tmp <- parse_beast_tracelog_file(paste0(main_path, sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_ATL2AM","Data.ld.ln" ,"Full.ld.ln")
      tmp$Sample <- as.numeric(tmp$Sample)
      tmp.ess <- calc_esses(tmp, sample_interval = 1000)
      rownames(tmp.ess) <- c(paste0('chr_',log.files[[i]]))
      ess.all[[i]] <- tmp.ess
    }

  } else if (mig.mod == 'BiDir'){
    main_path <- './3_GPhoCS_output/4_ATLAM_Bidir'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate ESS:
    ess.all <- list()
    for (i in 1:length(logs)){
      tmp <- parse_beast_tracelog_file(paste0(main_path, sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-11]
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_ATL2AM","m_AM2ATL","Data.ld.ln" ,"Full.ld.ln")
      tmp$Sample <- as.numeric(tmp$Sample)
      tmp.ess <- calc_esses(tmp, sample_interval = 1000)
      rownames(tmp.ess) <- c(paste0('chr_',log.files[[i]]))
      ess.all[[i]] <- tmp.ess
    }
  }
  ess.all <- do.call(rbind, ess.all)
  return(ess.all)
}

### Calculate Ne from theta:
############################
## Formula: Ne = theta / (4* 1.99e-09)
calc_Ne_from_theta <- function(mig.mod = c('No', 'AM2ATL', 'ATL2AM','BiDir'),
                               sub_sample = c('1_Subset1','2_Subset2','3_Subset3'),
                               chrs = NULL){
  print(" Calculating the Ne from the estimated theta values using Ne = theta / (4* 1.99e-09).
  IMPORTANT: empty log files should be removed from the log-file directory - otherwise the function will give an error message.
  IMPORTANT: !! Only select chromosomes that show ESS > 200 for all inferred theta values!!

INPUT:
  mig.mod     - Migration model that was used OPTIONS: No, AM2ATL, ATL2AM, BiDir
  sub_sample  - certain sub-sample of individuals that was used in the G-PhoCS run
  chrs        - vector of chromosome with ESS > 200 for all theta values.
  ")

  ## Set paths for the respective model that was used
  if(mig.mod == 'No') {
    main_path <- './3_GPhoCS_output/1_NoMigration/'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate Ne:
    thetas.all <- list()
    for (i in chrs){
      tmp <- parse_beast_tracelog_file(paste0(main_path,sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-9]
      tmp$chr <- c(rep(paste0(i),nrow(tmp)))
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "Data.ld.ln" ,"Full.ld.ln",'chrom_group')

      tmp.theta <- tmp[,c(2:6,11)]
      for (cl in 1:5){
        tmp.theta[,cl] <- tmp.theta[,cl]/(4 * (1.99e-09))
      }
      tmp.theta <- (c(mean(tmp.theta$theta_WAM), mean(tmp.theta$theta_CAM), mean(tmp.theta$theta_ATL),
                      mean(tmp.theta$theta_AM), mean(tmp.theta$theta_root), paste0(i)))
      thetas.all[[i]] <- tmp.theta
    }
    thetas.all <- do.call(rbind, thetas.all)
    colnames(thetas.all) <- c("theta_WAM" , "theta_CAM", "theta_ATL","theta_AM", "theta_root",'chrom_group')
    thetas.all <- as.data.frame(thetas.all)
    ## Conversion to long format:
    thetas.all.long <- gather(thetas.all, Ne_pop, Ne_value, c(theta_WAM,theta_CAM,theta_ATL,theta_AM), factor_key = T)
    thetas.all.long$Ne_value <- as.numeric(thetas.all.long$Ne_value)

  } else if (mig.mod == 'AM2ATL') {
    main_path <- './3_GPhoCS_output/2_AM2ATL/'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate Ne:
    thetas.all <- list()
    for (i in chrs){
      tmp <- parse_beast_tracelog_file(paste0(main_path,sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      tmp$chr <- c(rep(paste0(i),nrow(tmp)))
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_AM2ATL", "Data.ld.ln" ,"Full.ld.ln",'chrom_group')

      tmp.theta <- tmp[,c(2:6,12)]
      for (cl in 1:5){
        tmp.theta[,cl] <- tmp.theta[,cl]/(4 * (1.99e-09))
      }
      tmp.theta <- (c(mean(tmp.theta$theta_WAM), mean(tmp.theta$theta_CAM), mean(tmp.theta$theta_ATL),
                      mean(tmp.theta$theta_AM), mean(tmp.theta$theta_root), paste0(i)))
      thetas.all[[i]] <- tmp.theta
    }
    thetas.all <- do.call(rbind, thetas.all)
    colnames(thetas.all) <- c("theta_WAM" , "theta_CAM", "theta_ATL","theta_AM", "theta_root",'chrom_group')
    thetas.all <- as.data.frame(thetas.all)
    ## Conversion to long format:
    thetas.all.long <- gather(thetas.all, Ne_pop, Ne_value, c(theta_WAM,theta_CAM,theta_ATL,theta_AM), factor_key = T)
    thetas.all.long$Ne_value <- as.numeric(thetas.all.long$Ne_value)

  } else if (mig.mod == 'ATL2AM'){
    main_path <- './3_GPhoCS_output/3_ATL2AM'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate Ne:
    thetas.all <- list()
    for (i in chrs){
      tmp <- parse_beast_tracelog_file(paste0(main_path,sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      tmp$chr <- c(rep(paste0(i),nrow(tmp)))
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_AM2ATL", "Data.ld.ln" ,"Full.ld.ln",'chrom_group')

      tmp.theta <- tmp[,c(2:6,12)]
      for (cl in 1:5){
        tmp.theta[,cl] <- tmp.theta[,cl]/(4 * (1.99e-09))
      }
      tmp.theta <- (c(mean(tmp.theta$theta_WAM), mean(tmp.theta$theta_CAM), mean(tmp.theta$theta_ATL),
                      mean(tmp.theta$theta_AM), mean(tmp.theta$theta_root), paste0(i)))
      thetas.all[[i]] <- tmp.theta
    }
    thetas.all <- do.call(rbind, thetas.all)
    colnames(thetas.all) <- c("theta_WAM" , "theta_CAM", "theta_ATL","theta_AM", "theta_root",'chrom_group')
    thetas.all <- as.data.frame(thetas.all)
    ## Conversion to long format:
    thetas.all.long <- gather(thetas.all, Ne_pop, Ne_value, c(theta_WAM,theta_CAM,theta_ATL,theta_AM), factor_key = T)
    thetas.all.long$Ne_value <- as.numeric(thetas.all.long$Ne_value)

  } else if (mig.mod == 'BiDir'){
    main_path <- './3_GPhoCS_output/4_ATLAM_Bidir'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate Ne:
    thetas.all <- list()
    for (i in chrs){
      tmp <- parse_beast_tracelog_file(paste0(main_path,sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      tmp$chr <- c(rep(paste0(i),nrow(tmp)))
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_AM2ATL", "Data.ld.ln" ,"Full.ld.ln",'chrom_group')

      tmp.theta <- tmp[,c(2:6,12)]
      for (cl in 1:5){
        tmp.theta[,cl] <- tmp.theta[,cl]/(4 * (1.99e-09))
      }
      tmp.theta <- (c(mean(tmp.theta$theta_WAM), mean(tmp.theta$theta_CAM), mean(tmp.theta$theta_ATL),
                      mean(tmp.theta$theta_AM), mean(tmp.theta$theta_root), paste0(i)))
      thetas.all[[i]] <- tmp.theta
    }
    thetas.all <- do.call(rbind, thetas.all)
    colnames(thetas.all) <- c("theta_WAM" , "theta_CAM", "theta_ATL","theta_AM", "theta_root",'chrom_group')
    thetas.all <- as.data.frame(thetas.all)
    ## Conversion to long format:
    thetas.all.long <- gather(thetas.all, Ne_pop, Ne_value, c(theta_WAM,theta_CAM,theta_ATL,theta_AM), factor_key = T)
    thetas.all.long$Ne_value <- as.numeric(thetas.all.long$Ne_value)
  }
  return(thetas.all.long)
}


### Plot Ne per sample sub-set:
###############################

### Plot Ne with all subsets - facetted:
########################################

### Plot Ne with all subsets - averaged:
########################################

### Create summary table - mean Ne + std error:
###############################################

