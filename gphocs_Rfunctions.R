##########################################################
#####           Ne estimation - G-PhoCS analysis     #####
#####              Leishmania Braziliensis           #####
#####            Genomic population structure        #####
#####                    South America               #####
##########################################################

### Packages:
library(tracerer); library(tidyr)
library(ggplot2); library(dplyr)

#####################
##### FUNCTIONS #####
#####################

### Calculate ESS:
##################
calc_ess <- function(mig.mod = c('No', 'AM2ATL', 'ATL2AM','BiDir'), sub_sample = c('1_Subset1','2_Subset2','3_Subset3')){
  print(" Function to Calculate the ESS (effective Sample size) for each G-PhoCS run for a particular sample subset.")
  print("IMPORTANT: empty log files should be removed from the log-file directory - otherwise the function will give an error message")
  print("INPUT:")
  print("mig.mod     - Migration model that was used OPTIONS: No, AM2ATL, ATL2AM, BiDir")
  print("sub_sample  - certain sub-sample of individuals that was used in the G-PhoCS run")

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
      rownames(tmp.ess) <- c(paste0('chr_',logs[[i]]))
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
      rownames(tmp.ess) <- c(paste0('chr_',logs[[i]]))
      ess.all[[i]] <- tmp.ess
    }

  } else if (mig.mod == 'ATL2AM'){
    main_path <- './3_GPhoCS_output/3_ATL2AM/'
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
      rownames(tmp.ess) <- c(paste0('chr_',logs[[i]]))
      ess.all[[i]] <- tmp.ess
    }

  } else if (mig.mod == 'BiDir'){
    main_path <- './3_GPhoCS_output/4_ATLAM_Bidir/'
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
      rownames(tmp.ess) <- c(paste0('chr_',logs[[i]]))
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
  print(" Calculating the Ne from the estimated theta values using Ne = theta / (4* 1.99e-09).")
  print("IMPORTANT: empty log files should be removed from the log-file directory - otherwise the function will give an error message.")
  print("IMPORTANT: !! Only select chromosomes that show ESS > 200 for all inferred theta values!!")
  print("INPUT:")
  print("mig.mod     - Migration model that was used OPTIONS: No, AM2ATL, ATL2AM, BiDir")
  print("sub_sample  - certain sub-sample of individuals that was used in the G-PhoCS run")
  print("chrs        - vector of chromosome with ESS > 200 for all theta values.")

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
    main_path <- './3_GPhoCS_output/3_ATL2AM/'
    logs <- list.files(paste0(main_path,sub_sample), pattern = 'gphocs.log')

    ## Calculate Ne:
    thetas.all <- list()
    for (i in chrs){
      tmp <- parse_beast_tracelog_file(paste0(main_path,sub_sample,'/',logs[[i]]))
      if (nrow(tmp) == 0) {next}
      tmp <- tmp[,-10]
      tmp$chr <- c(rep(paste0(i),nrow(tmp)))
      colnames(tmp) <- c( "Sample" ,"theta_WAM" , "theta_CAM" , "theta_ATL" , "theta_AM" ,  "theta_root", "tau_AM" , "tau_root",
                          "m_ATL2AM", "Data.ld.ln" ,"Full.ld.ln",'chrom_group')

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
    main_path <- './3_GPhoCS_output/4_ATLAM_Bidir/'
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
Neplot_subset <- function(theta.data = NULL){
  print("Function to plot the estimated Ne for a certain model for a certain sample subset")
  print("INPUT:")
  print("theta.data - theta dataframe in long format. output of calc_Ne_from_theta")

  ## Set the order:
  level_order <- c("theta_WAM", "theta_CAM", "theta_AM" , "theta_ATL")
  ## Plot:
  subset_plot <- ggplot(theta.data, aes(x=factor(Ne_pop, level = level_order), y=Ne_value)) +
    geom_violin(data=theta.data, aes(col=Ne_pop, fill=Ne_pop), alpha=0.2, width=0.7)+
    scale_fill_manual(values = c("#A2662C", "#4E87B8", "#723957", "gray50")) +
    #geom_jitter(data=theta.data, aes(col = Ne_pop ), alpha=0.9, width=0.15, size = 3, stroke=1, shape=16) +
    scale_color_manual(values = c("#A2662C", "#4E87B8", "#723957", "gray50")) +
    geom_point(stat="summary", fun.y="mean_se", size = 3, shape=21, stroke = 1.5) +
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0, linewidth=1) +
    labs(x="Populations", y="Effective population size (No. inds.)") +
    theme_bw() +
    theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5))

  return(subset_plot)
}

### Plot Ne with all subsets - facetted/averaged:
########################################
Neplot_all <- function(theta.data.combined = NULL, plot.fmt = c('facet','average')){
  print("Function to plot the estimated Ne for a certain model for a certain sample subset")
  print("INPUT:")
  print("theta.data.combined - theta dataframe in long format. output of calc_Ne_from_theta")

  ## Set the order:
  level_order <- c("theta_WAM", "theta_CAM", "theta_AM" , "theta_ATL")

  ## Facet Plot:
  if (plot.fmt == 'facet'){
    facet_plot <- ggplot(theta.data.combined, aes(x=sbset, y=Ne_value)) +
      geom_violin(aes(color = Ne_pop, fill=Ne_pop), alpha=0.2)+
      scale_fill_manual(values =  c("#A2662C", "#4E87B8", "#723957", "gray50")) +
      facet_wrap(~ factor(Ne_pop, level = level_order), ncol=4) +
      scale_color_manual(values = c("#A2662C", "#4E87B8", "#723957", "gray50"))+
      #geom_jitter(aes(color = Ne_pop ), alpha=0.8, width=0.2, size = 3, stroke=1, shape=16) +
      geom_point(stat="summary", fun.y="mean_se", size = 3, shape=21, stroke = 1.5, position = position_dodge2(.5)) +
      geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0, linewidth=1, position = position_dodge2(.5)) +
      labs(x="Populations", y="Effective population size (No. inds.)") +
      theme_bw() +
      theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5))
    return(facet_plot)

    ## Averaged plot:
  } else if (plot.fmt == 'average'){
    average_plot <- ggplot(theta.data.combined, aes(x=factor(Ne_pop, level = level_order), y=Ne_value)) +
      geom_violin(aes(color = Ne_pop, fill=Ne_pop), alpha=0.2)+
      scale_fill_manual(values =  c("#A2662C", "#4E87B8", "#723957", "gray50")) +
      scale_color_manual(values = c("#A2662C", "#4E87B8", "#723957", "gray50"))+
      #geom_jitter(aes(color = Ne_pop ), alpha=0.8, width=0.2, size = 3, stroke=1, shape=16) +
      geom_point(stat="summary", fun.y="mean_se", size = 3, shape=21, stroke = 1.5, position = position_dodge2(.5)) +
      geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0, linewidth=1, position = position_dodge2(.5)) +
      labs(x="Populations", y="Effective population size (No. inds.)") +
      theme_bw() +
      theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5))
    return(average_plot)
  }
}

### Create summary table - mean Ne + std error:
###############################################
get_summary_Ne <- function(theta.data.combined = NULL){
  print("Function to extract summary statistics from Ne estimations for a certain migration model and its different subsets")
  print("INPUT:")
  print("theta.data.combined - theta dataframe in long format. output of calc_Ne_from_theta")

  tmp.all <- list()
  lvl <- levels(theta.data.combined[,3])
  for (i in 1:length(levels(theta.data.combined[,3]))){
    tmp.model <- subset(theta.data.combined, theta.data.combined[,3] == levels(theta.data.combined[,3])[[i]])
    tmp.summary <- tmp.model %>%
      group_by(sbset) %>%
      summarise(av = mean(Ne_value), se = sd(Ne_value)/sqrt(length(Ne_value)))
    tmp.summary$mod <- c(rep(lvl[i],nrow(tmp.summary)))
    tmp.all[[i]] <- tmp.summary
  }
  tmp.all <- do.call(rbind, tmp.all)
}
