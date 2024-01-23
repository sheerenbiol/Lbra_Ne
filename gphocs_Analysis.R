##########################################################
#####           Ne estimation - G-PhoCS analysis     #####
#####              Leishmania Braziliensis           #####
#####            Genomic population structure        #####
#####                    South America               #####
##########################################################

### Import functions:
#####################
setwd('~/Dropbox/Projects/Repositories/Lbra_Ne/')
source('./gphocs_Rfunctions.R')

### packages:
#############
library(car)

### Re-set working dir:
#######################
setwd('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/')


#################################
##### MODEL 1: NO MIGRATION #####
#################################
### Calculate ESS:
##################
ess.no.migration <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_ess(mig.mod = 'No', sub_sample = i)
  ess.no.migration[[i]] <- tmps
}

## ESS should be > 200
## Proceeding with chromosomes:
## subset 1: 1-9, 11-14, 16, 17, 19, 21, 22, 24 --> # order in log files for chrs  1-9, 11-14, 16, 17, 19, 21, 22, 24
## subset 2: 1-5, 7-9, 23.                      --> # order in log files for chrs  1-5, 7-9, 19
## subset 3: 1-5, 8, 9, 21, 24                  --> # order in log files for chrs  1-5, 8, 9, 20, 23
## -->
chr.selection <- list()
chr.selection[['1_Subset1']] <- c(1:9, 11:14, 16, 17, 19, 21, 22, 24)
chr.selection[['2_Subset2']] <- c(1:5, 7:9, 19)
chr.selection[['3_Subset3']] <- c(1:5, 8, 9, 20, 23)

### Convert theta values to Ne:
###############################
Ne.no.migration <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_Ne_from_theta(mig.mod = 'No', sub_sample = i, chrs = chr.selection[[i]])
  Ne.no.migration[[i]] <- tmps
}

### Plot Ne per subset:
#######################
Ne.no.migration.plots <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  tmps <- Neplot_subset(theta.data = Ne.no.migration[[i]])
  Ne.no.migration.plots[[i]] <- tmps
}
Ne.no.migration.plots[[1]]
Ne.no.migration.plots[[2]]
Ne.no.migration.plots[[3]] # --> two extreme outliers!!

## Check for outlier chromosomes in Subset3 (potentially remove them):
fit <- aov(Ne.no.migration[[3]][,4] ~ Ne.no.migration[[3]][,3])
outlierTest(fit); influenceIndexPlot(fit, vars = c("Studentized", "Bonf"))
## --> Remove observation 26, 35 in subset 3:
Ne.no.migration[[3]] <- Ne.no.migration[[3]][-c(26,35),]

## Replace subset3 plot:
Ne.no.migration.plots[[3]] <- Neplot_subset(theta.data = Ne.no.migration[[3]])

### Plot Subsets together:
##########################
Ne.no.migration.combined <- do.call(rbind, Ne.no.migration)
Ne.no.migration.combined$sbset <- c(rep('Set 1',nrow(Ne.no.migration[[1]])), rep('Set 2',nrow(Ne.no.migration[[2]])), rep('Set 3',nrow(Ne.no.migration[[3]])))

Ne.no.migration.facetted <- Neplot_all(theta.data.combined = Ne.no.migration.combined, plot.fmt = 'facet'); Ne.no.migration.facetted

### Summary table:
##################
Ne.no.migration.summary <- get_summary_Ne(theta.data.combined = Ne.no.migration.combined)
View(Ne.no.migration.summary)


########################################
##### MODEL 2: Amazon --> Atlantic #####
########################################
### Calculate ESS:
##################
ess.AM2ATL <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_ess(mig.mod = 'AM2ATL', sub_sample = i)
  ess.AM2ATL[[i]] <- tmps
}

## ESS should be > 200
## Proceeding with chromosomes:
## subset 1: 1-5, 8,9,11,12,13,14,17,18,22 --> # order in log files for chrs 1-5, 8,9,11,12,13,14,17,18,22
## subset 2: 1-5, 9, 10.                   --> # order in log files for chrs 1-5, 9, 10
## subset 3: 1-5                           --> # order in log files for chrs 1-5
## -->
chr.selection <- list()
chr.selection[['1_Subset1']] <- c(1:5, 8,9,11,12,13,14,17,18,22)
chr.selection[['2_Subset2']] <- c(1:5, 9, 10)
chr.selection[['3_Subset3']] <- c(1:5)

### Convert theta values to Ne:
###############################
Ne.AM2ATL <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_Ne_from_theta(mig.mod = 'AM2ATL', sub_sample = i, chrs = chr.selection[[i]])
  Ne.AM2ATL[[i]] <- tmps
}

### Plot Ne per subset:
#######################
Ne.AM2ATL.plots <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  tmps <- Neplot_subset(theta.data = Ne.AM2ATL[[i]])
  Ne.AM2ATL.plots[[i]] <- tmps
}
Ne.AM2ATL.plots[[1]]
Ne.AM2ATL.plots[[2]]
Ne.AM2ATL.plots[[3]]

### Plot Subsets together:
##########################
Ne.AM2ATL.combined <- do.call(rbind, Ne.AM2ATL)
Ne.AM2ATL.combined$sbset <- c(rep('Set 1',nrow(Ne.AM2ATL[[1]])), rep('Set 2',nrow(Ne.AM2ATL[[2]])), rep('Set 3',nrow(Ne.AM2ATL[[3]])))

Ne.AM2ATL.facetted <- Neplot_all(theta.data.combined = Ne.AM2ATL.combined, plot.fmt = 'facet'); Ne.AM2ATL.facetted

### Summary table:
##################
Ne.AM2ATL.summary <- get_summary_Ne(theta.data.combined = Ne.AM2ATL.combined)
View(Ne.AM2ATL.summary)


########################################
##### MODEL 3: Atlantic --> Amazon #####
########################################
### Calculate ESS:
##################
ess.ATL2AM <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_ess(mig.mod = 'ATL2AM', sub_sample = i)
  ess.ATL2AM[[i]] <- tmps
}

## ESS should be > 200
## Proceeding with chromosomes:
## subset 1: 1,3,5,7,13,15,16,18,19,21,24,25 --> # order in log files for chrs 1,3,5,7,13,15,16,18,19,21,24,25
## subset 2: 1,4,5,9,15                      --> # order in log files for chrs 1,4,5,9,15
## subset 3: 1-6, 8, 9, 14, 22               --> # order in log files for chrs 1-6, 8, 9, 14 ,22
## -->
chr.selection <- list()
chr.selection[['1_Subset1']] <- c(1,3,5,7,13,15,16,18,19,21,24,25)
chr.selection[['2_Subset2']] <- c(1,4,5,9,15)
chr.selection[['3_Subset3']] <- c(1:6, 8, 9, 14 ,22)

### Convert theta values to Ne:
###############################
Ne.ATL2AM <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_Ne_from_theta(mig.mod = 'ATL2AM', sub_sample = i, chrs = chr.selection[[i]])
  Ne.ATL2AM[[i]] <- tmps
}

### Plot Ne per subset:
#######################
Ne.ATL2AM.plots <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  tmps <- Neplot_subset(theta.data = Ne.ATL2AM[[i]])
  Ne.ATL2AM.plots[[i]] <- tmps
}
Ne.ATL2AM.plots[[1]]
Ne.ATL2AM.plots[[2]]
Ne.ATL2AM.plots[[3]]

### Plot Subsets together:
##########################
Ne.ATL2AM.combined <- do.call(rbind, Ne.ATL2AM)
Ne.ATL2AM.combined$sbset <- c(rep('Set 1',nrow(Ne.ATL2AM[[1]])), rep('Set 2',nrow(Ne.ATL2AM[[2]])), rep('Set 3',nrow(Ne.ATL2AM[[3]])))

Ne.ATL2AM.facetted <- Neplot_all(theta.data.combined = Ne.ATL2AM.combined, plot.fmt = 'facet'); Ne.ATL2AM.facetted

### Summary table:
##################
Ne.ATL2AM.summary <- get_summary_Ne(theta.data.combined = Ne.ATL2AM.combined)
View(Ne.ATL2AM.summary)


##################################
##### MODEL 4: AM <--> ATL   #####
##################################
### Calculate ESS:
##################
ess.Bidir <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_ess(mig.mod = 'BiDir', sub_sample = i)
  ess.Bidir[[i]] <- tmps
}

## ESS should be > 200
## Proceeding with chromosomes:
## subset 1: 1,2,4-6,8,17,23 --> # order in log files for chromosomes 1,2,4-6,8,17,23
## subset 2: 3-5, 9, 21.     --> # order in log files for chromosomes 3-5, 9, 20
## subset 3: 1-4, 8, 9, 22.  --> # order in log files for chromosomes 1-4, 8, 9, 22
## -->
chr.selection <- list()
chr.selection[['1_Subset1']] <- c(1,2,4:6,8,17,23)
chr.selection[['2_Subset2']] <- c(3:5, 9, 20)
chr.selection[['3_Subset3']] <- c(1:4, 8, 9, 22)

### Convert theta values to Ne:
###############################
Ne.Bidir <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  #print(i)
  tmps <- calc_Ne_from_theta(mig.mod = 'BiDir', sub_sample = i, chrs = chr.selection[[i]])
  Ne.Bidir[[i]] <- tmps
}

### Plot Ne per subset:
#######################
Ne.Bidir.plots <- list()
for (i in c('1_Subset1','2_Subset2','3_Subset3')){
  tmps <- Neplot_subset(theta.data = Ne.Bidir[[i]])
  Ne.Bidir.plots[[i]] <- tmps
}
Ne.Bidir.plots[[1]] # --> two extreme outliers!!
Ne.Bidir.plots[[2]]
Ne.Bidir.plots[[3]]

## Check for outlier chromosomes in Subset1 (potentially remove them):
fit <- aov(Ne.Bidir[[1]][,4] ~ Ne.Bidir[[1]][,3])
outlierTest(fit); influenceIndexPlot(fit, vars = c("Studentized", "Bonf"))
cooks.distance(fit); influenceIndexPlot(fit, vars = "Cook")
## --> Remove observation 24, 32 in subset 1:
Ne.Bidir[[1]] <- Ne.Bidir[[1]][-c(24,32),]

## Replace subset3 plot:
Ne.Bidir.plots[[1]] <- Neplot_subset(theta.data = Ne.Bidir[[1]])


### Plot Subsets together:
##########################
Ne.Bidir.combined <- do.call(rbind, Ne.Bidir)
Ne.Bidir.combined$sbset <- c(rep('Set 1',nrow(Ne.Bidir[[1]])), rep('Set 2',nrow(Ne.Bidir[[2]])), rep('Set 3',nrow(Ne.Bidir[[3]])))

Ne.Bidir.facetted <- Neplot_all(theta.data.combined = Ne.Bidir.combined, plot.fmt = 'facet'); Ne.Bidir.facetted

### Summary table:
##################
Ne.Bidir.summary <- get_summary_Ne(theta.data.combined = Ne.Bidir.combined)
View(Ne.Bidir.summary)

##################
### All PLOTS: ###
##################
Ne.no.migration.facetted
Ne.AM2ATL.facetted
Ne.ATL2AM.facetted
Ne.Bidir.facetted

library(ggpubr)
ggarrange(Ne.no.migration.facetted,
          Ne.AM2ATL.facetted,
          Ne.ATL2AM.facetted,
          Ne.Bidir.facetted, nrow = 4)
