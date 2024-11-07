##########################
##### MSMC2 ANALYSES #####
##########################

### Packages:
library(ggplot2);library(ggforce)

### Checking MSMC results for all-population-in-one input file:
###############################################################
### Checking Log-Likelihood:
log.files <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/', pattern = 'EM500.MSMC2.out.loop.txt')

## --> Showing that Log Likelihood is stable at 1000 iterations
par(mfrow=c(3,3))
for (file in 1:length(log.files)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/',log.files[[file]]), header = F)
  plot(x=rownames(tmp), y=tmp$V2,
       #main=paste0(log.files[[file]]),
       xlab='EM Iterations', ylab='Log Likelihood', main=paste0(log.files[[file]]),
       xlim=c(0,500), ylim=c(min(tmp$V2)-1000,max(tmp$V2)+1000))
}
par(mfrow=c(1,1))

### Data visualization
data.files <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/', pattern = 'EM500.MSMC2.out.final.txt')
data.combined <- list()
for ( file in 1:9){ #1:length(data.files)
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/',data.files[[file]]), header = T)
  if (file == 1) { tmp$pop <- c(rep('ATL-1',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 2) {tmp$pop <- c(rep('ATL-2',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 3) {tmp$pop <- c(rep('ATL-3',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 4) {tmp$pop <- c(rep('CAM-1',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 5) {tmp$pop <- c(rep('CAM-2',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 6) {tmp$pop <- c(rep('CAM-3',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 7) {tmp$pop <- c(rep('WAM-1',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 8) {tmp$pop <- c(rep('WAM-2',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 9) {tmp$pop <- c(rep('WAM-3',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  #else if (file == 10) {tmp$pop <- c(rep('WAM-1',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  #else if (file == 11) {tmp$pop <- c(rep('WAM-2',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  #else if (file == 12) {tmp$pop <- c(rep('WAM-3',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  
  data.combined[[file]] <- tmp
}

mu <- 1.99e-09 #unit bp per generation

## plot subset 1:
plot(data.combined[[4]]$left_time_boundary/mu, log10((1/data.combined[[4]]$lambda)/(2*mu)), #ylim=c(5,11),
     xlim=c(0,10000000),
     type="n", xlab="Generations", ylab="effective population size")
#lines(data.combined[[1]]$left_time_boundary/mu, log10((1/data.combined[[1]]$lambda)/(2*mu)), type="s",lwd=3, col="gray50")
lines(data.combined[[4]]$left_time_boundary/mu, log10((1/data.combined[[4]]$lambda)/(2*mu)), type="s",lwd=3, col="#723957")
lines(data.combined[[7]]$left_time_boundary/mu, log10((1/data.combined[[7]]$lambda)/(2*mu)), type="s",lwd=3, col="#4E87B8")
lines(data.combined[[10]]$left_time_boundary/mu, log10((1/data.combined[[10]]$lambda)/(2*mu)), type="s",lwd=3, col="#A2662C")
legend("bottomright",legend=c("ATL", "CAM","WAM"), col=c("#723957", "#4E87B8","#A2662C"), lwd=c(3,3,3), lty=c(1,1,1))

## plot subset 2:
plot(data.combined[[2]]$left_time_boundary/mu, log10((1/data.combined[[2]]$lambda)/(2*mu)), ylim=c(5,10.5),
     xlim=c(0,10000000),
     type="n", xlab="Generations", ylab="effective population size")
#lines(data.combined[[2]]$left_time_boundary/mu, log10((1/data.combined[[2]]$lambda)/(2*mu)), type="s",lwd=3, col="gray50")
lines(data.combined[[5]]$left_time_boundary/mu, log10((1/data.combined[[5]]$lambda)/(2*mu)), type="s",lwd=3, col="#723957")
lines(data.combined[[8]]$left_time_boundary/mu, log10((1/data.combined[[8]]$lambda)/(2*mu)), type="s",lwd=3, col="#4E87B8")
lines(data.combined[[11]]$left_time_boundary/mu, log10((1/data.combined[[11]]$lambda)/(2*mu)), type="s",lwd=3, col="#A2662C")
legend("bottomright",legend=c("ATL", "CAM","WAM"), col=c("#723957", "#4E87B8","#A2662C"), lwd=c(3,3,3), lty=c(1,1,1))

## plot subset 3:
plot(data.combined[[3]]$left_time_boundary/mu, log10((1/data.combined[[3]]$lambda)/(2*mu)), ylim=c(5,10.5),
     xlim=c(0,10000000),
     type="n", xlab="Generations", ylab="effective population size")
#lines(data.combined[[3]]$left_time_boundary/mu, log10((1/data.combined[[3]]$lambda)/(2*mu)), type="s",lwd=3, col="gray50")
lines(data.combined[[6]]$left_time_boundary/mu, log10((1/data.combined[[6]]$lambda)/(2*mu)), type="s",lwd=3, col="#723957")
lines(data.combined[[9]]$left_time_boundary/mu, log10((1/data.combined[[9]]$lambda)/(2*mu)), type="s",lwd=3, col="#4E87B8")
lines(data.combined[[12]]$left_time_boundary/mu, log10((1/data.combined[[12]]$lambda)/(2*mu)), type="s",lwd=3, col="#A2662C")
legend("bottomright",legend=c("ATL", "CAM","WAM"), col=c("#723957", "#4E87B8","#A2662C"), lwd=c(3,3,3), lty=c(1,1,1))

## Nicer plots with ggplot:
data.combined <- do.call(rbind, data.combined)

mu <- 1.99e-09 #unit bp per generation
#gen <- 0.10 # unit years per generation
#Gen time t. brucei --> 7 - 10 gens per year --> 0.14 t.em. 0.10 years /gen

time <- (data.combined$left_time_boundary/mu)
data.combined$time <- time
pop.size <- (1/data.combined$lambda)/(2*mu)
data.combined$popsize <- pop.size


ggplot(data=data.combined[which(data.combined$time <= 1e07),], aes(x=(time), y=log10(popsize), fill=pop2, col=pop2))+
  #geom_line(linewidth = 1.5, alpha=0.75) +
  geom_smooth(method = "loess",formula = y~x, se = F, span = 0.3 ) +
  scale_color_manual(values = c("gray50","#723957","#4E87B8","#A2662C","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  stat_cor(method = "pearson", label.x = 3, label.y = 8)+
  stat_summary_bin(fun.data='mean_cl_boot', size=0.75)+
  xlab('No. generations') + ylab('Log. Effective population size') #+ 
#xlim(-1,8)

ggplot(data=subset(data.combined, data.combined$pop2 != 'AM'), aes(x=(time), y=log10(popsize), fill=pop, col=pop))+
  #geom_line(linewidth = 1.5, alpha=0.75) +
  geom_step(subset(data.combined, data.combined$pop2 != 'AM'), mapping=aes(x=time, y=log10(popsize)), direction = 'hv', linetype=1, linewidth=1.5)+
  #geom_smooth(method = "loess",formula = y~x, se = F, span = 0.3 ) +
  scale_color_manual(values = c("#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #scale_color_manual(values = c("#723957","#4E87B8","#A2662C","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.75)+
  xlab('No. generations') + ylab('Log. Effective population size') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
  facet_zoom(xlim=c(0,2.5e07))
#xlim(-1,8)

ggplot(data=subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 2.5e07),], aes(x=(time), y=log10(popsize), fill=pop, col=pop))+
  #geom_line(linewidth = 1.5, alpha=0.75) +
  geom_step(subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 2.5e07),], mapping=aes(x=time, y=log10(popsize)), direction = 'hv', linetype=1, linewidth=1.5)+
  scale_color_manual(values = c("#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #scale_color_manual(values = c("#723957","#4E87B8","#A2662C","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.5)+
  xlab('No. generations') + ylab('Log. Effective population size') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
#facet_zoom(xlim=c(0,1e07))
#facet_zoom(xlim = c(0,1.25e06))
#facet_zoom(xlim = c(1.25e06,5e06))
#facet_zoom(xlim = c(5e06,10e06))
facet_zoom(xlim = c(0,5e05))
#xlim(-1,8)

#######=========================================================================
### Checking MSMC results for population-separate input files:
##############################################################
### Checking Log-Likelihood:
#data <- read.table('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/Pop_smallgenomeP/ATL.subset3.ALLsmallgen.EM1000.MSMC2.out.loop.txt', header = F)
#head(data)
#plot(x=rownames(data), y=data$V2)
log.files <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/Pop_smallgenomeP/', pattern = 'EM1000_MSMC2_out.loop.txt')

## --> Showing that Log Likelihood is stabel at 1000 iterationns
par(mfrow=c(3,3))
for (file in 1:length(log.files)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/Pop_smallgenomeP/',log.files[[file]]), header = F)
  plot(x=rownames(tmp), y=tmp$V2,
       #main=paste0(log.files[[file]]),
       xlab='EM Iterations', ylab='Log Likelihood',
       xlim=c(0,1000), ylim=c(min(tmp$V2)-1000,max(tmp$V2)+1000))
}
par(mfrow=c(1,1))

### Data visualization
data.files <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/Pop_smallgenomeP/', pattern = 'EM1000_MSMC2_out.final.txt')
data.combined <- list()
for ( file in 1:length(data.files)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/Pop_smallgenomeP/',data.files[[file]]), header = T)
  if (file == 1) { tmp$pop <- c(rep('ATL-1',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 2) {tmp$pop <- c(rep('ATL-2',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 3) {tmp$pop <- c(rep('ATL-3',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 4) {tmp$pop <- c(rep('CAM-1',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 5) {tmp$pop <- c(rep('CAM-2',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 6) {tmp$pop <- c(rep('CAM-3',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 7) {tmp$pop <- c(rep('WAM-1',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 8) {tmp$pop <- c(rep('WAM-2',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 9) {tmp$pop <- c(rep('WAM-3',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  
  data.combined[[file]] <- tmp
}

mu <- 1.99e-09 #unit bp per generation

data.combined <- do.call(rbind, data.combined)

time <- (data.combined$left_time_boundary/mu)
data.combined$time <- time
pop.size <- (1/data.combined$lambda)/(2*mu)
data.combined$popsize <- pop.size

ggplot(data=subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 1e07),], aes(x=(time), y=log10(popsize), fill=pop2, col=pop2))+
  #geom_line(linewidth = 1.5, alpha=0.75) +
  geom_step(subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 1e07),], mapping=aes(x=time, y=log10(popsize)), direction = 'hv', linetype=1, linewidth=1.5)+
  #scale_color_manual(values = c("#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  scale_color_manual(values = c("#723957","#4E87B8","#A2662C","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.5)+
  xlab('No. generations') + ylab('Log. Effective population size') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
#facet_zoom(xlim=c(0,1e07))
facet_zoom(xlim = c(0,1.25e06))
#facet_zoom(xlim = c(1.25e06,5e06))
#facet_zoom(xlim = c(5e06,10e06))
#facet_zoom(xlim = c(0,1e06))
#xlim(-1,8)


### Data visualization
data.files <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/', recursive = T, pattern = '.final.txt')
data.files <- data.files[-c(1:3,13:18)]
data.combined <- list()
for ( file in 1:length(data.files)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/',data.files[[file]]), header = T)
  if (file == 1) { tmp$pop <- c(rep('ATL-1',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 2) {tmp$pop <- c(rep('ATL-2',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 3) {tmp$pop <- c(rep('ATL-3',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 4) {tmp$pop <- c(rep('CAM-1',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 5) {tmp$pop <- c(rep('CAM-2',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 6) {tmp$pop <- c(rep('CAM-3',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 7) {tmp$pop <- c(rep('WAM-1',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 8) {tmp$pop <- c(rep('WAM-2',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 9) {tmp$pop <- c(rep('WAM-3',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 10) { tmp$pop <- c(rep('ATL-4',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 11) {tmp$pop <- c(rep('ATL-5',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 12) {tmp$pop <- c(rep('ATL-6',nrow(tmp))); tmp$pop2 <- c(rep('ATL',nrow(tmp)))}
  else if (file == 13) {tmp$pop <- c(rep('CAM-4',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 14) {tmp$pop <- c(rep('CAM-5',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 15) {tmp$pop <- c(rep('CAM-6',nrow(tmp))); tmp$pop2 <- c(rep('CAM',nrow(tmp)))}
  else if (file == 16) {tmp$pop <- c(rep('WAM-4',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 17) {tmp$pop <- c(rep('WAM-5',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  else if (file == 18) {tmp$pop <- c(rep('WAM-6',nrow(tmp))); tmp$pop2 <- c(rep('WAM',nrow(tmp)))}
  
  data.combined[[file]] <- tmp
}

mu <- 1.99e-09 #unit bp per generation

data.combined <- do.call(rbind, data.combined)

time <- (data.combined$left_time_boundary/mu)
data.combined$time <- time
pop.size <- (1/data.combined$lambda)/(2*mu)
data.combined$popsize <- pop.size

ggplot(data=subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 1e07),], aes(x=(time), y=log10(popsize), fill=pop, col=pop))+
  #geom_line(linewidth = 1.5, alpha=0.75) +
  geom_step(subset(data.combined, data.combined$pop2 != 'AM')[which(subset(data.combined, data.combined$pop2 != 'AM')$time <= 1e07),], mapping=aes(x=time, y=log10(popsize)), direction = 'hv', linetype=1, linewidth=1.5)+
  scale_color_manual(values = c("#723957","#723957","#723957","#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","#A2662C","#A2662C", "#A2662C"))+
  #scale_color_manual(values = c("#723957","#4E87B8","#A2662C","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.5)+
  xlab('No. generations') + ylab('Log. Effective population size') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
  #facet_zoom(xlim=c(0,1e07))
  facet_zoom(xlim = c(0,1.5e06))
#facet_zoom(xlim = c(1.25e06,5e06))
#facet_zoom(xlim = c(5e06,10e06))
#facet_zoom(xlim = c(0,1e06))
#xlim(-1,8)

##############################
### Cross-Population analyses:
##############################
## WAM-CAM
wam.cam <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/CrossPopulation/', pattern = 'combined_WAM-CAM')
wam.cam.combined <- list()
for ( file in 1:length(wam.cam)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/CrossPopulation/',wam.cam[[file]]), header = T)
  if (file == 1) { tmp$set <- c(rep('1',nrow(tmp)))}
  else if (file == 2) {tmp$set <- c(rep('2',nrow(tmp)))}
  else if (file == 3) {tmp$set <- c(rep('3',nrow(tmp)))}
  wam.cam.combined[[file]] <- tmp
}

wam.cam.combined <- do.call(rbind, wam.cam.combined)

time <- (wam.cam.combined$left_time_boundary/mu)
wam.cam.combined$time <- time

ggplot(data=wam.cam.combined[which(wam.cam.combined$time <= 3e07),], aes(x=(time), y=((lambda_01)/((lambda_00+lambda_11)/2)),fill=set, col=set))+
  geom_step(wam.cam.combined[which(wam.cam.combined$time <= 3e07),], mapping=aes(x=(time), y=((lambda_01)/((lambda_00+lambda_11)/2))), direction = 'hv', linetype=1, linewidth=1.5)+
  #scale_color_manual(values = c("#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  scale_color_manual(values = c("gray40", "gray60","gray80"))+
  geom_hline(yintercept = 1.00, col='black', linetype='dashed')+
  geom_hline(yintercept = 0.5, col='darkred', linetype='dashed')+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.5, bins = 100)+
  xlab('No. generations') + ylab('Relative Cross-Coalescence Rate') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
  facet_zoom(xlim = c(0,2e06))
  #facet_zoom(xlim = c(3.3e05,5e05))


## ==> for WAM-CAM: start of separation between 90,000 and 500,000 generations ago 

## AM-ATL:
am.atl <- list.files('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/CrossPopulation/', pattern = 'combined_AM-ATL')
am.atl.combined <- list()
for ( file in 1:length(am.atl)){
  tmp <- read.table(paste0('~/Dropbox/Projects/Lbraz_genome_paper/Lbraziliensis/3_LBRA/2_Popgen/2_Popgen_with_PER186/10_G-PhoCS/4_MSMC2/3_msmc2_output/3_AllSNPs_EM500/CrossPopulation/',am.atl[[file]]), header = T)
  if (file == 1) { tmp$set <- c(rep('1',nrow(tmp)))}
  else if (file == 2) {tmp$set <- c(rep('2',nrow(tmp)))}
  else if (file == 3) {tmp$set <- c(rep('3',nrow(tmp)))}
  am.atl.combined[[file]] <- tmp
}

am.atl.combined <- do.call(rbind, am.atl.combined)

time <- (am.atl.combined$left_time_boundary/mu)
am.atl.combined$time <- time

ggplot(data=am.atl.combined[which(wam.cam.combined$time <= 3e07),], aes(x=(time), y=((lambda_01)/((lambda_00+lambda_11)/2)),fill=set, col=set))+
  geom_step(am.atl.combined[which(wam.cam.combined$time <= 3e07),], mapping=aes(x=(time), y=((lambda_01)/((lambda_00+lambda_11)/2))), direction = 'hv', linetype=1, linewidth=1.5)+
  #scale_color_manual(values = c("#723957","#723957","#723957","#4E87B8","#4E87B8","#4E87B8","#A2662C","#A2662C", "#A2662C","gray50","gray50","gray50"))+
  scale_color_manual(values = c("gray40", "gray60","gray80"))+
  geom_hline(yintercept = 1.00, col='black', linetype='dashed')+
  geom_hline(yintercept = 0.5, col='darkred', linetype='dashed')+
  #stat_summary_bin(fun.data='mean_cl_boot', size=0.5, bins = 100)+
  xlab('No. generations') + ylab('Relative Cross-Coalescence Rate') + 
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0, hjust=0.5)) +
  facet_zoom(xlim = c(0,6e06))
#facet_zoom(xlim = c(3.3e05,5e05))


