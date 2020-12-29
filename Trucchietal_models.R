#R code to replicate results in Trucchi et al.,2021, MBE
#R versions 4.0.3

############## load libraries ##############

library(ape)
library(adegenet)
library(phangorn)
library(MCMCglmm)
library(brms)
library(ggplot2)
library(tidyverse)
library(hierfstat)
library(corrplot)

############## no lockdown filter - sample allele frequency covariance on all sequences ######################

load("workspace0.RData")

m4 <- brm(
  type ~ day + (1+day|tips)+(1+day|state), data = seqall_model,
  family = bernoulli(), cov = list(tips = covariance_mat),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m4,file="m4.RData")

############## no lockdown filter - sample allele frequency covariance on genomes carrying G ######################

load("workspace2.RData")

covariance_mat<-diag(0.00001,nrow(as.matrix(covariance_mat)))+as.matrix(covariance_mat)
m4G <- brm(
  type ~ day + (1+day|tips)+(1+day|state), data = seqall_model,
  family = bernoulli(), cov = list(tips = covariance_mat),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m4G,file="m4G.RData")

############## lockdown filter - sample allele frequency covariance on genomes carrying G ######################

load("workspace3.RData")

covariance_mat<-diag(0.0000001,nrow(as.matrix(covariance_mat)))+as.matrix(covariance_mat)
m4l1G <- brm(
  type ~ day + (1+day|tips)+(1+day|state), data = seqall_model,
  family = bernoulli(), cov = list(tips = covariance_mat),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m4l1G,file="m4l1G.RData")

############## no lockdown filter - Fst-based covariance on all sequences ######################

load("workspace0.RData")

phylo <- ape::read.nexus("seqall_nofilter_nexus_Fst.nex")
A <- ape::vcv.phylo(phylo)
m6 <- brm(
  type ~ day + (1+day|gr(tips,cov=A))+(1+day|state), data = seqall_model,
  family = bernoulli(), data2 = list(A = A),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m6,file="m6.RData")

############## no lockdown filter - Fst-based covariance on genomes carrying G ######################

load("workspace2.RData")

phylo <- ape::read.nexus("seqall_nofilter_nexus_Fst.nex")
A <- ape::vcv.phylo(phylo)
m6G <- brm(
  type ~ day + (1+day|gr(tips,cov=A))+(1+day|state), data = seqall_model,
  family = bernoulli(), data2 = list(A = A),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m6G,file="m6G.RData")

############## lockdown filter - Fst-based covariance on genomes carrying G ######################

load("workspace3.RData")

phylo <- ape::read.nexus("seqall_nofilter_nexus_Fst.nex")
A <- ape::vcv.phylo(phylo)
m6l1G <- brm(
  type ~ day + (1+day|gr(tips,cov=A))+(1+day|state), data = seqall_model,
  family = bernoulli(), data2 = list(A = A),iter=4000, warmup=2000,save_all_pars=TRUE,cores=4, control=list(adapt_delta=0.98,max_treedepth = 25)
)
save(m6l1G,file="m6l1G.RData")


############## no lockdown filter - spatial model on all sequences ######################

library(spaMM)

load("workspace0.RData")

coords<-read.table("territories_pop_centers.txt",header=T)
names(coords)[1]<-"state"
coords$region<-coords$state
coords$region[grep("USA",coords$state)]<-"USA"
coords$region[grep("China",coords$state)]<-"China"

coords$region[! (coords$region %in% c("China","USA","Australia","Canada","CostaRica","Mexico","Russia"))]<-"Europe"
coords$region[coords$region %in% c("Mexico","CostaRica")]<-"CentralAmerica"
seqall_model <- left_join(seqall_model,coords)

m1geo_spaMM_null <- fitme(type~1+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1| lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1geo_spaMM_null,file="m1geo_spaMM_null.RData")

m1geo_spaMM <- fitme(type~1+day+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1 | lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1geo_spaMM,file="m1geo_spaMM.RData")

############## no lockdown filter - spatial model on genomes carrying G ######################

load("workspace2.RData")

coords<-read.table("territories_pop_centers.txt",header=T)
names(coords)[1]<-"state"
coords$region<-coords$state
coords$region[grep("USA",coords$state)]<-"USA"
coords$region[grep("China",coords$state)]<-"China"

coords$region[! (coords$region %in% c("China","USA","Australia","Canada","CostaRica","Mexico","Russia"))]<-"Europe"
coords$region[coords$region %in% c("Mexico","CostaRica")]<-"CentralAmerica"
seqall_model <- left_join(seqall_model,coords)

m1Ggeo_spaMM_null <- fitme(type~1+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1| lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1Ggeo_spaMM_null,file="m1Ggeo_spaMM_null.RData")

m1Ggeo_spaMM <- fitme(type~1+day+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1 | lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1Ggeo_spaMM,file="m1Ggeo_spaMM.RData")

############## lockdown filter - spatial model on genomes carrying G ######################

load("workspace3.RData")

coords<-read.table("territories_pop_centers.txt",header=T)
names(coords)[1]<-"state"
coords$region<-coords$state
coords$region[grep("USA",coords$state)]<-"USA"
coords$region[grep("China",coords$state)]<-"China"

coords$region[! (coords$region %in% c("China","USA","Australia","Canada","CostaRica","Mexico","Russia"))]<-"Europe"
coords$region[coords$region %in% c("Mexico","CostaRica")]<-"CentralAmerica"
seqall_model <- left_join(seqall_model,coords)

m1l1geo_spaMM <- fitme(type~1+day+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1 | lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1l1geo_spaMM,file="m1l1geo_spaMM.RData")

m1l1geo_spaMM_null <- fitme(type~1+(1+day|state)+(1+day|region)+Matern(0+day | lon + lat)+Matern(1| lon + lat),
data=seqall_model,family=binomial(),method="ML")
save(m1l1geo_spaMM_null,file="m1l1geo_spaMM_null.RData")


############## source sink model ##############

load("workspace0.RData")

load("m4.RData")
library(tidyverse)

head(m4$data)

datapredict<-m4$data
datapredict$day<-62
datapredict$type<-0
datapredict<- datapredict %>% filter(!duplicated(datapredict))
predicted<-predict(m4,newdata=datapredict)
predicted<-cbind(datapredict,predicted)
predicted[order(predicted$Estimate),]

Nm<-((1-matFst)/matFst)[colnames(matFst)=="Iceland",]
predicted$Iceland_model<-sapply(names(Nm), function(x) (predicted$Estimate[predicted$state=="Iceland"]-predicted$Estimate[predicted$state==x])*Nm[names(Nm)==x])
predicted$Iceland_model[predicted$Iceland_model<0 | is.na(predicted$Iceland_model)]<-0.000001
Nm<-((1-matFst)/matFst)[colnames(matFst)=="Italy",]
predicted$Italy_model<-sapply(names(Nm), function(x) (1-predicted$Estimate[predicted$state==x])*Nm[names(Nm)==x])
predicted$Italy_model[predicted$Italy_model<0 | is.na(predicted$Italy_model)]<-0.000001
Nm<-((1-matFst)/matFst)[colnames(matFst)=="Belgium",]
predicted$Belgium_model<-sapply(names(Nm), function(x) (predicted$Estimate[predicted$state=="Belgium"]-predicted$Estimate[predicted$state==x])*Nm[names(Nm)==x])
predicted$Belgium_model[predicted$Belgium_model<0 | is.na(predicted$Belgium_model)]<-0.000001
predicted$type<-NULL
predicted$day<-NULL
predicted$tips<-NULL

seqall_model<-left_join(seqall_model,predicted,by="state")

seqall_model$Iceland_model.day<-seqall_model$Iceland_model*seqall_model$day
seqall_model$Italy_model.day<-seqall_model$Italy_model*seqall_model$day
seqall_model$Belgium_model.day<-seqall_model$Belgium_model*seqall_model$day

seqall_model2<-seqall_model
seqall_model2$Iceland_model<-seqall_model2$Iceland_model/max(seqall_model2$Iceland_model)

write.table(seqall_model,file="seqall_model_filter0_sourcesink.tab",quote=FALSE,col.names=TRUE,row.names=FALSE)


library(brms)
seqall_model<-read.table("seqall_model_filter0_sourcesink.tab",header=TRUE)
m.iceland <- brm(
  type ~ day + (1+day+Iceland_model.day|state), data = seqall_model, 
  family = bernoulli(), save_all_pars=TRUE,cores=4,iter=2000
)
save(m.iceland,file="m.iceland.RData")

m.belgium <- brm(
  type ~ day + (1+day+Belgium_model.day|state), data = seqall_model, 
  family = bernoulli(), save_all_pars=TRUE,cores=4,iter=2000
)
save(m.belgium,file="m.belgium.RData")

m.italy <- brm(
  type ~ day + (1+day+Italy_model.day|state), data = seqall_model, 
  family = bernoulli(), save_all_pars=TRUE,cores=4,iter=2000
)
save(m.italy,file="m.italy.RData")

m.iceland <- brm(
  type ~ day + (1+day+Belgium_model.day|state), data = seqall_model, 
  family = binomial("logit"), save_all_pars=TRUE,cores=4,iter=2000
)
save(m.iceland,file="m.belgium.RData")


