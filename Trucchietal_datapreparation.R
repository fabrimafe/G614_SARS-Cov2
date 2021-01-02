#R code to replicate results in Trucchi et al.,2021, MBE, Population dynamics and structural effects at short and long range support the hypothesis of the selective advantage of the G614 SARS-Cov2 spike variant
#R versions 4.0.3

######################
## DATA PREPARATION ##
######################
## In this script data are parsed and filtered. For each filter covariance matrices are calculated and saved in .RData that can be load when running the different models, which use different matrices and filters as indicated in Trucchietal_models.R

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

############## align and extract polymorphisms ##############

#These steps can be avoided by dowloading directly all.aln.varfreqmin0.005.txt from https://doi.org/10.6084/m9.figshare.13498428.v1

#Alignment to the Wuhan reference sequence with mafft v7.471
system("mafft --keeplength --addfragments allSeq.fa ref_wuhan.fa > allSeq.aln.fa")

#Convert the fasta interleaved alignment to nexus with Aliview 1.26 (GUI). The file is allSeq.aln.nex and can be found on https://doi.org/10.6084/m9.figshare.13493178.v1

#Selection of the variants based on their min freq in the alignment
system("python3 variants_in_alignment.py -i allSeq.aln.nex -d 2020-09-18 -f 0.005")

#the output is a file called all.aln.varfreqmin0.005.txt, which can be downloaded from https://doi.org/10.6084/m9.figshare.13498428.v1

############## import data ##############
seqall.metadata<-read.table("all.aln.varfreqmin0.005.txt",header=TRUE,stringsAsFactors=FALSE)
seqall.metadata$state<-seqall.metadata$country
seqall.metadata$state[seqall.metadata$country=="China"]<-paste0("China-",seqall.metadata$subcountry[seqall.metadata$country=="China"])
seqall.metadata$state[seqall.metadata$country=="USA"]<-paste0("USA-",seqall.metadata$subcountry[seqall.metadata$country=="USA"])
max(seqall.metadata$daysbefore1809) #269 -> data already in December -> starting data is 1st December 2019 -> 6*31+29+2*30+18 #293
seqall.metadata$day<-293-seqall.metadata$daysbefore1809

## Function to create fasta containing only informative sites. This allows you to download only the polymorphisms file all.aln.varfreqmin0.005.txt.
## If the original fasta files are available, steps including the function table2fa can be skipped.
pos_cols<-grep("^X",names(seqall.metadata))
table2fa<-function(mydata,pos_cols,prefix="X",fileoutput="output.fa")
{
    file.create(fileoutput)
    fileConn<-file(fileoutput)
    xpos<-names(mydata)[pos_cols]
    xpos<-as.numeric(gsub("X","",xpos))
    s0<-rep("n",xpos[length(xpos)]+1)
    st_ar<-c()
    for (i in 1:nrow(mydata))
        {        
        print(i)
        st<-s0; 
        st[xpos]<-unlist(mydata[i,pos_cols]); 
        st<-paste0(st,collapse="");
        st_ar<-c(st_ar,paste(">",mydata$ID[i]))
        st_ar<-c(st_ar,st)
        }
    writeLines(st_ar,con=fileConn);
    close(fileConn)
}

ordered_table<-function(x,xallele) { if (xallele %in% names(table(x))) table(x)[[xallele]] else return(0) } 

#### load function to compute sample covariance for allele frequency, similarly to the one used by Treemix. The correction for branch length is adjusted for haploid genomes, but note that it makes almost no effect ~10^-9.
## Note that sometimes sample covariance matrices ( as our case, if many correlations) are not positive definite but semipositive definite, and in these cases Cholesky decomposition fails (so brms!). There are two solutions. To solve this the function implement the common practice to add a tiny number to the diagonal that does not affect the estimates. See also https://www.value-at-risk.net/non-positive-definite-covariance-matrices/

covariance_alleles_treemix<-function(xdata,pos_cols,withcorrectiontreemix=TRUE)
{
temp_hi<-rep(0,length(unique(xdata$state)))
for (itt in 1:length(pos_cols))
#it is the 66 that gives problem!!
#for (itt in 1:66)
    {
    xxdata<-xdata[,pos_cols[itt]]
    xxdata[xxdata=="r" | xxdata=="?"]<-NA
    xxdata<-tibble(allele=xxdata,state=xdata$state)
    major_allele<-names(table(xxdata$allele))[table(xxdata$allele)==max(table(xxdata$allele))][1]
    xxdata<- xxdata %>% group_by(state) %>% summarise(af=ordered_table(allele,major_allele)/sum(table(allele)),n=sum(table(allele))) 
    temp_mat<-sapply(xxdata$af,function(x) (xxdata$af-mean(xxdata$af,na.rm=T))*(x-mean(xxdata$af,na.rm=T)))
    temp_hi<-temp_hi+xxdata$af*(1-xxdata$af)*xxdata$n/(xxdata$n-1)    
    #covariance matrix
    # --- to explore if state is weird ---
    #    print(c(itt,mean(xxdata$af,na.rm=T),xxdata$af[xxdata$state=="USA-WA"],xxdata$af[xxdata$state=="USA-CA"]))
    # --- ---------------------------- ---
    if ( itt == 1 ) { covariance_mat<-temp_mat } else { covariance_mat<-covariance_mat+temp_mat }
    }
    temp_hi<-temp_hi/length(pos_cols)/2/table(xdata$state) #Bi according to treemix paper S1 (for haploid)    
    for (i in 1:nrow(covariance_mat))
        {
        for (j in 1:nrow(covariance_mat))
            {
            if (withcorrectiontreemix)
                {
                covariance_mat[i,j]<-covariance_mat[i,j]+
                temp_hi[i]/length(unique(xdata$state))+
                temp_hi[j]/length(unique(xdata$state))-
                sum(temp_hi)/length(unique(xdata$state))^2;
                if (i==j) {
                    covariance_mat[i,j]<-covariance_mat[i,j]-temp_hi[i]
                    }
                }
            if (covariance_mat[i,j]!=covariance_mat[j,i]){ covariance_mat[j,i]<-covariance_mat[i,j] }
            }
        }
    covariance_mat<-covariance_mat/length(pos_cols)
    covariance_mat<-as.matrix(covariance_mat)
    rownames(covariance_mat)<-xxdata$state
    colnames(covariance_mat)<-xxdata$state
    covariance_mat
}



##############################################
##    CREATE FILTERS FOR DIFFERENT MODELS   ##
##############################################

##MINIMUM FILTER (all sequences)

filter0<-names(table(seqall.metadata$country))[table(seqall.metadata$country)>30]
filter0<-names(table(seqall.metadata$state))[table(seqall.metadata$state)>30]
filter0<-filter0[filter0!="mink"] 
filter0<-filter0[filter0!="US-UN"]
xdata<-seqall.metadata[seqall.metadata$state %in% filter0,]

covariance_mat<-covariance_alleles_treemix(xdata,pos_cols)
corr_matrix<-cov2cor(covariance_mat)

##this is the list of geographical regions
write.table(filter0,file="filter0.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

## to be run but slow
## table2fa(seqall.metadata,pos_cols,fileoutput="seqall.fa")

seqall = read.dna('seqall.fa', format='fasta')

 pol = DNAbin2genind(seqall, pop = as.character(seqall.metadata$state))

 pol.sel = repool(seppop(pol)$Australia, seppop(pol)$Austria,seppop(pol)$Belgium,
 seppop(pol)[["Canada"]],
 seppop(pol)[["China-Beijing"]],
 seppop(pol)[["China-Guangdong"]],
 seppop(pol)[["China-Guangzhou"]],
 seppop(pol)[["China-Hangzhou"]],
 seppop(pol)[["China-Shanghai"]],
 seppop(pol)[["China-Sichuan"]],
 seppop(pol)[["China-Wuhan"]],
 seppop(pol)[["CostaRica"]],
 seppop(pol)[["CzechRepublic"]],
 seppop(pol)[["Denmark"]],
 seppop(pol)[["England"]],
 seppop(pol)[["Finland"]],
 seppop(pol)[["France"]],
 seppop(pol)[["Germany"]],
 seppop(pol)[["Greece"]],
 seppop(pol)[["Hungary"]],
 seppop(pol)[["Iceland"]],
 seppop(pol)[["Ireland"]],
 seppop(pol)[["Italy"]],
 seppop(pol)[["Latvia"]],
 seppop(pol)[["Lithuania"]],
 seppop(pol)[["Luxembourg"]],
 seppop(pol)[["Mexico"]],
 seppop(pol)[["Netherlands"]],
 seppop(pol)[["NorthernIreland"]],
 seppop(pol)[["NorthMacedonia"]],
 seppop(pol)[["Norway"]],
 seppop(pol)[["Poland"]],
 seppop(pol)[["Portugal"]],
 seppop(pol)[["Romania"]],
 seppop(pol)[["Russia"]],
 seppop(pol)[["Scotland"]],
 seppop(pol)[["Serbia"]],
 seppop(pol)[["Spain"]],
 seppop(pol)[["Sweden"]],
 seppop(pol)[["Switzerland"]],
 seppop(pol)[["Turkey"]],
 seppop(pol)[["Ukraine"]],
 seppop(pol)[["USA-AK"]],
 seppop(pol)[["USA-AR"]],
 seppop(pol)[["USA-AZ"]],
 seppop(pol)[["USA-CA"]],
 seppop(pol)[["USA-CO"]],
 seppop(pol)[["USA-CT"]],
 seppop(pol)[["USA-DC"]],
 seppop(pol)[["USA-FL"]],
 seppop(pol)[["USA-GA"]],
 seppop(pol)[["USA-IA"]],
 seppop(pol)[["USA-ID"]],
 seppop(pol)[["USA-IL"]],
 seppop(pol)[["USA-LA"]],
 seppop(pol)[["USA-MA"]],
 seppop(pol)[["USA-MD"]],
 seppop(pol)[["USA-ME"]],
 seppop(pol)[["USA-MI"]],
 seppop(pol)[["USA-MN"]],
 seppop(pol)[["USA-NE"]],
 seppop(pol)[["USA-NJ"]],
 seppop(pol)[["USA-NM"]],
 seppop(pol)[["USA-NV"]],
 seppop(pol)[["USA-NY"]],
 seppop(pol)[["USA-OR"]],
 seppop(pol)[["USA-PA"]],
 seppop(pol)[["USA-SC"]],
 seppop(pol)[["USA-TX"]],
 seppop(pol)[["USA-UN"]],
 seppop(pol)[["USA-UT"]],
 seppop(pol)[["USA-VA"]],
 seppop(pol)[["USA-WA"]],
 seppop(pol)[["USA-WI"]],
 seppop(pol)[["USA-WY"]],
 seppop(pol)[["Wales"]])

mathierfstat<-genind2hierfstat(pol.sel)
matFst <- pairwise.WCfst(mathierfstat,diploid=FALSE)

## heatmaps of Fst and covariance matrices
heatmap(matFst,cexCol=0.5,,cexRow=0.5)
heatmap(1-corr_matrix,cexCol=0.5,,cexRow=0.5)

pol.genpop = genind2genpop(pol.sel)
Dall = dist.genpop(pol.genpop, method = 1)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="complete")
phylotree = as.phylo(hc)

plot(phylotree, type="fan")
dev.off()
write.nexus(phylotree, file="seqall_nofilter_nexus.nex")
phylotree<-read.nexus(file="seqall_nofilter_nexus.nex")

Dall<-as.dist(matFst)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="ward.D2")
phylotreeFst = as.phylo(hc)
A <- ape::vcv.phylo(phylotreeFst)
is.positive.definite(A)

plot(phylotreeFst, type="fan")
dev.off()
write.nexus(phylotreeFst, file="seqall_nofilter_nexus_Fst.nex")
phylotreeFst<-read.nexus(file="seqall_nofilter_nexus_Fst.nex")

corrplot(corr_matrix,method="shade",order="AOE",tl.cex=0.3,shade.col=NA,tl.col="black",cl.pos="n")

plot(phylotree.cov, type="fan",use.edge.length=FALSE)
dev.off()
write.nexus(phylotree.cov, file="seqall_nofilter_nexus_cov.nex")
phylotree.cov<-read.nexus(file="seqall_nofilter_nexus_cov.nex")

#construct the inverse of A matrix of the phylogeny
Ainv.phyloall<-inverseA(phylotree,nodes="TIPS")$Ainv

seqall_model<-seqall.metadata[seqall.metadata$state %in% filter0,]
seqall_model$tips<-seqall_model$state

#polymorphism 85 is G614. However, the same analyses can be performed for all variants
res_l<-list()
i=85 
seqall_model$type<-as.character(seqall_model[[i]])
seqall_model$type[seqall_model$type!="a" & seqall_model$type!="c" & seqall_model$type!="g" & seqall_model$type!="t"]<-NA
var1<-"g"
seqall_model$type[seqall_model$type==var1]<-1
seqall_model$type[seqall_model$type!=var1 & seqall_model$type!=1]<-0
seqall_model$type <- as.integer(seqall_model$type)

save.image(file="workspace0.RData")
load("workspace0.RData")

heatmap(covariance_mat,cexCol=0.5,cexRow=0.5,scale="none")
A <- ape::vcv.phylo(phylotreeFst)
heatmap(A,cexCol=0.5,cexRow=0.5,scale="none")


#######################################
############### FILTER LOCKDOWN #######
#######################################

filter_lockdown1<-seqall.metadata$day>=(31+31+29+31) & seqall.metadata$day<(31+31+29+31+30)
states_in_filter_lockdown1<-table(seqall.metadata$state[filter_lockdown1])[table(seqall.metadata$state[filter_lockdown1])>30]
states_in_filter_lockdown1<-names(states_in_filter_lockdown1)
xdata.april<-seqall.metadata[seqall.metadata$state %in% states_in_filter_lockdown1 & filter_lockdown1,]

pos_cols<-grep("^X",names(xdata.april))
table2fa(xdata.april,pos_cols,fileoutput="seqall_april.fa")


covariance_mat<-covariance_alleles_treemix(xdata.april,pos_cols)
corr_matrix<-cov2cor(covariance_mat)

seqall = read.dna('seqall_april.fa', format='fasta')

pol = DNAbin2genind(seqall, pop = as.character(xdata.april$state))

 pol.sel = repool(seppop(pol)$Australia, seppop(pol)$Austria,seppop(pol)$Belgium,
 seppop(pol)[["Canada"]],
 seppop(pol)[["Denmark"]],
 seppop(pol)[["England"]],
 seppop(pol)[["Finland"]],
 seppop(pol)[["France"]],
 seppop(pol)[["Germany"]],
 seppop(pol)[["Hungary"]],
 seppop(pol)[["Italy"]],
 seppop(pol)[["Lithuania"]],
 seppop(pol)[["Luxembourg"]],
 seppop(pol)[["Netherlands"]],
 seppop(pol)[["NorthernIreland"]],
 seppop(pol)[["Poland"]],
 seppop(pol)[["Portugal"]],
 seppop(pol)[["Romania"]],
 seppop(pol)[["Russia"]],
 seppop(pol)[["Scotland"]],
 seppop(pol)[["Spain"]],
 seppop(pol)[["Sweden"]],
 seppop(pol)[["Switzerland"]],
 seppop(pol)[["Turkey"]],
 seppop(pol)[["USA-CA"]],
 seppop(pol)[["USA-CO"]],
 seppop(pol)[["USA-CT"]],
 seppop(pol)[["USA-FL"]],
 seppop(pol)[["USA-IL"]],
 seppop(pol)[["USA-LA"]],
 seppop(pol)[["USA-MA"]],
 seppop(pol)[["USA-MD"]],
 seppop(pol)[["USA-MI"]],
 seppop(pol)[["USA-MN"]],
 seppop(pol)[["USA-NM"]],
 seppop(pol)[["USA-NV"]],
 seppop(pol)[["USA-NY"]],
 seppop(pol)[["USA-OR"]],
 seppop(pol)[["USA-TX"]],
 seppop(pol)[["USA-UN"]],
 seppop(pol)[["USA-UT"]],
 seppop(pol)[["USA-VA"]],
 seppop(pol)[["USA-WA"]],
 seppop(pol)[["USA-WI"]],
 seppop(pol)[["Wales"]])


pol.genpop = genind2genpop(pol.sel)
Dall = dist.genpop(pol.genpop, method = 1)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="complete")
phylotree = as.phylo(hc)
plot(phylotree, type="fan")
write.nexus(phylotree, file="seqall_lockdown1_nexus.nex")
phylotree<-read.nexus(file="seqall_lockdown1_nexus.nex")
#construct the inverse of A matrix of the phylogeny
Ainv.phyloall<-inverseA(phylotree,nodes="TIPS")$Ainv

seqall_model<-xdata.april
seqall_model$tips<-xdata.april$state

res_l<-list()
i=85
seqall_model$type<-as.character(seqall_model[[i]])
seqall_model$type[seqall_model$type!="a" & seqall_model$type!="c" & seqall_model$type!="g" & seqall_model$type!="t"]<-NA
var1<-"g"
seqall_model$type[seqall_model$type==var1]<-1
seqall_model$type[seqall_model$type!=var1 & seqall_model$type!=1]<-0
seqall_model$type <- as.integer(seqall_model$type)

save.image(file="workspace1.RData")
load(file="workspace1.RData")

heatmap(covariance_mat,cexCol=0.5,cexRow=0.5,scale="none")
phylo <- ape::read.nexus("seqall_nofilter_nexus_Fst.nex")
A <- ape::vcv.phylo(phylo)
heatmap(A,cexCol=0.5,cexRow=0.5,scale="none")

#######################################
############## FILTER G ###############
#######################################
## filter including only sequences carrying the G allele. This is done to estimate migration exclusively based on the relationships
## between G and to avoid biases introduced by changes in allele frequencies between the source and recipient populations.
filter0<-names(table(seqall.metadata$country))[table(seqall.metadata$country)>30]
filter0<-names(table(seqall.metadata$state))[table(seqall.metadata$state)>30]
filter0<-filter0[filter0!="mink"] #removed because not real states
filter0<-filter0[filter0!="USA-UN"] #removed because not real states
filterG<-seqall.metadata[,85]=="g"
filterA<-seqall.metadata[,85]=="a"
states_in_filter_G<-table(seqall.metadata$state[filterG])[table(seqall.metadata$state[filterG])>3]
states_in_filter_A<-table(seqall.metadata$state[filterA])[table(seqall.metadata$state[filterA])>0]

names(states_in_filter_G) %in% filter0 #not totally overlapping
filter0 %in% names(states_in_filter_G) #not totally overlapping
filter0 %in% names(states_in_filter_A) #not totally overlapping but pretty OK (only 5 not)
states_in_filter_G<-filter0[filter0 %in% names(states_in_filter_G) & filter0 %in% names(states_in_filter_A)]

xdata.G<-seqall.metadata[seqall.metadata$state %in% states_in_filter_G & filterG,]

pos_cols<-grep("^X",names(xdata.G))
table2fa(xdata.G,pos_cols,fileoutput="seqG.fa")


covariance_mat<-covariance_alleles_treemix(xdata.G,pos_cols)
corr_matrix<-cov2cor(covariance_mat)
corrplot(corr_matrix,method="shade",order="AOE",tl.cex=0.3,shade.col=NA,tl.col="black",cl.pos="n")


seqall = read.dna('seqG.fa', format='fasta')

seqall$state
pol = DNAbin2genind(seqall, pop = as.character(xdata.G$state))

 pol.sel = repool(seppop(pol)$Australia, seppop(pol)$Austria,seppop(pol)$Belgium,
 seppop(pol)[["Canada"]],
 seppop(pol)[["China-Beijing"]],
 seppop(pol)[["China-Hangzhou"]],
 seppop(pol)[["CostaRica"]],
 seppop(pol)[["CzechRepublic"]],
 seppop(pol)[["Denmark"]],
 seppop(pol)[["England"]],
 seppop(pol)[["Finland"]],
 seppop(pol)[["France"]],
 seppop(pol)[["Germany"]],
 seppop(pol)[["Greece"]],
 seppop(pol)[["Hungary"]],
 seppop(pol)[["Iceland"]],
 seppop(pol)[["Ireland"]],
 seppop(pol)[["Italy"]],
 seppop(pol)[["Latvia"]],
 seppop(pol)[["Lithuania"]],
 seppop(pol)[["Luxembourg"]],
 seppop(pol)[["Mexico"]],
 seppop(pol)[["Netherlands"]],
 seppop(pol)[["NorthernIreland"]],
 seppop(pol)[["Norway"]],
 seppop(pol)[["Poland"]],
 seppop(pol)[["Portugal"]],
 seppop(pol)[["Russia"]],
 seppop(pol)[["Scotland"]],
 seppop(pol)[["Spain"]],
 seppop(pol)[["Sweden"]],
 seppop(pol)[["Switzerland"]],
 seppop(pol)[["Turkey"]],
 seppop(pol)[["Ukraine"]],
 seppop(pol)[["USA-AK"]],
 seppop(pol)[["USA-AZ"]],
 seppop(pol)[["USA-CA"]],
 seppop(pol)[["USA-CO"]],
 seppop(pol)[["USA-CT"]],
 seppop(pol)[["USA-DC"]],
 seppop(pol)[["USA-FL"]],
 seppop(pol)[["USA-GA"]],
 seppop(pol)[["USA-ID"]],
 seppop(pol)[["USA-IL"]],
 seppop(pol)[["USA-LA"]],
 seppop(pol)[["USA-MA"]],
 seppop(pol)[["USA-MD"]],
 seppop(pol)[["USA-ME"]],
 seppop(pol)[["USA-MI"]],
 seppop(pol)[["USA-MN"]],
 seppop(pol)[["USA-NE"]],
 seppop(pol)[["USA-NJ"]],
 seppop(pol)[["USA-NM"]],
 seppop(pol)[["USA-NV"]],
 seppop(pol)[["USA-NY"]],
 seppop(pol)[["USA-OR"]],
 seppop(pol)[["USA-PA"]],
 seppop(pol)[["USA-SC"]],
 seppop(pol)[["USA-TX"]],
 seppop(pol)[["USA-UT"]],
 seppop(pol)[["USA-VA"]],
 seppop(pol)[["USA-WA"]],
 seppop(pol)[["USA-WI"]],
 seppop(pol)[["USA-WY"]],
 seppop(pol)[["Wales"]])
 
pol.genpop = genind2genpop(pol.sel)
Dall = dist.genpop(pol.genpop, method = 1)
mathierfstat<-genind2hierfstat(pol.sel)
matFst <- pairwise.WCfst(mathierfstat,diploid=FALSE)

heatmap(matFst,cexCol=0.5,,cexRow=0.5)
heatmap(1-corr_matrix,cexCol=0.5,,cexRow=0.5)
mypca<-pcoa(Dall, correction="none", rn=NULL)
biplot(mypca,cex=0.1)

temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="complete")
phylotree = as.phylo(hc)
plot(phylotree, type="fan")
write.nexus(phylotree, file="seqall_G_nexus.nex")
phylotree<-read.nexus(file="seqall_G_nexus.nex")
#construct the inverse of A matrix of the phylogeny
Ainv.phyloall<-inverseA(phylotree,nodes="TIPS")$Ainv

Dall<-as.dist(matFst)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="ward.D2")
phylotreeFst = as.phylo(hc)
A <- ape::vcv.phylo(phylotreeFst)
is.positive.definite(A)
plot(phylotreeFst, type="fan")
write.nexus(phylotreeFst, file="seqall_G_nexus_Fst.nex")
phylotreeFst<-read.nexus(file="seqall_G_nexus_Fst.nex")

#library(hierfstat)
#matFst <- pairwise.fst(pol.sel,res.type="matrix")

seqall_model<-seqall.metadata[seqall.metadata$state %in% states_in_filter_G,]
seqall_model$tips<-seqall_model$state


res_l<-list()
i=85
seqall_model$type<-as.character(seqall_model[[i]])
seqall_model$type[seqall_model$type!="a" & seqall_model$type!="c" & seqall_model$type!="g" & seqall_model$type!="t"]<-NA
var1<-"g"
seqall_model$type[seqall_model$type==var1]<-1
seqall_model$type[seqall_model$type!=var1 & seqall_model$type!=1]<-0
seqall_model$type <- as.integer(seqall_model$type)

heatmap(covariance_mat,cexCol=0.5,cexRow=0.5,scale="none")
A <- ape::vcv.phylo(phylotreeFst)
heatmap(A,cexCol=0.5,cexRow=0.5,scale="none")

save.image(file="workspace2.RData")
load("workspace2.RData")

####################################
##COMBINING FILTER LOCKDOWN + G ####
####################################

filter_lockdown1<-seqall.metadata$day>=(31+31+29+31) & seqall.metadata$day<(31+31+29+31+30)
states_in_filter_lockdown1<-table(seqall.metadata$state[filter_lockdown1])[table(seqall.metadata$state[filter_lockdown1])>30]
states_in_filter_lockdown1<-names(states_in_filter_lockdown1)
xdata.april<-seqall.metadata[seqall.metadata$state %in% states_in_filter_lockdown1 & filter_lockdown1,]
pos_cols<-grep("^X",names(xdata.april))

filter0<-names(table(xdata.april$country))[table(xdata.april$country)>30]
filter0<-names(table(xdata.april$state))[table(xdata.april$state)>30]
filter0<-filter0[filter0!="mink"] 
filter0<-filter0[filter0!="USA-UN"]
filterG<-xdata.april[,85]=="g"
filterA<-xdata.april[,85]=="a"
states_in_filter_G<-table(xdata.april$state[filterG])[table(xdata.april$state[filterG])>3]
states_in_filter_A<-table(xdata.april$state[filterA])[table(xdata.april$state[filterA])>0]

states_in_filter_G<-filter0[filter0 %in% names(states_in_filter_G) & filter0 %in% names(states_in_filter_A)]
xdata.G<-xdata.april[xdata.april$state %in% states_in_filter_G & filterG,]

pos_cols<-grep("^X",names(xdata.G))
table2fa(xdata.G,pos_cols,fileoutput="seqG.april.fa")

covariance_mat<-covariance_alleles_treemix(xdata.G,pos_cols)
corr_matrix<-cov2cor(covariance_mat)

seqall = read.dna('seqG.april.fa', format='fasta')
sort(unique(xdata.G$state))

pol = DNAbin2genind(seqall, pop = as.character(xdata.G$state))

 pol.sel = repool(seppop(pol)$Australia, seppop(pol)$Austria,seppop(pol)$Belgium,
 seppop(pol)[["Canada"]],
 seppop(pol)[["Denmark"]],
 seppop(pol)[["England"]],
 seppop(pol)[["Finland"]],
 seppop(pol)[["Germany"]],
 seppop(pol)[["Lithuania"]],
 seppop(pol)[["Luxembourg"]],
 seppop(pol)[["Netherlands"]],
 seppop(pol)[["NorthernIreland"]],
 seppop(pol)[["Poland"]],
 seppop(pol)[["Portugal"]],
 seppop(pol)[["Scotland"]],
 seppop(pol)[["Spain"]],
 seppop(pol)[["Sweden"]],
 seppop(pol)[["Switzerland"]],
 seppop(pol)[["Turkey"]],
 seppop(pol)[["USA-CA"]],
 seppop(pol)[["USA-CT"]],
 seppop(pol)[["USA-FL"]],
 seppop(pol)[["USA-IL"]],
 seppop(pol)[["USA-LA"]],
 seppop(pol)[["USA-MA"]],
 seppop(pol)[["USA-MD"]],
 seppop(pol)[["USA-MI"]],
 seppop(pol)[["USA-MN"]],
 seppop(pol)[["USA-NM"]],
 seppop(pol)[["USA-NV"]],
 seppop(pol)[["USA-NY"]],
 seppop(pol)[["USA-OR"]],
 seppop(pol)[["USA-TX"]],
 seppop(pol)[["USA-UT"]],
 seppop(pol)[["USA-VA"]],
 seppop(pol)[["USA-WA"]],
 seppop(pol)[["USA-WI"]],
 seppop(pol)[["Wales"]])
 
pol.genpop = genind2genpop(pol.sel)

mathierfstat<-genind2hierfstat(pol.sel)
matFst <- pairwise.WCfst(mathierfstat,diploid=FALSE)

heatmap(matFst,cexCol=0.5,,cexRow=0.5)
heatmap(1-corr_matrix,cexCol=0.5,,cexRow=0.5)

Dall = dist.genpop(pol.genpop, method = 1)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="complete")
phylotree = as.phylo(hc)
plot(phylotree, type="fan")
write.nexus(phylotree, file="seqall_aprilG_nexus.nex")
phylotree<-read.nexus(file="seqall_aprilG_nexus.nex")
#construct the inverse of A matrix of the phylogeny
Ainv.phyloall<-inverseA(phylotree,nodes="TIPS")$Ainv

Dall<-as.dist(matFst)
temp <- as.data.frame(as.matrix(Dall))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)
tre <- upgma(Dall)
hc = hclust(Dall, method="ward.D2")
phylotreeFst = as.phylo(hc)
A <- ape::vcv.phylo(phylotreeFst)
is.positive.definite(A)
plot(phylotreeFst, type="fan")
write.nexus(phylotreeFst, file="seqall_aprilG_nexus_Fst.nex")
phylotreeFst<-read.nexus(file="seqall_aprilG_nexus_Fst.nex")

#library(hierfstat)
#matFst <- pairwise.fst(pol.sel,res.type="matrix")

seqall_model<-xdata.april[xdata.april$state %in% states_in_filter_G,]
seqall_model$tips<-seqall_model$state


res_l<-list()
i=85
seqall_model$type<-as.character(seqall_model[[i]])
seqall_model$type[seqall_model$type!="a" & seqall_model$type!="c" & seqall_model$type!="g" & seqall_model$type!="t"]<-NA
var1<-"g"
seqall_model$type[seqall_model$type==var1]<-1
seqall_model$type[seqall_model$type!=var1 & seqall_model$type!=1]<-0
seqall_model$type <- as.integer(seqall_model$type)

rm(seqall)
save.image(file="workspace3.RData")
load("workspace3.RData")

heatmap(covariance_mat,cexCol=0.5,cexRow=0.5,scale="none")

A <- ape::vcv.phylo(phylotreeFst)
heatmap(A,cexCol=0.5,cexRow=0.5,scale="none")
