
### AUTHOR: Flora Jay ###


#################
### LOAD DATA ###
#################

# Replace by the proper path to file, and change the number of rows you want to load:
snp = read.table("../GdP-material/chr22.CG_54genomes_shapeit_phased.txt.haps",
                 nrows=100000)  # at first load a smaller nb of rows, eg 1000 then try witt 100000
sample.info = read.table("../GdP-material/CG_54genomes_indiv.txt",
                         header=T,stringsAsFactors = F)

head(snp)
head(sample.info)

info = snp[,1:5]
snp = snp[,-(1:5)]


#### Setting up color palettes:
upop = unique(sample.info[,"POP"])
ureg = unique(sample.info[,"REGION"])
# if possible :
require(RColorBrewer)
colpop =  brewer.pal(length(upop),"Paired")    #rainbow(length(upop))
names(colpop) = upop
sample.colpop = colpop[sample.info[,"POP"]]
colreg = brewer.pal(length(ureg),"Set1")    #rainbow(length(ureg))
names(colreg) = ureg
sample.colreg = colreg[sample.info[,"REGION"]]



######################
#### Heterozygotie ###
######################
# Un individu a 2 alleles (colonnes successives dans le tableau)

# Astuce
# les nombres impairs
odd = seq(1,ncol(snp),by=2)

# Example for one individual
column=5
het5=snp[,column] != snp[, column+1]
mean(het5)
het=NULL
## DO NOT RUN
for (indiv in odd) {
#  het= c(het, ???) # replace ??? by the proper expression
}
##

# FJ version with sapply:
het = sapply(odd, FUN=function(i) mean(snp[,i]!=snp[,i+1]) )

# Plot
plot(het,pch=19,type="n")
text(het,sample.info[odd,"POP"],col=sample.colpop[odd])

# Plot ordered het.
rank=order(het)
plot(het[rank],pch=19,type="n")
text(het[rank],sample.info[odd,"POP"][rank],col=sample.colreg[odd][rank])
legend("topleft",names(colreg),pch=19,col=colreg,ncol=5)


#############
### THETA ###
#############


### Estimateur de theta (Watterson)
theta_watterson_est = function(snp) {
  allele.counts = rowSums(snp)
  fixed = allele.counts %in% c(0,ncol(snp))  # sites that are no polymorphic (if all 0 -> sum=0 if all 1 -> sum=n)
  S = sum(!fixed) # observed number of segragating sites
  theta.wat = S / sum( 1/1:(ncol(snp)-1)) 
  return(theta.wat)
}

### Estimateur de theta (Tajima)
# loop version is easier to understand at first:
theta_tajima_est_loop = function(snp) {
  dij = NULL
  for (i in 2:ncol(snp)) {
    for (j in 1:(i-1)) {
      dij = c(dij,
              sum(snp[,i] != snp[,j])  # nb of differences between chromosome i and chromosome j
      )
    }
  }
  theta.taj=mean(dij)     # = sum(dij)/(n(n-1)/2)
  return(theta.taj)
}


# Theta Tajima, combn version  
# ?combn if needed
theta_tajima_est = function(snp) {
  dij = combn(ncol(snp),m = 2,FUN = function(pair){ sum(snp[,pair[1]] != snp[,pair[2]])})
  theta.taj=mean(dij)   
  return(theta.taj)
}


# Results: 

mu_chrom = 35e06 * 2.5e-08  # full dataset loaded
# Cange to 0.27 if you have loaded only 10000 rows

theta.taj = theta_tajima_est(snp[,sample.info[,"REGION"]=="AFR"])
N.taj = theta.taj / (4*mu_chrom)
cat("AFR",theta.taj, N.taj)

theta.taj = theta_tajima_est(snp[,sample.info[,"REGION"]=="EUR"])
N.taj = theta.taj / (4*mu_chrom)
cat("EUR",theta.taj, N.taj)

theta.wat = theta_watterson_est(snp[,sample.info[,"REGION"]=="AFR"])
N.wat = theta.wat / (4*mu_chrom)
cat("AFR",theta.wat, N.wat)

theta.wat = theta_watterson_est(snp[,sample.info[,"REGION"]=="EUR"])
N.wat = theta.wat / (4*mu_chrom)
cat("EUR",theta.wat, N.wat)



###############################
### Site-Frequency-Spectrum ### 
###############################

# this function returns the number of derived alleles ('1') at each polymorphic site
compute_allele_counts = function(snp) {
  allele.counts = rowSums(snp)
  nsamples=ncol(snp)
  fixed = allele.counts %in% c(0,nsamples)  # sites that are not polymorphic
  return(allele.counts[!fixed])
}



# Loop for plotting SFS for each Region 
for (selected in ureg) {      
  mask = sample.info[,"REGION"] %in% selected  
  n = sum(mask)  # nb of selected individuals
  AC = compute_allele_counts(snp[,mask])
  AC_tab = table(AC)  # absolute number of sites in each class (?table)
  afs = as.numeric(AC_tab)/sum(AC_tab)  # proportion of ...
  
  constant = sum( 1/ (1:(n-1)) )  # the constant term in the expected SFS formula (cf lecture)
  expected = sapply(1:(n-1), function(j) (1/j) /constant )   # expected value for each class from 1 derived to n-1 derived  #equiv: (1/(1:(n-1))) /constant
  plot(expected,xlim=c(0,n),ylim=c(0,0.6),type="b",
       xlab="count",ylab="Proportion",main=selected)
  lines(as.numeric(names(AC_tab)), afs ,type="l", col=colreg[selected])
  legend("top",c("Expected SFS","Observed SFS"),col=c("black",colreg[selected]),lty=1)
}

# In the previsous code, change ureg to upop, and sample.info[,"REGION"] to
# sample.info[,"POP"]  to plot populations' SFS


#############
#### ACP ####
#############

allele.counts = rowSums(snp)
nsamples=ncol(snp)
fixed = allele.counts %in% c(0,nsamples)  # sites that are no polymorphic
snp.p = snp[!fixed,]   # removing fixed sites
pc = prcomp(t(snp.p),center=T,scale=T)
pc = prcomp(t(snp.p),center=T) # if previous line does not work


# plotting first 2 axes of PCA

plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=sample.colreg,sample.info[,"REGION"])

plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=sample.colpop,sample.info[,"POP"])

