
## GENETIC DRIFT ##
genetic_drift_memalloc = function(n,p,gen,chain=rep(NA,gen+1)) {
  #recursive+memory alloc at beginning
  chain[gen+1]=p
  if (gen==0) {
    return(rev(chain))
  } else {
    p1=rbinom(1,size=n,prob=p)/n  #1 is for 1 sample
    genetic_drift_memalloc(n,p1,gen-1,chain)
  }
  
}

genetic_drift = function(n,p,gen,chain=NULL) {
  #recursive version
  chain=c(chain,p)
  if (gen==0) {
    return(chain)
  } else {
    p1=rbinom(1,size=n,prob=p)/n  #1 is for 1 sample
    genetic_drift(n,p1,gen-1,chain)
  }
  
}

genetic_drift_loop = function(n,p,gen) {
  chain=c(p,rep(NA,gen))
  for (g in 1:gen) {
    p1 = rbinom(1,size=n,prob=p)/n  #1 is for 1 sample
    chain[g+1] = p1
  }
  return(chain)
}

n=1000
plot( genetic_drift(n=n,p=0.5,gen=500), type="l" ,ylim=c(0,1),
      main=paste("Size",n), ylab="Allele A frequency",xlab="Time (generations)")
for (rep in 1:20) {
  lines( genetic_drift(n=n,p=0.5,gen=500))
}




## PROBA COAL n=2 ##
prob_coal_at = function(t,N) {
  return((1-1/(2*N))**(t-1)*1/(2*N))
}

prob_coal_after = function(t,N) {
  return((1-1/(2*N))**t)
}
plot(1:1000,prob_coal_after(1:1000,N=100))
plot(1:1000,prob_coal_after(1:1000,N=1000))
abline(h=.1)
abline(h=.5)

seuil = function(prob,N) {
  log(prob)/log(1-1/(2*N))
}
xaxis=seq(1,100000,by=500)
plot(xaxis,prob_coal_after(xaxis,N=10000))
abline(h=.1)
abline(v=seuil(.1,10000),col="red")


## Approximation ##
xx=seq(-1,1,len=100)
t=1
plot(xx,(1-xx)**t)
lines(xx,exp(-xx*t))

##################
##### CG data ####
##################

snp = read.table("../GdP-material/chr22.CG_54genomes_shapeit_phased.txt.haps")#,nrows=100000)  # at first load a smaller nb of rows, eg 1000
sample.info = read.table("../GdP-material/CG_54genomes_indiv.txt",header=T,stringsAsFactors = F)
head(snp)
head(sample.info)
info = snp[,1:5]
snp = snp[,-(1:5)]
nsamples = ncol(snp)

cat("dim:", dim(snp))

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

select.reg=c("EUR","SAS","EAS")
mask = sample.info[,"REGION"]%in%select.reg

odd = seq(1,ncol(snp),by=2)
het = sapply(odd, FUN=function(i) mean(snp[,i]!=snp[,i+1]) )
plot(het,pch=19,type="n")
text(het,sample.info[odd,"POP"],col=sample.colpop[odd])

#plot(het/nrow(snp),pch=19,type="n")
#text(het/nrow(snp),sample.info[odd,"POP"],col=sample.colpop[odd])

rank=order(het)
plot(het[rank],pch=19,type="n")
text(het[rank],sample.info[odd,"POP"][rank],col=sample.colreg[odd][rank])
legend("topleft",names(colreg),pch=19,col=colreg,ncol=5)

#############
### THETA ###
#############

### Estimateur de theta (Tajima)

theta_tajima_est = function(snp) {
  allele.counts = rowSums(snp)
  fixed = allele.counts %in% c(0,ncol(snp))  # sites that are no polymorphic
  dij = apply(combn(ncol(snp),2),FUN=function(pair){ sum(snp[!fixed,pair[1]] != snp[!fixed,pair[2]])} , MARG=2)
  theta.taj=mean(dij)     # sum(dij)/(n(n-1)/2)
return(theta.taj)
}
  # ou une astuce:
  #apply(combn(ncol(snp),2)[,1:10],FUN=function(pair){ sum(rowSums(snp[,pair])%%2==1)} , MARG=2)


mu_chrom = 35e06 * 2.5e-08 
# Change to 0.27 if you have loaded only 10000 rows

theta.taj = theta_tajima_est(snp[,sample.info[,"REGION"]=="AFR"])
N.taj = theta.taj / (4*mu_chrom)
cat("AFR",theta.taj, N.taj)

theta.taj = theta_tajima_est(snp[,sample.info[,"REGION"]=="EUR"])
N.taj = theta.taj / (4*mu_chrom)
cat("EUR",theta.taj, N.taj)

### Estimateur de theta (Watterson)

theta_watterson_est = function(snp) {
  allele.counts = rowSums(snp)
  fixed = allele.counts %in% c(0,ncol(snp))  # sites that are no polymorphic
  S = sum(!fixed) # observed number of segragating sites
  theta.wat = S / sum( 1/1:(ncol(snp)-1)) 
  return(theta.wat)
}

theta.wat = theta_watterson_est(snp[,sample.info[,"REGION"]=="AFR"])
N.wat = theta.wat / (4*mu_chrom)
cat("AFR",theta.wat, N.wat)

theta.wat = theta_watterson_est(snp[,sample.info[,"REGION"]=="EUR"])
N.wat = theta.wat / (4*mu_chrom)
cat("EUR",theta.wat, N.wat)

# difference betw estimators due to chance, but if too big, due to selection (Tajima D, Antoine)

dmat=matrix(NA,ncol(snp),ncol(snp))
allele.counts = rowSums(snp)
fixed = allele.counts %in% c(0,ncol(snp))  # sites that are no polymorphic
for (i in 2:ncol(snp)) {
  for (j in 1:i) {
    dij = sum(snp[!fixed,i] != snp[!fixed,j])
    dmat[i,j] =  dij
    dmat[j,i] =  dij
  }
}

heatmap(dmat)
dmat2=dmat
diag(dmat2)=NA
image(dmat2)

heatmap(dmat2, Rowv = NA, Colv = NA, scale = "column",
        main = "", symm=T)




### ALLELE FREQUENCIES
allele.counts = rowSums(snp)
fixed = allele.counts %in% c(0,nsamples)  # sites that are no polymorphic
hist(allele.counts[!fixed],breaks = 1:nsamples,main="Allele counts",prob=T)


compute_allele_counts = function(snp) {
  allele.counts = rowSums(snp)
  nsamples=ncol(snp)
  fixed = allele.counts %in% c(0,nsamples)  # sites that are not polymorphic
  return(allele.counts[!fixed])
}



### SFS for each Region ###
for (selected in ureg) {
  mask = sample.info[,"REGION"] %in% selected
  n = sum(mask)  # nb selected individuals
  AC = compute_allele_counts(snp[,mask])
  constant = sum( 1/ (1:(n-1)) )
  expected = sapply(1:(n-1), function(j) (1/j) /constant )   # equiv: (1/(1:(n-1))) /constant
  plot(expected,xlim=c(0,n),ylim=c(0,0.6),type="b",
       xlab="count",ylab="Proportion",main=selected)
  afs = table(AC)
  lines(as.numeric(names(afs)), as.numeric(afs)/length(AC),type="l", col=colreg[selected])
  legend("top",c("Expected SFS","Observed SFS"),col=c("black",colreg[selected]),lty=1)
}



for (selected in c("CEU","YRI","ASW","MXL")) {   #upop) {
  mask = sample.info[,"POP"] %in% selected
  n = sum(mask)  # nb selected individuals
  AC = compute_allele_counts(snp[,mask])
  constant = sum( 1/ (1:(n-1)) )
  expected = sapply(1:(n-1), function(j) (1/j) /constant )   # equiv: (1/(1:(n-1))) /constant
  plot(expected,xlim=c(0,n),ylim=c(0,0.4),type="b",
       xlab="count",ylab="Proportion",main=selected)
  afs = table(AC)
  lines(as.numeric(names(afs)), as.numeric(afs)/length(AC),type="l", col=colpop[selected],lwd=2)
  legend("top",c("Expected SFS","Observed SFS"),col=c("black",colpop[selected]),lty=1)
}




#### ACP ####
allele.counts = rowSums(snp)
fixed = allele.counts %in% c(0,nsamples)  # sites that are no polymorphic
snp.p = snp[!fixed,]
pc = prcomp(t(snp.p),center=T,scale=T)
pc = prcomp(t(snp.p),center=T)
plot(pc$x[,1:2])

  
# plotting first 2 axes of PCA
plot(pc$x[,1:2],col=sample.colpop,pch=19)
legend("topleft",names(colpop),pch=19,col=colpop,ncol=5)
plot(pc$x[,1:2],col=sample.colreg,pch=19)
legend("topleft",names(colreg),pch=19,col=colreg,ncol=5)

plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=sample.colreg,sample.info[,"REGION"])

plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=sample.colreg,sample.info[,"ID"])
legend("topleft",names(colreg),pch=19,col=colreg,ncol=5)

plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=sample.colpop,sample.info[,"POP"])

pcaxes=c(1,3)  #1:2
plot(pc$x[,pcaxes],pch=19,type="n")
text(pc$x[,pcaxes],col=sample.colreg,sample.info[,"ID"])

sample.info[sample.info[,1]=="NA19670",]

select.reg=c("EUR","SAS","EAS")
mask = sample.info[,"REGION"]%in%select.reg
pc = prcomp(t(snp.p[,mask]),center=T,scale=F)
plot(pc$x[,1:2],pch=19,type="n")
text(pc$x[,1:2],col=colpop[sample.info[mask,"POP"]],sample.info[mask,"POP"])

