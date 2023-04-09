#########################
# Packages and software #
#########################
#you need to have 'JAGS - Just Another Gibbs Sampler' installed:
#https://sourceforge.net/projects/mcmc-jags/
require(BEST)  #this package requires JAGS
require(compositions)
require(robCompositions)
require(MASS)
require(metRology)

################
# Data loading #
################
data = read.csv2("metabolites.csv", dec=".") #MCADD data set
data = data[,2:length(data[1,])] 
colndata = colnames(data) = round(as.numeric(substr(colnames(data), 2, 11)),4)

d = length(data[1,]) #number of metabolites d=278
n1 = n2 = 25 #number of samples in the groups
cont = data[1:n1,] #control group
pat = data[n1+1:n2,] #patient group

clr_pat = clr(constSum(pat))
clr_cont = clr(constSum(cont))

#####################
# Testing normality #
#####################
shap_p=c(1:d) 
for (i in 1:d)
  shap_p[i] = shapiro.test(clr_pat[,i])$p.value
length(which(shap_p<0.05))  #normality rejected for 113 metabolites in the group P

shap_c=c(1:d)
for (i in 1:d)
  shap_c[i] = shapiro.test(clr_cont[,i])$p.value
length(which(shap_c<0.05)) #normality rejected for 207 metabolites in the group C
# --> let's choose non-parametric approach

##################
# Wilcoxon tests #
##################
alfa=0.05/d #Bonferroni correction
wil_test=c(1:d)
med=c(1:d)
for (i in 1:d)
{  
  wil_test[i] = wilcox.test(clr_pat[,i], clr_cont[,i])$p.value
  med[i] = median(clr_pat[,i])-median(clr_cont[,i])
}

#####################################################################
# Traditional Volcano plot for Wilcoxon multiple hypotheses testing #
#####################################################################
pdf("Volcano_Wilcoxon.pdf",height=5,width=9)
plot(med,-log10(wil_test),pch=16,cex=0.5,xlab="Difference of medians",ylab="-log10(p-value)")
barevne=which(((wil_test < alfa) & (med > 0.5)) | ((wil_test < alfa) & (med < -0.6)))
text(med[barevne],-log10(wil_test[barevne]),labels=colnames(data[barevne]),pos=4,cex=0.7)
points(med[barevne],-log10(wil_test[barevne]),col="red",pch=17,cex=0.7)
abline(v=0,col="gray")
title("Volcano plot for Wilcoxon tests")
dev.off()


######################
# Bayesian t-testing #
######################
rsh = numeric()
tbval1 = numeric()
tbval2 = numeric()
ebval1 = numeric()
mk   = matrix(ncol=2,nrow=d)

pdf("Bayes_MCMC_histograms.pdf",height=10,width=32)
par(mfrow=c(1,2))
for (i in 1:d)
{
  b_t_test = BESTmcmc(clr_pat[,i], clr_cont[,i]) 
  diffmus = b_t_test$mu1 - b_t_test$mu2
  fitt = fitdistr(diffmus, "t")  #fitting theoretical t-distribution because of potential zero b-values
  k = as.numeric(fitt$estimate)  #parameters of fitted t-distribution
  tbval1[i] = c(1-pt.scaled(0,mean=k[1],df=k[3],sd=k[2])) #controls
  tbval2[i] = c(pt.scaled(0,mean=k[1],df=k[3],sd=k[2]))  #patients
  rsh[i] = mean(diffmus)
  mk[i, ] = as.numeric(quantile(diffmus,c(0.025,0.975)))
  ebval1[i] = sum(diffmus>0)/length(diffmus) 

  plot(b_t_test, title = "Difference of means in clr coef.s")
  plot(b_t_test, which="sd", title = "Difference od std. dev.s 
     in clr coef.s")
  plotPostPred(b_t_test) 
  #summary(b_t_test)
  #print(b_t_test)
}
dev.off()

tbval = pmin(tbval1,tbval2)   #"theoretical" b-values
ebval = pmin(ebval1,1-ebval1)   #"empirical" b-values
#most10t=c(which(tbval1 %in% sort(tbval1)[1:10]),which(tbval2 %in% sort(tbval2)[1:10])) #10 most significant patient and 10 control biomarkers
#most10e=c(which(ebval1 %in% sort(ebval1)[1:10]),which((1-ebval1) %in% sort(1-ebval1)[1:10]))  #most significant patient and control biomarkers

saveRDS(mk, file="HDI_boudaries.RDS")
saveRDS(tbval, file="Theoretical_b-values.RDS")
saveRDS(ebval, file="Empirical_b-values.RDS")
saveRDS(rsh, file="MPDs.RDS")
saveRDS(colnames(data), file="Names_of_metabolites.RDS")
#save.image(file="Data_for_Volcano.RData")


#########################
# Bayesian volcano plot #
#########################
mk = readRDS("HDI_boudaries.RDS")
tbval = readRDS("Theoretical_b-values.RDS")
ebval = readRDS("Empirical_b-values.RDS")
rsh = readRDS("MPDs.RDS")
colndata = readRDS("Names_of_metabolites.RDS")
#load("Data_for_Volcano.RData")

HDI0   = which(mk[,1]<0 & mk[,2]>0) #HDI contains zero
HDI0.5 = which((abs(mk[,1]-0)<0.5 | abs(mk[,2]-0)<0.5) & !(mk[,1]<0 & mk[,2]>0)) #HDI has the boundary close to zero
HDI1   = which((abs(mk[,1]-0)<1 | abs(mk[,2]-0)<1) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<0.5 | abs(mk[,2]-0)<0.5)) #less close... 
HDI1.5 = which((abs(mk[,1]-0)<1.5 | abs(mk[,2]-0)<1.5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<1 | abs(mk[,2]-0)<1)) 
HDI2   = which((abs(mk[,1]-0)<2 | abs(mk[,2]-0)<2) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<1.5 | abs(mk[,2]-0)<1.5)) 
HDI2.5 = which((abs(mk[,1]-0)<2.5 | abs(mk[,2]-0)<2.5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<2 | abs(mk[,2]-0)<2)) 
HDI3   = which((abs(mk[,1]-0)<3 | abs(mk[,2]-0)<3) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<2.5 | abs(mk[,2]-0)<2.5)) 
HDI3.5 = which((abs(mk[,1]-0)<3.5 | abs(mk[,2]-0)<3.5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<3 | abs(mk[,2]-0)<3)) 
HDI4   = which((abs(mk[,1]-0)<4 | abs(mk[,2]-0)<4) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<3.5 | abs(mk[,2]-0)<3.5)) 
HDI4.5 = which((abs(mk[,1]-0)<4.5 | abs(mk[,2]-0)<4.5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<4 | abs(mk[,2]-0)<4)) 
HDI5   = which((abs(mk[,1]-0)<5 | abs(mk[,2]-0)<5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<4.5 | abs(mk[,2]-0)<4.5)) 
HDI5.5 = which((abs(mk[,1]-0)<5.5 | abs(mk[,2]-0)<5.5) & !(mk[,1]<0 & mk[,2]>0) & !(abs(mk[,1]-0)<5 | abs(mk[,2]-0)<5)) 
#length(HDI0)+length(HDI0.5)+length(HDI1)+length(HDI1.5)+length(HDI2)+length(HDI2.5)+length(HDI3)+length(HDI3.5)+length(HDI4)+length(HDI4.5)+length(HDI5)+length(HDI5.5)  #needs to be equal to d

if (min(ebval)==0) {
  #using the "theoretical" b-values
  bval=tbval
} else {
  #using the "empirical" b-values
  bval=ebval
}


pdf("Bayesian_volcano.pdf",height=5,width=9)
xl=ceiling(max(abs(rsh)))
yl=ceiling(max(-log10(bval)))
plot(rsh,-log10(bval),pch=16,cex=0.5,xlim=c(-xl,xl+2),ylim=c(0,yl),xlab="Difference of means",ylab="-log10(b-value)",yaxt="n",xaxt="n",col="white")

points(rsh[HDI0],-log10(bval[HDI0]),col="firebrick",pch=16,cex=0.5)
points(rsh[HDI0.5],-log10(bval[HDI0.5]),col="firebrick1",pch=17,cex=0.5)
points(rsh[HDI1],-log10(bval[HDI1]),col="darkorange",pch=15,cex=0.5)
points(rsh[HDI1.5],-log10(bval[HDI1.5]),col="gold",pch=16,cex=0.5)
points(rsh[HDI2],-log10(bval[HDI2]),col="yellowgreen",pch=17,cex=0.5)
points(rsh[HDI2.5],-log10(bval[HDI2.5]),col="seagreen",pch=15,cex=0.5)
points(rsh[HDI3],-log10(bval[HDI3]),col="mediumturquoise",pch=16,cex=0.5)
points(rsh[HDI3.5],-log10(bval[HDI3.5]),col="deepskyblue1",pch=17,cex=0.5)
points(rsh[HDI4],-log10(bval[HDI4]),col="dodgerblue3",pch=15,cex=0.5)
points(rsh[HDI4.5],-log10(bval[HDI4.5]),col="dodgerblue3",pch=15,cex=0.5)
points(rsh[HDI5],-log10(bval[HDI5]),col="darkslateblue",pch=15,cex=0.5)
points(rsh[HDI5.5],-log10(bval[HDI5.5]),col="darkslateblue",pch=16,cex=0.5)
#text(rsh[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)],-log10(bval)[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)],labels=colnames(data[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)]),pos=2,cex=0.5)
text(rsh[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)],-log10(bval)[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)],labels=colndata[c(HDI5.5,HDI5,HDI4.5,HDI4,HDI3.5,HDI3,HDI2.5,HDI2,HDI1.5)],pos=2,cex=0.5)

abline(v=0,col="gray")
legend("topright", legend = c("HDI contains 0", "|HDI-0|<0.5", "|HDI-0|<1", "|HDI-0|<1.5", "|HDI-0|<2", "|HDI-0|<2.5", "|HDI-0|<3", "|HDI-0|<3.5", "|HDI-0|<4","|HDI-0|<4.5", "|HDI-0|<5", "|HDI-0|<5.5"), 
       pch = c(16, 17, 15, 16, 17, 15, 16, 17, 15, 15, 15, 16), col = c("firebrick", "firebrick1", "darkorange", "gold","yellowgreen","seagreen","mediumturquoise","deepskyblue1","dodgerblue3","white","white","darkslateblue"),cex=0.7)    
title("Bayesian Volcano plot with HDI zones")
axis(side = 1, at = c((-xl+1):(xl-1)), labels = c((-xl+1):(xl-1)))
axis(side = 2, at = c(0:(yl-1)), labels = c(0:(yl-1)))
dev.off()
  







