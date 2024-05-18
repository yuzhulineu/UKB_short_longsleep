% the code was referred to http://www.bigc.online/BrainMR/Browse/1.2-Code.html
% Guo, J., et al., Nature Neuroscience, 2022. 25(11): p. 1519-1527. https://doi.org/10.1038/s41593-022-01174-7

library(RadialMR)

data1<-read.table(input1,header=T)
data2<-read.table(input2,header=T)

data3<-format_radial(data1$B,data2$B,data1$SE,data2$SE,data1$SNP)

res1<-ivw_radial(data3,0.05,1,0.0001)
res2<-egger_radial(data3,0.05,1)
d1<-res1$data
d2<-res2$data
d1<-as.data.frame(d1)
d2<-as.data.frame(d2)
write.table(d1,output1,quote=F,row.names=F,col.names=T,sep="\t")
write.table(d2,output2,quote=F,row.names=F,col.names=T,sep="\t")
a1<-cbind(res1$outliers[1])
a2<-cbind(res2$outliers[1])

library(TwoSampleMR)

exp_dat <- read_exposure_data( filename = input1,sep = "\t",snp_col = "SNP",beta_col ="B",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "MAF",pval_col = "P",samplesize_col = "N")

outcome_dat <- read_outcome_data(snps = exp_dat$SNP,filename = input2,sep = "\t",snp_col = "SNP",beta_col = "B",se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "MAF",pval_col = "P",samplesize_col = "N")

dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat = outcome_dat,action=1)

tsmr1<-mr(dat, method_list=c("mr_weighted_median"))
tsmr2<-mr(dat, method_list=c("mr_ivw"))
tsmr3<-mr(dat, method_list=c("mr_egger_regression"))


