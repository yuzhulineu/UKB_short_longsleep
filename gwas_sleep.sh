/home1/plink2 
--bfile /home/UKB/gene/v3/QCed/ukb_imp_chr${chr}_v3 
--freq 
--linear firth-fallback 
hide-covar 
--covar-variance-standardize 
--pheno /home1/sleep_trait.txt 
--input-missing-phenotype -9 
--covar /home1/covsleep_data.txt  
--out /home1/gwassleep_resultnew/result.chr$chr

%% cal the gwas of sleep duration utlizing plink2
