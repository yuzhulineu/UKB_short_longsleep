Rscript /home1/liyz/Project/PHESANT/PHESANT-master/WAS/phenomeScan.r \
--phenofile="/home1/UKB_pheno4.csv" \
--traitofinterestfile="/home1/sleep_data.csv" \
--confounderfile="/home1/cov_data_use.csv" \
--variablelistfile="/home/PHESANT-master/variable-info/outcome-info.tsv" \
--datacodingfile="/home/PHESANT-master/variable-info/data-coding-ordinal-info.txt" \
--traitofinterest="pheno_binary7" \
--resDir="/home1//result_binary_short/" \
--sensitivity \
--userId="userId" \
--genetic=FALSE \
--mincase=400\

%% details could refer to Millard LAC, et al. Software Application Profile: PHESANT: a tool for performing automated phenome scans in UK Biobank. International Journal of Epidemiology (2017).
%% https://github.com/MRCIEU/PHESANT
