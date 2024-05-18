% this is a matlab code utilized for the calculation between sleep and blood biomarkers
clc,clear

load sleepnewall
subjid=table2array(sleepnewall(:,1));
sleep0=table2array(sleepnewall(:,2));
sleep0_data=sleep0(sleep0>0);
sleeprepid0=subjid(sleep0>0);


index1=find(sleep0_data<=7);
index2=find(sleep0_data>=7);

sleeprepid7= sleepid0(index1);
sleeprepid10= sleepid0(index2);

sleepreport7= sleep0_data(index1);
sleepreport10=sleep0_data(index2);

load cov_ukb_int0
cov_ukb0=table2array(cov_ukb_int0(:,2:end));
cov_id0=table2array(cov_ukb_int0(:,1));


load blood_data
bioid=table2array(blood_data(:,1));


C_value_bio=nan(size(blood_data,2)-1,1);
P_value_bio=nan(size(blood_data,2)-1,1);

for ii=2:size(blood_data,2)

    biodata=table2array(blood_data(:,ii));
    phenodataori= biodata;
    phenodata = phenodataori(~isnan(phenodataori));
    pheno_id = bioid(~isnan(phenodataori));

    final_id=intersect_multi({sleeprepid7;cov_id0;pheno_id});

    [~,index]=intersect(cov_id0,final_id);

    cov_1=cov_ukb0(index,:);
    cov_final1=nan(size(cov_1));
    for kk=1:size(cov_1,2)
        cov1=cov_1(:,kk);
        cov1(isnan(cov1))=nanmean(cov1);
        cov_final1(:,kk)=cov1;
    end

    [~,index]=intersect(sleeprepid7,final_id);
    sleepmeasure_final=sleepreport7(index,:);

    [~,index]=intersect(pheno_id,final_id);
    pheno_final=phenodata(index,:);

    Covariate = cov_final1;


    [C_value_bio(ii-1),P_value_bio(ii-1)]=partialcorr(sleepmeasure_final,pheno_final,Covariate,'rows','pairwise');

end
bio_name=blood_data.Properties.VariableNames(2:end)';
result_sleepblood7=[cell2table(bio_name),array2table(C_value_bio), array2table(P_value_bio)];

%% cal the blood biomarkers with sleep duration >=7
C_value_bio10=nan(size(blood_data,2)-1,1);
P_value_bio10=nan(size(blood_data,2)-1,1);
for ii=2:size(blood_data,2)

    biodata=table2array(blood_data(:,ii));
    phenodataori= biodata;
    phenodata = phenodataori(~isnan(phenodataori));
    pheno_id = bioid(~isnan(phenodataori));

    final_id=intersect_multi({sleeprepid10;cov_id0;pheno_id});

    [~,index]=intersect(cov_id0,final_id);
    cov_1=cov_ukb0(index,:);
    cov_final1=nan(size(cov_1));
    for kk=1:size(cov_1,2)
        cov1=cov_1(:,kk);
        cov1(isnan(cov1))=nanmean(cov1);
        cov_final1(:,kk)=cov1;
    end

    [~,index]=intersect(sleeprepid10,final_id);
    sleepmeasure_final2=sleepreport10(index,:);

    [~,index]=intersect(pheno_id,final_id);
    pheno_final2=phenodata(index,:);

    Covariate = cov_final1;

    [C_value_bio10(ii-1),P_value_bio10(ii-1)]=partialcorr(sleepmeasure_final2,pheno_final2,Covariate,'rows','pairwise');

end
result_sleepbio10=[cell2table(bio_name),array2table(C_value_bio10), array2table(P_value_bio10)];

result_sleepbio={result_sleepblood7,result_sleepbio10};
%%  age stratify analysis sleep<=7

C_value_bio=nan(size(bloodchem_data,2)-1,1);
P_value_bio=nan(size(bloodchem_data,2)-1,1);

C_value_bio_old=nan(size(bloodchem_data,2)-1,1);
P_value_bio_old=nan(size(bloodchem_data,2)-1,1);

for ii=2:size(bloodchem_data,2)

    biodata=table2array(bloodchem_data(:,ii));
    phenodataori= biodata;
    phenodata = phenodataori(~isnan(phenodataori));
    pheno_id = bioid(~isnan(phenodataori));

    final_id=intersect_multi({sleeprepid7;cov_id0;pheno_id});

    [~,index]=intersect(cov_id0,final_id);
    %cov_1=cov_ukb0(index,1:2);%age sex townsend ethic bmi edu
    cov_1=cov_ukb0(index,:);
    cov_final1=nan(size(cov_1));
    for kk=1:size(cov_1,2)
        cov1=cov_1(:,kk);
        cov1(isnan(cov1))=nanmean(cov1);
        cov_final1(:,kk)=cov1;
    end

    [~,index]=intersect(sleeprepid7,final_id);
    sleepmeasure_final=sleepreport7(index,:);

    [~,index]=intersect(pheno_id,final_id);
    pheno_final=phenodata(index,:);


    index_young=find(cov_final1(:,1)<58);

    cov_young=cov_final1(index_young,:);
    sleepmeasure_young=sleepmeasure_final(index_young,:);
    pheno_young=pheno_final(index_young,:);

    [C_value_bio(ii-1),P_value_bio(ii-1)]=partialcorr(sleepmeasure_young,pheno_young,cov_young,'rows','pairwise');


    index_old=find(cov_final1(:,1)>=58);

    cov_old=cov_final1(index_old,:);
    sleepmeasure_old=sleepmeasure_final(index_old,:);
    pheno_old=pheno_final(index_old,:);

    [C_value_bio_old(ii-1),P_value_bio_old(ii-1)]=partialcorr(sleepmeasure_old,pheno_old,cov_old,'rows','pairwise');
end
bio_name=blood_data.Properties.VariableNames(2:end)';
result_sleepbio7_young=[cell2table(bio_name),array2table(C_value_bio), array2table(P_value_bio)];
result_sleepbio7_old=[cell2table(bio_name),array2table(C_value_bio_old), array2table(P_value_bio_old)];
%% age stratify analysis sleep>=7
C_value_bio_young10=nan(size(blood_data,2)-1,1);
P_value_bio_young10=nan(size(blood_data,2)-1,1);

C_value_bio_old10=nan(size(blood_data,2)-1,1);
P_value_bio_old10=nan(size(blood_data,2)-1,1);
for ii=2:size(blood_data,2)

    biodata=table2array(blood_data(:,ii));
    phenodataori= biodata;
    phenodata = phenodataori(~isnan(phenodataori));
    pheno_id = bioid(~isnan(phenodataori));

    final_id=intersect_multi({sleeprepid10;cov_id0;pheno_id});

    [~,index]=intersect(cov_id0,final_id);
    cov_1=cov_ukb0(index,:);%age sex townsend ethic bmi edu
    cov_final1=nan(size(cov_1));
    for kk=1:size(cov_1,2)
        cov1=cov_1(:,kk);
        cov1(isnan(cov1))=nanmean(cov1);
        cov_final1(:,kk)=cov1;
    end

    [~,index]=intersect(sleeprepid10,final_id);
    sleepmeasure_final2=sleepreport10(index,:);

    [~,index]=intersect(pheno_id,final_id);
    pheno_final2=phenodata(index,:);

    Covariate = cov_final1;

    index_young10=find(cov_final1(:,1)<58);

    cov_young=cov_final1(index_young10,:);
    sleepmeasure_young=sleepmeasure_final2(index_young10,:);
    pheno_young=pheno_final2(index_young10,:);

    [C_value_bio_young10(ii-1),P_value_bio_young10(ii-1)]=partialcorr(sleepmeasure_young,pheno_young,cov_young,'rows','pairwise');


    index_old10=find(cov_final1(:,1)>=58);

    cov_old=cov_final1(index_old10,:);
    sleepmeasure_old=sleepmeasure_final2(index_old10,:);
    pheno_old=pheno_final2(index_old10,:);

    [C_value_bio_old10(ii-1),P_value_bio_old10(ii-1)]=partialcorr(sleepmeasure_old,pheno_old,cov_old,'rows','pairwise');
end
result_sleepbio_young10=[cell2table(bio_name),array2table(C_value_bio_young10), array2table(P_value_bio_young10)];
result_sleepbio_old10=[cell2table(bio_name),array2table(C_value_bio_old10), array2table(P_value_bio_old10)];

