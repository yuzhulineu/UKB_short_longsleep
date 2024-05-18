# cal the association between sleep and brain structure
clc,clear

load sleepnewall
 subjid=table2array( sleepnewall(:,1));
 sleep2=table2array(sleepnewall(:,4));
 sleepreport=sleep2(sleep2>0);
 sleeprepid=subjid(sleep2>0);
 
 index1=find(sleepreport<=7);
 index2=find(sleepreport>=7);
 sleeprepid7=sleeprepid(index1);
 sleepreport7=sleepreport(index1);
 sleeprepid10=sleeprepid(index2);
 sleepreport10=sleepreport(index2);
 

 load cov_ukb_int2
 cov_ukb2=str2double(cellstr(table2array(cov_ukb_int2(:,2:end))));
 cov_id2=str2double(cellstr(table2array(cov_ukb_int2(:,1))));


load UKB_aparcnew
newimageid=table2array(UKB_aparc_Area(:,1));


UKB_aseg_Volume=readtable('UKB_2_0_aseg_Volume.csv');
imageid_aseg=table2array(UKB_aseg_Volume(:,1));

load site_cov_id.mat


final_id7 = intersect_multi({cov_id2;sleeprepid7;newimageid;imageid_aseg;site_id});


[~,index]=intersect(cov_id2,final_id7);
cov_final= cov_ukb2(index,:);
cov_final1=nan(size(cov_final));
for kk=1:size(cov_final,2)
    cov1=cov_final(:,kk);
    cov1(isnan(cov1))=nanmean(cov1);
    cov_final1(:,kk)=cov1;
end

[~, index] = intersect(site_id, final_id7);
 site_cov_final = site_cov(index,:);


[~,index]=intersect(imageid_aseg,final_id7);
UKB_aseg_final=double(cell2mat(table2cell(UKB_aseg_Volume(index,2:46))));
totalvolume_final=double(cell2mat(table2cell(UKB_aseg_Volume(index,67))));%etiv
asegname=UKB_aseg_Volume.Properties.VariableNames(2:46);

cov_final1 = [cov_final1, site_cov_final, totalvolume_final];
%adjust for imaging sites, etiv and other cov

[~,index]=intersect(sleeprepid7,final_id7);
sleep7_final=sleepreport7(index,:);

[~,index]=intersect(newimageid,final_id7);
UKB_area_final=UKB_aparc_Area(index,2:end);
UKB_thickness_final=UKB_aparc_Thickness(index,2:end);
UKB_volume_final=UKB_aparc_Volume(index,2:end); 
UKB_aparc=[UKB_area_final UKB_thickness_final UKB_volume_final];
UKB_aparc_final=double(cell2mat(table2cell(UKB_aparc)));
aparcname=UKB_aparc.Properties.VariableNames;


C_value_aparctotal = nan(size(UKB_aparc_final,1),1);
P_value_aparctotal= nan(size(UKB_aparc_final,1),1);

for i=1:size(UKB_aparc,2)
    
    Covariate = cov_final1;
   
    [C_value_aparctotal(i,:), P_value_aparctotal(i,:)] = partialcorr(UKB_aparc_final(:,i),sleep7_final,Covariate,'rows','complete');
end

result_sleep7aparccov=[aparcname' num2cell(C_value_aparctotal) num2cell(P_value_aparctotal)];

C_value_asegsleep7 = nan(size(UKB_aseg_final,1),1);
P_value_asegsleep7= nan(size(UKB_aseg_final,1),1);
for i=1:size(UKB_aseg_final,2)

Covariate = cov_final1;
[C_value_asegsleep7(i,1), P_value_asegsleep7(i,1)] = partialcorr(UKB_aseg_final(:,i),sleep7_final,Covariate,'rows','complete');
end
result_aseg_icv7=[asegname' num2cell(C_value_asegsleep7) num2cell(P_value_asegsleep7)];%not regress tIV

save result_aseg_icv7 result_aseg_icv7
%% age stratify

C_value_aparctotal_young = nan(size(UKB_aparc_final,1),1);
P_value_aparctotal_young= nan(size(UKB_aparc_final,1),1);

C_value_aparctotal_old = nan(size(UKB_aparc_final,1),1);
P_value_aparctotal_old= nan(size(UKB_aparc_final,1),1);

for i=1:size(UKB_aparc,2)
    
    index_young=find(cov_final1(:,1)<64);
    
    cov_young=cov_final1(index_young,:);
    UKB_aparc_young=UKB_aparc_final(index_young,:);
    sleep7_young=sleep7_final(index_young);
    
   
    [C_value_aparctotal_young(i,:), P_value_aparctotal_young(i,:)] = partialcorr(UKB_aparc_young(:,i),sleep7_young,cov_young,'rows','complete');
    
    
    index_old=find(cov_final1(:,1)>=64);
    
    cov_old=cov_final1(index_old,:);
    UKB_aparc_old=UKB_aparc_final(index_old,:);
    sleep7_old=sleep7_final(index_old);
    
     [C_value_aparctotal_old(i,:), P_value_aparctotal_old(i,:)] = partialcorr(UKB_aparc_old(:,i),sleep7_old,cov_old,'rows','complete');
    
end

result_sleep7aparc_young=[aparcname' num2cell(C_value_aparctotal_young) num2cell(P_value_aparctotal_young)];
result_sleep7aparc_old=[aparcname' num2cell(C_value_aparctotal_old) num2cell(P_value_aparctotal_old)];


C_value_asegsleep7_young = nan(size(UKB_aseg_final,1),1);
P_value_asegsleep7_young= nan(size(UKB_aseg_final,1),1);

C_value_asegsleep7_old = nan(size(UKB_aseg_final,1),1);
P_value_asegsleep7_old= nan(size(UKB_aseg_final,1),1);

for i=1:size(UKB_aseg_final,2)

    
    index_young=find(cov_final1(:,1)<64);
    
    cov_young=cov_final1(index_young,:);
    UKB_aseg_young=UKB_aseg_final(index_young,:);
    sleep7_young=sleep7_final(index_young);
    

[C_value_asegsleep7_young(i,1), P_value_asegsleep7_young(i,1)] = partialcorr(UKB_aseg_young(:,i),sleep7_young,cov_young,'rows','complete');

    index_old=find(cov_final1(:,1)>=64);
    
    cov_old=cov_final1(index_old,:);
    UKB_aseg_old=UKB_aseg_final(index_old,:);
    sleep7_old=sleep7_final(index_old);
    
    [C_value_asegsleep7_old(i,1), P_value_asegsleep7_old(i,1)] = partialcorr(UKB_aseg_old(:,i),sleep7_old,cov_old,'rows','complete');
end
result_aseg_icv7_young=[asegname' num2cell(C_value_asegsleep7_young) num2cell(P_value_asegsleep7_young)];%not regress tIV
result_aseg_icv7_old=[asegname' num2cell(C_value_asegsleep7_old) num2cell(P_value_asegsleep7_old)];%not regress tIV


save result_aseg_icv7_young result_aseg_icv7_young
save result_aseg_icv7_old result_aseg_icv7_old
%%
final_id10 = intersect_multi({cov_id2;sleeprepid10;newimageid;imageid_aseg;site_id});

[~,index]=intersect(cov_id2,final_id10);
cov_final= cov_ukb2(index,:);
cov_final1=nan(size(cov_final));
for kk=1:size(cov_final,2)
    cov1=cov_final(:,kk);
    cov1(isnan(cov1))=nanmean(cov1);
    cov_final1(:,kk)=cov1;
end

[~, index] = intersect(site_id, final_id10);
 site_cov_final = site_cov(index,:);


[~,index]=intersect(imageid_aseg,final_id10);
UKB_aseg_final=double(cell2mat(table2cell(UKB_aseg_Volume(index,2:46))));
totalvolume_final=double(cell2mat(table2cell(UKB_aseg_Volume(index,67))));
asegname=UKB_aseg_Volume.Properties.VariableNames(2:46);

cov_final1 = [cov_final1, site_cov_final, totalvolume_final];

[~,index]=intersect(sleeprepid10,final_id10);
sleep10_final=sleepreport10(index,:);

[~,index]=intersect(newimageid,final_id10);
UKB_area_final=UKB_aparc_Area(index,2:end);
UKB_thickness_final=UKB_aparc_Thickness(index,2:end);
UKB_volume_final=UKB_aparc_Volume(index,2:end); 
UKB_aparc=[UKB_area_final UKB_thickness_final UKB_volume_final];
UKB_aparc_final=double(cell2mat(table2cell(UKB_aparc)));
aparcname=UKB_aparc.Properties.VariableNames;


C_value_aparcsleep10_young = nan(size(UKB_aparc,2),1);
P_value_aparcsleep10_young= nan(size(UKB_aparc,2),1);

C_value_aparcsleep10_old= nan(size(UKB_aparc,2),1);
P_value_aparcsleep10_old= nan(size(UKB_aparc,2),1);

for i=1:size(UKB_aparc,2)
    
    index_young1=find(cov_final1(:,1)<64);
    cov_young=cov_final1(index_young1,:);
    
    UKB_aparc_young=UKB_aparc_final(index_young1,:);
    sleep10_young=sleep10_final(index_young1);
    
   
    [C_value_aparcsleep10_young(i,:), P_value_aparcsleep10_young(i,:)] = partialcorr(UKB_aparc_young(:,i),sleep10_young,cov_young,'rows','complete');
    
    
    index_old1=find(cov_final1(:,1)>=64);
    
    cov_old=cov_final1(index_old1,:);
    UKB_aparc_old=UKB_aparc_final(index_old1,:);
    sleep10_old=sleep10_final(index_old1);
    
     [C_value_aparcsleep10_old(i,:), P_value_aparcsleep10_old(i,:)] = partialcorr(UKB_aparc_old(:,i),sleep10_old,cov_old,'rows','complete');
    
end
result_sleep10aparc_young=[aparcname' num2cell(C_value_aparcsleep10_young) num2cell(P_value_aparcsleep10_young)];
result_sleep10aparc_old=[aparcname' num2cell(C_value_aparcsleep10_old) num2cell(P_value_aparcsleep10_old)];

result_sleep10_thickness_young=result_sleep10aparc_young(69:136,:);
result_sleep10_thickness_old=result_sleep10aparc_old(69:136,:);

result_sleep10_volume_young=result_sleep10aparc_young(137:end,:);
result_sleep10_volume_old=result_sleep10aparc_old(137:end,:);

result_sleep10_volume=result_sleep10aparccov(137:204,:);
result_sleep10_area=result_sleep10aparccov(1:68,:);
result_sleep10_thickness=result_sleep10aparccov(69:136,:);


C_value_asegsleep10_young = nan(size(UKB_aseg_final,2),1);
P_value_asegsleep10_young= nan(size(UKB_aseg_final,2),1);

C_value_asegsleep10_old = nan(size(UKB_aseg_final,2),1);
P_value_asegsleep10_old= nan(size(UKB_aseg_final,2),1);

for i=1:size(UKB_aseg_final,2)


    index_young=find(cov_final1(:,1)<64);

    cov_young=cov_final1(index_young,:);
    UKB_aseg_young=UKB_aseg_final(index_young,:);
    sleep10_young=sleep10_final(index_young);


    [C_value_asegsleep10_young(i,1), P_value_asegsleep10_young(i,1)] = partialcorr(UKB_aseg_young(:,i),sleep10_young,cov_young,'rows','complete');

    index_old=find(cov_final1(:,1)>=64);

    cov_old=cov_final1(index_old,:);
    UKB_aseg_old=UKB_aseg_final(index_old,:);
    sleep10_old=sleep10_final(index_old);

    [C_value_asegsleep10_old(i,1), P_value_asegsleep10_old(i,1)] = partialcorr(UKB_aseg_old(:,i),sleep10_old,cov_old,'rows','complete');
end
result_aseg_icv10_young=[asegname' num2cell(C_value_asegsleep10_young) num2cell(P_value_asegsleep10_young)];
result_aseg_icv10_old=[asegname' num2cell(C_value_asegsleep10_old) num2cell(P_value_asegsleep10_old)];
