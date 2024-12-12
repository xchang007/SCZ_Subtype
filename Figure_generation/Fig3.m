%% 1.Sustain subtype comparison clozapine 
clc; clear; close all;
load('R3_sustain_subtype_new.mat')
load('R6_sustain_medication.mat')

group = sustain_subtype.group==3;
subtype = sustain_subtype.newdata_ml_subtype(group);

%%%%% check if sustain subtypes differ in using clozapine %%%%%
aptype = 4; % clozapine (Leponex)
tf_CPZ_T1 = double(sustain_demo_T1.MF_APNAME_1 == aptype | sustain_demo_T1.MF_APNAME_2 == aptype |sustain_demo_T1.MF_APNAME_3 == aptype | sustain_demo_T1.MF_APNAME_4 == aptype);
tf_CPZ_T1(isnan(sustain_demo_T1.MF_APNAME_1)&isnan(sustain_demo_T1.MF_APNAME_2)&isnan(sustain_demo_T1.MF_APNAME_3)&isnan(sustain_demo_T1.MF_APNAME_4))=nan;
tf_CPZ_T1(sustain_demo_T1.MF_AP_STAT==0)=nan;
dose_CPZ_T1 = (tf_CPZ_T1==1).*sustain_demo_T1.MF_APCPZE_1;  % dosage of clozapine
dose_CPZ_T1(tf_CPZ_T1==1 & sustain_demo_T1.MF_APNAME_1~=4) = sustain_demo_T1.MF_APCPZE_2(tf_CPZ_T1==1 & sustain_demo_T1.MF_APNAME_1~=4);
dose_CPZ_T1(dose_CPZ_T1==0)=nan;

tf_CPZ_T2 = double(sustain_demo_T2.MF_APNAME_1 == aptype | sustain_demo_T2.MF_APNAME_2 == aptype |sustain_demo_T2.MF_APNAME_3 == aptype | sustain_demo_T2.MF_APNAME_4 == aptype);
tf_CPZ_T2(isnan(sustain_demo_T2.MF_APNAME_1)&isnan(sustain_demo_T2.MF_APNAME_2)&isnan(sustain_demo_T2.MF_APNAME_3)&isnan(sustain_demo_T2.MF_APNAME_4))=nan;
tf_CPZ_T2(sustain_demo_T1.MF_AP_STAT==0)=nan;
dose_CPZ_T2 = (tf_CPZ_T2==1).*sustain_demo_T2.MF_APCPZE_1;  % dosage of clozapine
dose_CPZ_T2(tf_CPZ_T2==1 & sustain_demo_T2.MF_APNAME_1~=4) = sustain_demo_T2.MF_APCPZE_2(tf_CPZ_T2==1 & sustain_demo_T2.MF_APNAME_1~=4);
dose_CPZ_T2(dose_CPZ_T2==0)=nan;

tf_CPZ_T3 = double(sustain_demo_T3.MF_APNAME_1 == aptype | sustain_demo_T3.MF_APNAME_2 == aptype |sustain_demo_T3.MF_APNAME_3 == aptype | sustain_demo_T3.MF_APNAME_4 == aptype);
tf_CPZ_T3(isnan(sustain_demo_T3.MF_APNAME_1)&isnan(sustain_demo_T3.MF_APNAME_2)&isnan(sustain_demo_T3.MF_APNAME_3)&isnan(sustain_demo_T3.MF_APNAME_4))=nan;
tf_CPZ_T3(sustain_demo_T1.MF_AP_STAT==0)=nan;
dose_CPZ_T3 = (tf_CPZ_T3==1).*sustain_demo_T3.MF_APCPZE_1;  % dosage of clozapine
dose_CPZ_T3(tf_CPZ_T3==1 & sustain_demo_T3.MF_APNAME_1~=4) = sustain_demo_T3.MF_APCPZE_2(tf_CPZ_T3==1 & sustain_demo_T3.MF_APNAME_1~=4);
dose_CPZ_T3(dose_CPZ_T3==0)=nan;

tf_CPZ_T1 = tf_CPZ_T1(group);     tf_CPZ_T2 = tf_CPZ_T2(group);     tf_CPZ_T3 = tf_CPZ_T3(group);
dose_CPZ_T1 = dose_CPZ_T1(group); dose_CPZ_T2 = dose_CPZ_T2(group); dose_CPZ_T3 = dose_CPZ_T3(group);
[nnz(tf_CPZ_T1==1),nnz(tf_CPZ_T2==1),nnz(tf_CPZ_T3==1)] %  21    27    23
% save('R6_sustain_medication.mat','tf_CPZ*','dose_CPZ*','-append')



%%%%% chisquare test %%%%%
tf_CPZ = zeros(length(subtype),1);
idx = (tf_CPZ_T1==1) | (tf_CPZ_T2==1) | (tf_CPZ_T3==1);       % any time use cloz
tf_CPZ(idx) = 1;
idx = isnan(tf_CPZ_T1) & isnan(tf_CPZ_T2) & isnan(tf_CPZ_T3); % all missing
tf_CPZ(idx) = nan;
% save('R6_sustain_medication.mat','tf_CPZ','-append')

[nnz(tf_CPZ==0) nnz(tf_CPZ==1) nnz(isnan(tf_CPZ))] %  101    42    31
[tbl,chi2,p,labels] = crosstab(subtype,tf_CPZ) %chi2 = 5.6213, p = 0.0177
%  CLO OtherAP Missing
% S1 34  61   19
% S2 8   40   12
OR = (34*40)/(8*61) % Odds Ratio 2.7869


[conttbl,chi2,p,labels] = crosstab(subtype,tf_CPZ_T1) % p = 0.0703, chi=3.2753
%     Other   CLO
% S1    69    18
% S2    36     3
OR = (18/69)/(3/36) % Odds Ratio 3.1304

[conttbl,chi2,p,labels] = crosstab(subtype,tf_CPZ_T2) % p = 0.0531, chi=3.7421
%     Other   CLO
% S1    29    20
% S2    27     7
OR = (20/29)/(7/27) % Odds Ratio 2.6601

[conttbl,chi2,p,labels] = crosstab(subtype,tf_CPZ_T3) % p = 0.0403, chi=4.2055
%     Other   CLO
% S1    22    18
% S2    20     5
OR = (18/22)/(5/20) % Odds Ratio 3.2727



%%%%% fisher exact test %%%%%
% among 18 sub T1 otherAP, T2/T3 cloz, 3 hipp, 15 ifg 
% among 21 sub T1 cloz, 3 hipp, 18 ifg

x = array2table ([15,61;3,40],'VariableNames',{'Cloz','OtherAP'},'RowNames',{'S1','S2'});
[h,p,stats] = fishertest(x,'Tail','right') % p=0.0505, OddsRatio: 3.2787    CI: [0.8916 12.0564]

x = array2table ([18,61;3,40],'VariableNames',{'Cloz','OtherAP'},'RowNames',{'S1','S2'});
[h,p,stats] = fishertest(x,'Tail','right') % p=0.0211, OddsRatio: 3.9344    CI: [1.0877 14.2310]

