%% 1.1 brain plots of patients and sibling subtypes
clc; clear; close all;
cd('E:\Projects\2022-GROUP\Scripts')
load('R3_sustain_subtype_new.mat','sustain_subtype')
load('R5_neuroimaging_regional_zscore_alltime.mat')
load('regionDescriptions.mat');

group = sustain_subtype.group==3;                                 % select:patients (3) or siblings (1)
subtype_volume = sustain_volume_T1_zscore(group,:);                % select:T1, T2, T3
pic_folder = 'R5_neuroimaging_brain_plot/volume_T1_SZ_zscore/';  % save brain plots to ...

subtype = sustain_subtype.newdata_ml_subtype(group);
subtype_volume_mean = grpstats(subtype_volume,subtype,'mean');      % grpstats ignores missing values


%%%%% brain plots save to R5_neuroimaging_brain_plot/volume_T*_SZ/
cmap = flipud(cbrewer('div','RdBu',256));
cmap(cmap>1)=1; cmap(cmap<0)=0;
clims = [-2, 2];       % zscore

for i=1:2 % two subtypes
    tmp_plot = subtype_volume_mean(i,:)';
    tmp_str = ['subtype_',num2str(i)];
    plotBrain(regionDescriptions.aparc_aseg, tmp_plot, cmap,'atlas','aparc_aseg','limits',clims ,...
        'savePath',[pic_folder,tmp_str]);
end



%% 1.2 brain plots of healthy controls
clc; clear; close all;
cd('E:\Projects\2022-GROUP\Scripts')
load('R3_sustain_subtype_new.mat','sustain_subtype')
load('R5_neuroimaging_regional_zscore_alltime.mat')
load('regionDescriptions.mat');

group = sustain_subtype.group==2;                                % select:controls (2)
subtype_volume = sustain_volume_T3_zscore(group,:);                % select:T1, T2, T3
pic_folder = 'R5_neuroimaging_brain_plot/volume_T3_HC_zscore/';  % save brain plots to ...

subtype_volume_mean = nanmean(subtype_volume);      % grpstats ignores missing values

cmap = flipud(cbrewer('div','RdBu',256));
cmap(cmap>1)=1; cmap(cmap<0)=0;

clims = [-2,2];         % zscore
plotBrain(regionDescriptions.aparc_aseg, subtype_volume_mean, cmap,'atlas','aparc_aseg' ,'limits',clims,...
    'savePath',[pic_folder,'HC']);



%% 1.3 Jonckheere-Terpstra Test between control, sibling and patient subtypes 
clc;clear;close all
cd E:\Projects\2022-GROUP\Scripts
load('R3_sustain_subtype_new.mat', 'sustain_demo_T1','sustain_subtype')
load('R4_region_label_plot.mat')
load('R5_neuroimaging_regional_zscore_alltime.mat','sustain_volume_T1_zscore')
load('R21_siblings_info.mat','subtype_sib');

group = sustain_demo_T1.MF_GROUP;
subtype = sustain_subtype.newdata_ml_subtype;

regional_vol = sustain_volume_T1_zscore;    

% roi = [31,65,40,22,77,59]; % ifg  lowest 6 t-stats
roi = [5,12,6,13,1,7];     % hipp lowest 6 t-stats

Nperm = 10000;
pTtj = nan(length(roi),3);
Ttj = nan(length(roi),3);
for i=1:length(roi)
    i
    %%%% con, sib S1, scz S1 %%%% 
    con = regional_vol(group==0,roi(i));
    sib = regional_vol(group==1 & subtype_sib==0, roi(i));
    scz = regional_vol(group==3 & subtype==0, roi(i));
    [pTtj(i,1),Ttj(i,1)] = jonckheereterpstra_perm_function(con, sib, scz, 2, 1, Nperm);

    %%%% con, sib S2, scz S2 %%%%
    con = regional_vol(group==0,roi(i));
    sib = regional_vol(group==1 & subtype_sib==1, roi(i));
    scz = regional_vol(group==3 & subtype==1, roi(i));
    [pTtj(i,2),Ttj(i,2)] = jonckheereterpstra_perm_function(con, sib, scz, 2, 1, Nperm);
    clear con sib* scz

end

pTtj_zscore_T1_roi_hipp = array2table(pTtj,'VariableNames',{'p_sib_scz_s1', 'p_sib_scz_s2'},'RowNames',region_label_plot(roi));
Ttj_zscore_T1_roi_hipp = array2table(Ttj,'VariableNames',{'t_sib_scz_s1', 't_sib_scz_s2'},'RowNames',region_label_plot(roi));
[h,crit_p]=fdr_bh([pTtj_zscore_T1_roi_ifg.p_sib_scz_s1, pTtj_zscore_T1_roi_hipp.p_sib_scz_s2],0.05);


%%%% compare with wholebrain Ttj %%%%
roi = 1:82;
pTtj_hipp_wholebrain = pTtj;
Ttj_hipp_wholebrain = Ttj_empirical; clear Ttj_empirical

ifg_95perc  = prctile(Ttj_ifg_wholebrain, 95)  % 3.2706
hipp_95perc = prctile(Ttj_hipp_wholebrain, 95) % 1.3230


%%%% plot wholebrain Ttj %%%%
close all;
pic_folder = 'R7_neuroimaging_regional_fitlme_global_new_age_CON_SIB/'; 
cmap = [[255,127,0]; [152,78,163]]/255;
numBins = 12;alpha = 0.6;

figure; set(gcf,'Position',[100,100,400,400])
histogram(Ttj_ifg_wholebrain,numBins, 'FaceColor', cmap(1,:), 'EdgeColor', 'w', 'FaceAlpha', alpha);
title('Whole-brain J-T statistic'); 
xlim([-3.4,4.8]);ylim([0 20]);
set(gca,'FontSize',14); box off;
line([ifg_95perc ifg_95perc],[0 6.5],'Color','r','LineStyle','--','LineWidth',2)
saveas(gcf,[pic_folder,'IFG_Wholebrain_JT'],'tiffn')

figure; set(gcf,'Position',[100,100,400,400])
histogram(Ttj_hipp_wholebrain,numBins, 'FaceColor', cmap(2,:), 'EdgeColor', 'w', 'FaceAlpha', alpha);
title('Whole-brain J-T statistic'); 
xlim([-4.5,4.5]);ylim([0 30]);
set(gca,'FontSize',14); box off;
line([hipp_95perc hipp_95perc],[0 7],'Color','r','LineStyle','--','LineWidth',2)
saveas(gcf,[pic_folder,'Hipp_Wholebrain_JT'],'tiffn')
