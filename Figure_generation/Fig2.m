%% 1.1 fitlme in controls and patients subtypes 
% S_func_fitlme_grpwave_global: 'regional~MF_GROUP*Period+MF_AGE+MF_AGE_2+MF_GENDER+MF_ETHN2+ICV+(1|pseudoid)'

clc; clear; close all;
cd('E:\Projects\2022-GROUP\Scripts')
load('R3_sustain_subtype_new.mat')
load('R4_neuroimaging_global.mat', 'ICV*_excl')
load('R5_neuroimaging_regional.mat');
addpath('E:\Projects\2022-GROUP\Scripts\Utility')

tbl = [sustain_demo_T1(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'});...
    sustain_demo_T2(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'});...
    sustain_demo_T3(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'})];

MF_AGE = repmat(sustain_demo_T1.MF_AGE,3,1);
MF_AGE_2 = MF_AGE.^2;
tbl = [tbl, array2table([MF_AGE,MF_AGE_2],'VariableNames',{'MF_AGE','MF_AGE_2'})] ; % baseline age

regional_vol = [sustain_volume_T1;sustain_volume_T2;sustain_volume_T3];    
global_vol = [ICV_T1_excl; ICV_T1_excl; ICV_T1_excl];                               % baseline ICV  


%%%%%% brain volume fitlme, subtype: CONTROL+SUBTYPES(Siblings & Patients) %%%%%%

subtype = sustain_subtype.newdata_ml_subtype;
group = sustain_demo_T1.MF_GROUP;
% [nnz(group==0), nnz(group==1), nnz(group==3)] % 136   202   174
subtype_ifg = (group==0);                       % controls
subtype_ifg(subtype==0 & group==3) = 1;         % patients, IFG
subtype_ifg = repmat(subtype_ifg,3,1);

subtype_hipp = (group==0);                      % controls
subtype_hipp(subtype==1 & group==3) = 1;        % patients, hipp
subtype_hipp = repmat(subtype_hipp,3,1);

[p_fitlme_grpwave_vol_ifg,t_fitlme_grpwave_vol_ifg, coef_fitlme_grpwave_vol_ifg, regional_vol_reg_ifg]=...
    S_func_fitlme_grpwave_global(tbl(subtype_ifg==1,:), regional_vol(subtype_ifg==1,:), global_vol(subtype_ifg==1));

[p_fitlme_grpwave_vol_hipp,t_fitlme_grpwave_vol_hipp, coef_fitlme_grpwave_vol_hipp, regional_vol_reg_hipp]=...
    S_func_fitlme_grpwave_global(tbl(subtype_hipp==1,:), regional_vol(subtype_hipp==1,:), global_vol(subtype_hipp==1));

% save('R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ.mat','regional_*','tbl','*_fitlme_*','region_label')


%% 1.2 fitlme in controls and patients subtypes, excl tf_CPZ 

clc; clear; close all;
cd('E:\Projects\2022-GROUP\Scripts')
load('R3_sustain_subtype_new.mat','sustain_demo_T*','sustain_subtype')
load('R4_neuroimaging_global.mat', 'ICV*_excl')
load('R5_neuroimaging_regional.mat');
load('R6_sustain_medication.mat')

addpath('E:\Projects\2022-GROUP\Scripts\Utility')

tbl = [sustain_demo_T1(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'});...
    sustain_demo_T2(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'});...
    sustain_demo_T3(:,{'pseudoid','MF_GROUP','Period','MF_GENDER','MF_ETHN2'})];
MF_AGE = repmat(sustain_demo_T1.MF_AGE,3,1);
MF_AGE_2 = MF_AGE.^2;
tbl = [tbl, array2table([MF_AGE,MF_AGE_2],'VariableNames',{'MF_AGE','MF_AGE_2'})] ; % baseline age

regional_vol = [sustain_volume_T1;sustain_volume_T2;sustain_volume_T3];    
global_vol = [ICV_T1_excl; ICV_T1_excl; ICV_T1_excl];                               % baseline ICV  

Nsub = size(sustain_subtype,1);

%%%%%% brain volume fitlme, subtype: CONTROL+SUBTYPES(Siblings & Patients) %%%%%%
subtype = sustain_subtype.newdata_ml_subtype;
group = sustain_demo_T1.MF_GROUP;
% [nnz(group==0), nnz(group==1), nnz(group==3)] % 136   202   174

scz_tf_cpz = nan(Nsub,1);
scz_tf_cpz(group==3)=tf_CPZ;
[nnz(scz_tf_cpz==0), nnz(scz_tf_cpz==1)]        % 101    42

subtype_ifg = (group==0);                       % controls
subtype_ifg(subtype==0 & group==3) = 1;         % patients, IFG
subtype_ifg(scz_tf_cpz==1)=0; 
[nnz(subtype_ifg)-nnz(group==0)]                % excl cloz patients, 80=114-34
subtype_ifg = repmat(subtype_ifg,3,1);

subtype_hipp = (group==0);                      % controls
subtype_hipp(subtype==1 & group==3) = 1;        % patients, hipp
subtype_hipp(scz_tf_cpz==1)=0; 
[nnz(subtype_hipp)-nnz(group==0)]               % excl cloz patients, 52=60-8
subtype_hipp = repmat(subtype_hipp,3,1);


[p_fitlme_grpwave_vol_ifg,t_fitlme_grpwave_vol_ifg, coef_fitlme_grpwave_vol_ifg, regional_vol_reg_ifg]=...
    S_func_fitlme_grpwave_global(tbl(subtype_ifg==1,:), regional_vol(subtype_ifg==1,:), global_vol(subtype_ifg==1));

[p_fitlme_grpwave_vol_hipp,t_fitlme_grpwave_vol_hipp, coef_fitlme_grpwave_vol_hipp, regional_vol_reg_hipp]=...
    S_func_fitlme_grpwave_global(tbl(subtype_hipp==1,:), regional_vol(subtype_hipp==1,:), global_vol(subtype_hipp==1));

% save('R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ_excl_cloz.mat','regional_*','tbl','*_fitlme_*','region_label','subtype_*')


%% 2.Brain plot of Group, Time, Group*Time effects
clc; clear; close all;
load('R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ.mat')
% load('R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ_excl_cloz.mat')

cmap = flipud(cbrewer('div','RdBu',256));
cmap(cmap>1)=1; cmap(cmap<0)=0;
clims = [-5.0,5.0];
pic_folder = 'R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ/'; % save brain plots to ...

%%%% Plot subtype, FDR sig %%%%
[h, crit_p]=fdr_bh(p_fitlme_grpwave_vol_hipp.MF_GROUP); nnz(h)           % crit_p=1.9905e-04, 3 sig
[h, crit_p]=fdr_bh(p_fitlme_grpwave_vol_hipp.("MF_GROUP:Period")); nnz(h)% crit_p=0.0156, 28 sig
region_label(fdr_bh(p_fitlme_grpwave_vol_hipp.MF_GROUP))
plotBrain(region_label, (t_fitlme_grpwave_vol_hipp.MF_GROUP).*fdr_bh(p_fitlme_grpwave_vol_hipp.MF_GROUP), cmap,'atlas','aparc_aseg','limits',clims,...
            'savePath',[pic_folder,'Hipp_Group_tstat_FDRsig'])

region_label(fdr_bh(p_fitlme_grpwave_vol_hipp.("MF_GROUP:Period")))
plotBrain(region_label, (t_fitlme_grpwave_vol_hipp.("MF_GROUP:Period")).*fdr_bh(p_fitlme_grpwave_vol_hipp.("MF_GROUP:Period")), cmap,'atlas','aparc_aseg','limits',clims,...
            'savePath',[pic_folder,'Hipp_GroupPeriod_tstat_FDRsig'])

[h, crit_p]=fdr_bh(p_fitlme_grpwave_vol_ifg.MF_GROUP); nnz(h)            % crit_p=0.0195, 35 sig
region_label(fdr_bh(p_fitlme_grpwave_vol_ifg.MF_GROUP))
plotBrain(region_label, (t_fitlme_grpwave_vol_ifg.MF_GROUP).*fdr_bh(p_fitlme_grpwave_vol_ifg.MF_GROUP), cmap,'atlas','aparc_aseg','limits',clims,...
            'savePath',[pic_folder,'IFG_Group_tstat_FDRsig'])
region_label(fdr_bh(p_fitlme_grpwave_vol_ifg.("MF_GROUP:Period")))       % none



%% 3.Plot change across time, ifg and hipp overlap regions
clc; clear; close all;
load('R3_sustain_subtype_new.mat', 'sustain_demo_T1','sustain_subtype')
load('R4_region_label_plot.mat')
load('R5_neuroimaging_regional_regress_alltime.mat', 'sustain_volume_T*_reg_ICV')
load('R7_neuroimaging_regional_fitlme_global_new_age.mat','tbl')

region_label = region_label_plot;

load('R7_neuroimaging_regional_fitlme_global_new_age_CON_SZ.mat', 'p_fitlme_grpwave_vol_*')
tmp1=union(find(fdr_bh(p_fitlme_grpwave_vol_hipp.("MF_GROUP:Period"))),find(fdr_bh(p_fitlme_grpwave_vol_hipp.("MF_GROUP"))));% 30
tmp2=find(fdr_bh(p_fitlme_grpwave_vol_ifg.MF_GROUP));
roi = intersect(tmp1,tmp2);

timepoint = tbl.Period;
group = tbl.MF_GROUP;
subtype = repmat(sustain_subtype.newdata_ml_subtype,3,1);
regional_vol_reg = [sustain_volume_T1_reg_ICV; sustain_volume_T2_reg_ICV; sustain_volume_T3_reg_ICV];

col = zeros(length(group),1);
col(group==0)=1;
col(group==3 & subtype==0)=2;
col(group==3 & subtype==1)=3;


cmap = [[0 0 0];[0,180, 255]; [255,127,0]; [152,78,163]]/255; 

close all
clear g
figure('Position',[50 50 1200 750]);
for i=1:length(roi)
    tmp_roi = regional_vol_reg(:,roi(i));
    % g(1,i)=gramm('x',timepoint, 'y',tmp_roi, 'color',col, 'subset',group~=1);   %%%% change here, subset %%%%
    % g(1,i).stat_smooth(); g(1,i).set_title(region_label(roi(i)));g(1,i).set_color_options('map',cmap);g(1,i).no_legend();
    % g(1,i).axe_property('YLim',[nanmean(tmp_roi)*0.7 nanmean(tmp_roi)*1.3], 'XTick',[1,2,3],'XTickLabel',{'Baseline','3 Year','6 Year'});

    if i<=5
        g(1,i)=gramm('x',timepoint, 'y',tmp_roi, 'color',col, 'subset',group~=1);
        g(1,i).stat_smooth(); g(1,i).set_title(region_label(roi(i)));g(1,i).set_color_options('map',cmap,'n_lightness',1);g(1,i).no_legend();
        g(1,i).axe_property('YLim',[floor(prctile(tmp_roi,15)/1000)*1000 ceil(prctile(tmp_roi,85)/1000)*1000], 'XTick',[1,2,3],'XTickLabel',{'Baseline','3 Year','6 Year'});
    elseif i>5 & i<=10
        g(2,i-5)=gramm('x',timepoint, 'y',tmp_roi, 'color',col, 'subset',group~=1);   %%%% change here, subset %%%%
        g(2,i-5).stat_smooth(); g(2,i-5).set_title(region_label(roi(i)));g(2,i-5).set_color_options('map',cmap);g(2,i-5).no_legend();
        g(2,i-5).axe_property('YLim',[floor(prctile(tmp_roi,15)/1000)*1000 ceil(prctile(tmp_roi,85)/1000)*1000], 'XTick',[1,2,3],'XTickLabel',{'Baseline','3 Year','6 Year'});
    else
        g(3,i-10)=gramm('x',timepoint, 'y',tmp_roi, 'color',col, 'subset',group~=1);   %%%% change here, subset %%%%
        g(3,i-10).stat_smooth(); g(3,i-10).set_title(region_label(roi(i))); g(3,i-10).set_color_options('map',cmap); g(3,i-10).no_legend();
        g(3,i-10).axe_property('YLim',[floor(prctile(tmp_roi,15)/1000)*1000 ceil(prctile(tmp_roi,85)/1000)*1000], 'XTick',[1,2,3],'XTickLabel',{'Baseline','3 Year','6 Year'});
    end
end

g.set_names('color','Group','x','','y','');
g.draw();
g.set_text_options('font','Arial','base_size',16)


pic_folder = 'E:\Projects\2022-GROUP\Scripts\R7_regional_fitlme_subtype_plot_global_new_age_CON_SZ/'; % save brain plots to ...
saveas(gcf,[pic_folder,'IFG_Hipp_overlap_1'],'tiffn')
