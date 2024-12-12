%% 1.corrplot receptor data, hierarchical clustering
clear; clc; close all;
cd E:\Projects\2022-GROUP\Scripts
addpath(genpath('E:/Matlab/Toolbox_Fudan/'));
load('R14_PET_maps.mat')

[r_PET_corr_subcort,p_PET_corr_subcort]=corr(table2array(receptor_data(1:14,:)));
[r_PET_corr_cort,p_PET_corr_cort]=corr(table2array(receptor_data(15:end,:)));

cmap = flipud(cbrewer('div','RdBu',128)); cmap(cmap<0)=0; cmap(cmap>1)=1;
figure; set(gcf,'Position',[100,100,800,800])
h=heatmap(r_PET_corr_subcort,'Colormap',cmap,'ColorLimits',[-1.5,1.5],'CellLabelFormat','%.2f');
h.XDisplayLabels = receptor_names;
h.YDisplayLabels = receptor_names;

%%%% chiplot clustering %%%%
pic_folder = 'R14_PET_maps\PET_corrplot';
r_tbl = array2table(r_PET_corr_subcort,'VariableNames',receptor_names, 'RowNames',receptor_names);
writetable(r_tbl, fullfile(pic_folder,'r_PET_corr_subcort.tsv'), 'FileType', 'text','Delimiter', '\t', 'WriteRowNames', true);


%% 2.1 plot r_PET between subtypes, reorder by module, horizontal
clear; clc; close all;
load('R3_sustain_subtype_new.mat')
load('R14_PET_maps.mat')
load('R16_PET_maps_corr_zscore_individual.mat')

group = sustain_demo_T1.MF_GROUP;
subtype = sustain_subtype.newdata_ml_subtype;

r_PET_cort = r_PET_cort_individual_T3;       % change to T2,T3 here!
r_PET_subcort = r_PET_subcort_individual_T3; % change to T2,T3 here!


%%%% chiplot heatmap, reorder by module, horizontal %%%%
r_subtype_cort = [nanmean(r_PET_cort(group==3 & subtype==0,:))',...
    nanmean(r_PET_cort(group==3 & subtype==1,:))'];
r_subtype_subcort = [nanmean(r_PET_subcort(group==3 & subtype==0,:))',...
    nanmean(r_PET_subcort(group==3 & subtype==1,:))'];

r_subtype_cort = r_subtype_cort(PET_corr_cluster_cort.OriginalIndex,:);
r_subtype_subcort = r_subtype_subcort(PET_corr_cluster_subcort.OriginalIndex,:);

r_subtype_cort = array2table(r_subtype_cort','VariableNames',receptor_names(PET_corr_cluster_cort.OriginalIndex),...
    'RowNames',{'SCZ S1','SCZ S2'});
r_subtype_subcort = array2table(r_subtype_subcort','VariableNames',receptor_names(PET_corr_cluster_subcort.OriginalIndex),...
    'RowNames',{'SCZ S1','SCZ S2'});

pic_folder = 'Figures\Fig3_PET_fitlme\PET_corr_individual\r_SCZ_subtype';
writetable(r_subtype_cort, fullfile(pic_folder,'r_subtype_cort_T3_horz.tsv'), 'FileType', 'text','Delimiter', '\t', 'WriteRowNames', true);
writetable(r_subtype_subcort, fullfile(pic_folder,'r_subtype_subcort_T3_horz.tsv'), 'FileType', 'text','Delimiter', '\t', 'WriteRowNames', true);




%% 2.2 compare r_PET between controls and subtypes,  by module

clear; clc; close all;
load('R3_sustain_subtype_new.mat','sustain_demo_T1','sustain_subtype')
load('R14_PET_maps.mat')
load('R16_PET_maps_corr_zscore_individual.mat','r_PET_*cort_individual*')

covar = [sustain_demo_T1.MF_AGE, sustain_demo_T1.MF_AGE.^2, sustain_demo_T1.MF_GENDER, sustain_demo_T1.MF_ETHN2];

group = sustain_demo_T1.MF_GROUP;
subtype = sustain_subtype.newdata_ml_subtype;

group_CON_subtype = group;                    % 0:con, 1:sib, 3:pat
group_CON_subtype(group==3 & subtype==0) = 2; % 0:con, 1:sib, 2:pat-IFG, 3:pat-hipp
[nnz(group_CON_subtype==0),nnz(group_CON_subtype==1),nnz(group_CON_subtype==2),nnz(group_CON_subtype==3)]
% 136   202   114    60


r_PET_subcort = r_PET_subcort_individual; %%%%% choose here: T1,T2,T3 %%%%%
r_PET_cort = r_PET_cort_individual;       %%%%% choose here: T1,T2,T3 %%%%%

r_PET_subcort_reg = my_out_zscore_covar(r_PET_subcort,[],[], covar)+nanmean(r_PET_subcort);
r_PET_cort_reg = my_out_zscore_covar(r_PET_cort,[],[], covar)+nanmean(r_PET_cort);

z_PET_subcort = 0.5*log((1+r_PET_subcort_reg)./(1-r_PET_subcort_reg));
z_PET_cort = 0.5*log((1+r_PET_cort_reg)./(1-r_PET_cort_reg));

Nsub = size(sustain_demo_T1,1);
Nmodule = 4;

p_group_cort = nan(Nmodule,4);          % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
t_group_cort = nan(Nmodule,4);          % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
p_group_subcort = nan(Nmodule,4);       % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
t_group_subcort = nan(Nmodule,4);       % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
r_PET_module_cort=nan(Nsub, Nmodule);
r_PET_module_subcort=nan(Nsub, Nmodule);

for m=1:Nmodule
    %%%% cortical %%%%
    pet_idx = PET_corr_cluster_cort_dist.OriginalIndex(PET_corr_cluster_cort_dist.ClusterLabel==m);
    r_PET_module_cort(:,m) = nanmean(r_PET_cort(:,pet_idx),2);

    tmp_z = nanmean(z_PET_cort(:,pet_idx),2);
    [p,tbl]   = anova1(tmp_z(group_CON_subtype~=1), group_CON_subtype(group_CON_subtype~=1),'off');
    p_group_cort(m,1) = p; t_group_cort(m,1)   = tbl{2,5};

    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==0), tmp_z(group_CON_subtype==2));
    p_group_cort(m,2) = p; t_group_cort(m,2)   = stats.tstat;
    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==0), tmp_z(group_CON_subtype==3));
    p_group_cort(m,3) = p; t_group_cort(m,3)   = stats.tstat;
    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==2), tmp_z(group_CON_subtype==3));
    p_group_cort(m,4) = p; t_group_cort(m,4)   = stats.tstat;

    %%%% subcortical %%%%
    pet_idx = PET_corr_cluster_subcort_dist.OriginalIndex(PET_corr_cluster_subcort_dist.ClusterLabel==m);
    r_PET_module_subcort(:,m) = nanmean(r_PET_subcort(:,pet_idx),2);
   
    tmp_z = nanmean(z_PET_subcort(:,pet_idx),2);
    [p,tbl]   = anova1(tmp_z(group_CON_subtype~=1),group_CON_subtype(group_CON_subtype~=1),'off');
    p_group_subcort(m,1)= p; t_group_subcort(m,1)   = tbl{2,5};

    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==0), tmp_z(group_CON_subtype==2));
    p_group_subcort(m,2) = p; t_group_subcort(m,2)   = stats.tstat;
    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==0), tmp_z(group_CON_subtype==3));
    p_group_subcort(m,3) = p; t_group_subcort(m,3)   = stats.tstat;
    [~,p,~,stats] = ttest2(tmp_z(group_CON_subtype==2), tmp_z(group_CON_subtype==3));
    p_group_subcort(m,4) = p; t_group_subcort(m,4)   = stats.tstat;
end

p_group_cort_T1 = array2table(p_group_cort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
t_group_cort_T1 = array2table(t_group_cort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
p_group_subcort_T1 = array2table(p_group_subcort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
t_group_subcort_T1 = array2table(t_group_subcort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});

r_PET_module_cort_T1 = r_PET_module_cort;
r_PET_module_subcort_T1 = r_PET_module_subcort;

% save('R16_PET_maps_corr_zscore_module.mat','*_group_*cort_*','r_PET_module_*', '-append')

