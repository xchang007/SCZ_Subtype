%% 1.ENIGMA geneexp of BrainSeq DEG gene lists 

clc;clear;close all;
load('R18_DEG_gene_exp.mat')
load('E:\Matlab\Toolbox_Fudan\ENIGMA-2.0.0\R_practice_ENIGMA_toolbox.mat')

gene_symbol = genes_10.Properties.VariableNames(2:end)';  % 15634 genes
GG = table2array(genes_10(reg_reorder_idx, 2:end));       % reorder brain regions
GG = GG./ nanstd(GG);

[tf,idx_gs] = ismember(gene_symbol, DEG_genes_Caudate.Symbol);
nnz(tf)  % 2175/2495 genes found in gene_symbol
GE_ENIGMA_mean_Caudate =  nanmean(GG(:, tf), 2);               
genes = DEG_genes_Caudate.Symbol;
tf_up = ismember(gene_symbol, genes(DEG_genes_Caudate.t>0)); nnz(tf_up)   % 1029
tf_down =ismember(gene_symbol,genes(DEG_genes_Caudate.t<0)); nnz(tf_down) % 1146
GE_ENIGMA_mean_Caudate_up =  nanmean(GG(:, tf_up), 2);  
GE_ENIGMA_mean_Caudate_down =nanmean(GG(:, tf_down), 2);  

% save('R18_DEG_gene_exp.mat','GE_ENIGMA_mean_*','-append')

%%%% brain plot DEG gene exp %%%%
clc;close all
addpath(genpath('E:/Matlab/Toolbox_Fudan/'));
pic_folder = 'E:\Projects\2022-GROUP\Scripts\S16_BrainSeq\DEG_brainplot';

cort_idx = 15:82;
subcort_idx = [7	6	2	5	4	3	1	14	13	9	12	11	10	8];

GE_plot = GE_ENIGMA_mean_Caudate_up;
plot_fsa5 = parcel_to_surface(GE_plot(cort_idx), 'aparc_fsa5');
f = figure,
plot_cortical(plot_fsa5, 'surface_name', 'fsa5', 'cmap', 'Reds',...
    'color_range', [3.5 5.5],  'label_text', ['DEG up Cortical']);
saveas(gcf, fullfile(pic_folder,['DEG up Cortical']), 'tiffn');

f = figure,
plot_subcortical(GE_plot(subcort_idx), 'ventricles', 'False', 'cmap', 'Reds',...
    'color_range', [3.5 5.5],  'label_text',  ['DEG up Subcortical']);
saveas(gcf, fullfile(pic_folder,'DEG up Subcortical'), 'tiffn');

%% 2. DEG gene exp corr with receptor module

clc; clear; close all;
cd('E:\Projects\2022-GROUP\Scripts')
load('R3_sustain_subtype_new.mat','sustain_demo_T1','sustain_subtype')
load('R14_PET_maps.mat')
load('R18_DEG_gene_exp.mat','GE_ENIGMA_mean_Caudate_*')

cort_idx = 15:82;
subcort_idx = 1:14;


%%%% corr with receptor_data by module %%%%
Nmodule = 4;
Nroi = size(receptor_data,1);
receptor_module = nan(Nroi, Nmodule);
for m= 1:Nmodule
    %%%% subcortical %%%%
    pet_idx = PET_corr_cluster_subcort.OriginalIndex(PET_corr_cluster_subcort.ClusterLabel==m);
    receptor_names(pet_idx)
    receptor_module(subcort_idx,m) = mean(table2array(receptor_data(subcort_idx,pet_idx)),2);

    pet_idx = PET_corr_cluster_cort.OriginalIndex(PET_corr_cluster_cort.ClusterLabel==m);
    receptor_names(pet_idx)
    receptor_module(cort_idx,m)  = mean(table2array(receptor_data(cort_idx,pet_idx)),2);
end

r_DEG_PET_corr = nan(4,Nmodule); % 4rows: subc_up, cort_up, subc_down, cort_down
p_DEG_PET_corr = nan(4,Nmodule); % 4rows: subc_up, cort_up, subc_down, cort_down

figure; [r,p]= corrplot([GE_ENIGMA_mean_Caudate_up(subcort_idx), receptor_module(subcort_idx,:)],'testR','on'); title('DEG up subc & receptor module')
r_DEG_PET_corr(1,:) = r(1,2:5);   p_DEG_PET_corr(1,:) = p(1,2:5);
figure; [r,p]= corrplot([GE_ENIGMA_mean_Caudate_up(cort_idx),    receptor_module(cort_idx,:)],'testR','on'); title('DEG up cort & receptor module')
r_DEG_PET_corr(2,:) = r(1,2:5);   p_DEG_PET_corr(2,:) = p(1,2:5);
figure; [r,p]= corrplot([GE_ENIGMA_mean_Caudate_down(subcort_idx), receptor_module(subcort_idx,:)],'testR','on'); title('DEG down subc & receptor module')
r_DEG_PET_corr(3,:) = r(1,2:5);   p_DEG_PET_corr(3,:) = p(1,2:5);
figure; [r,p]= corrplot([GE_ENIGMA_mean_Caudate_down(cort_idx),    receptor_module(cort_idx,:)],'testR','on'); title('DEG down cort & receptor module')
r_DEG_PET_corr(4,:) = r(1,2:5);   p_DEG_PET_corr(4,:) = p(1,2:5);

r_DEG_PET_corr = array2table(r_DEG_PET_corr,'VariableNames',{'module1','module2','module3','module4'},'RowNames',{'subc_up', 'cort_up', 'subc_down','cort_down'});
p_DEG_PET_corr = array2table(p_DEG_PET_corr,'VariableNames',{'module1','module2','module3','module4'},'RowNames',{'subc_up', 'cort_up', 'subc_down','cort_down'});

[h_fdr, crit_p]=fdr_bh(table2array(p_DEG_PET_corr)) % 1, crit_p = 3.3461e-04

%%%%% spin test using ENIGMA toolbox %%%%%
Nperm = 5000;

% DEG up_cort & module 2 
[p_spin_cort, r_null_spin_cort] = spin_test(GE_ENIGMA_mean_Caudate_up(cort_idx), receptor_module(cort_idx, 2), 'surface_name', 'fsa5', ...
    'parcellation_name', 'aparc', 'n_rot', Nperm, 'type', 'pearson');
% p_spin_cort = 0.0387

% DEG down_subc & module 3 
[p_spin_subcort, r_null_spin_subcort] = shuf_test(GE_ENIGMA_mean_Caudate_down(subcort_idx), receptor_module(subcort_idx,3),...
    'n_rot', Nperm, 'type', 'pearson');
% p_spin_subcort = 3.0000e-04


%%%% Plot for paper: wholebrain t-maps correlation %%%%
pic_folder = 'E:\Projects\2022-GROUP\Scripts\S16_BrainSeq\DEG_PET_corr';

% DEG up_cort & module 2 
DEG_PET_plot = [[0:length(cort_idx)-1]', GE_ENIGMA_mean_Caudate_up(cort_idx), receptor_module(cort_idx, 2)];
DEG_PET_plot= array2table(DEG_PET_plot,'VariableNames',{'sample','DEG','PET'});
writetable(DEG_PET_plot, fullfile(pic_folder,'DEG_PET_plot_cort.tsv'), 'FileType', 'text','Delimiter', '\t', 'WriteRowNames', true);
% DEG down_subc & module 3 
DEG_PET_plot = [[0:length(subcort_idx)-1]', GE_ENIGMA_mean_Caudate_down(subcort_idx), receptor_module(subcort_idx,3)];
DEG_PET_plot= array2table(DEG_PET_plot,'VariableNames',{'sample','DEG','PET'});
writetable(DEG_PET_plot, fullfile(pic_folder,'DEG_PET_plot_subcort.tsv'), 'FileType', 'text','Delimiter', '\t', 'WriteRowNames', true);


%% 3.DEG gene exp corr with zscore
clc; clear; close all;
load('R3_sustain_subtype_new.mat','sustain_demo_T1','sustain_subtype')
load('R18_DEG_gene_exp.mat')

%%%%% group-avg z-score %%%%% 
load('R5_neuroimaging_regional_zscore_alltime.mat')
group = sustain_demo_T1.MF_GROUP;
subtype = sustain_subtype.newdata_ml_subtype;

cort_idx = 15:82;
subcort_idx = 1:14;

group_CON_subtype = group;                    % 0:con, 1:sib, 3:pat
group_CON_subtype(group==3 & subtype==0) = 2; % 0:con, 1:sib, 2:pat-IFG, 3:pat-hipp
 
volume_zscore_grp = grpstats(sustain_volume_T1_zscore,group_CON_subtype,'mean');   %  grpstats ignores missing values

load('R18_DEG_gene_exp.mat','GE_ENIGMA_mean_*')
cort_idx = 15:82;
subcort_idx = 1:14;
Nperm=1000;

GE_maps = [GE_ENIGMA_mean_Caudate_up,GE_ENIGMA_mean_Caudate_down,...
    GE_ENIGMA_mean_DLPFC_up,GE_ENIGMA_mean_DLPFC_down,...
    GE_ENIGMA_mean_HIPPO_up,GE_ENIGMA_mean_HIPPO_down ];
figure;corrplot([volume_zscore_grp(3:4,cort_idx)',GE_maps(cort_idx,:)],'rows','pairwise','testR','on');
figure;corrplot([volume_zscore_grp(3:4,subcort_idx)',GE_maps(subcort_idx,:)],'rows','pairwise','testR','on');
[p, r_null] = shuf_test(volume_zscore_grp(4,subcort_idx)', GE_maps(subcort_idx,2)',...
        'n_rot', Nperm, 'type', 'pearson');


%%%%% individual z-score %%%%% 
load('R3_sustain_subtype_new.mat','sustain_demo_T1','sustain_subtype')
load('R5_neuroimaging_regional_zscore_alltime.mat')
group = sustain_demo_T1.MF_GROUP;
subtype = sustain_subtype.newdata_ml_subtype;

cort_idx = 15:82;
subcort_idx = 1:14;

group_CON_subtype = group;                    % 0:con, 1:sib, 3:pat
group_CON_subtype(group==3 & subtype==0) = 2; % 0:con, 1:sib, 2:pat-IFG, 3:pat-hipp 
volume_zscore_individual = sustain_volume_T1_zscore;   %  select T1 T2 T3

[r_DEG_cort_individual(:,1),   p_DEG_cort_individual(:,1)]=corr(GE_ENIGMA_mean_Caudate_up(cort_idx),  volume_zscore_individual(:,cort_idx)','rows','pairwise');
[r_DEG_cort_individual(:,2),   p_DEG_cort_individual(:,2)]=corr(GE_ENIGMA_mean_Caudate_down(cort_idx),volume_zscore_individual(:,cort_idx)','rows','pairwise');
[r_DEG_subcort_individual(:,1),p_DEG_subcort_individual(:,1)]=corr(GE_ENIGMA_mean_Caudate_up(subcort_idx),volume_zscore_individual(:,subcort_idx)','rows','pairwise');
[r_DEG_subcort_individual(:,2),p_DEG_subcort_individual(:,2)]=corr(GE_ENIGMA_mean_Caudate_down(subcort_idx),volume_zscore_individual(:,subcort_idx)','rows','pairwise');

r_DEG_cort_mean = grpstats(r_DEG_cort_individual,group_CON_subtype,'mean');
r_DEG_subcort_mean = grpstats(r_DEG_subcort_individual,group_CON_subtype,'mean');

z_DEG_cort_individual = 0.5*log((1+r_DEG_cort_individual)./(1-r_DEG_cort_individual));
z_DEG_subcort_individual = 0.5*log((1+r_DEG_subcort_individual)./(1-r_DEG_subcort_individual));


%%%%% compare r_PET between controls and subtypes %%%%%
NPET = 2;                     % 1st col Caudate_up, 2nd col Caudate_down
p_group_cort = nan(NPET,4);   % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
t_group_cort = nan(NPET,4);   % anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
p_group_subcort = nan(NPET,4);% anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp
t_group_subcort = nan(NPET,4);% anova, con vs pat-IFG, con vs pat-hipp, pat-IFG vs pat-hipp

for i=1:NPET
    %%%% cortical %%%%
    [p,tbl]   = anova1(z_DEG_cort_individual(group_CON_subtype~=1,i), group_CON_subtype(group_CON_subtype~=1),'off');
    p_group_cort(i,1) = p; t_group_cort(i,1)   = tbl{2,5};

    [~,p,~,stats] = ttest2(z_DEG_cort_individual(group_CON_subtype==0,i), z_DEG_cort_individual(group_CON_subtype==2,i));
    p_group_cort(i,2) = p; t_group_cort(i,2)   = stats.tstat;
    [~,p,~,stats] = ttest2(z_DEG_cort_individual(group_CON_subtype==0,i), z_DEG_cort_individual(group_CON_subtype==3,i));
    p_group_cort(i,3) = p; t_group_cort(i,3)   = stats.tstat;
    [~,p,~,stats] = ttest2(z_DEG_cort_individual(group_CON_subtype==2,i), z_DEG_cort_individual(group_CON_subtype==3,i));
    p_group_cort(i,4) = p; t_group_cort(i,4)   = stats.tstat;

    %%%% subcortical %%%%
    [p,tbl]   = anova1(z_DEG_subcort_individual(group_CON_subtype~=1,i),group_CON_subtype(group_CON_subtype~=1),'off');
    p_group_subcort(i,1)= p; t_group_subcort(i,1)   = tbl{2,5};

    [~,p,~,stats] = ttest2(z_DEG_subcort_individual(group_CON_subtype==0,i), z_DEG_subcort_individual(group_CON_subtype==2,i));
    p_group_subcort(i,2) = p; t_group_subcort(i,2)   = stats.tstat;
    [~,p,~,stats] = ttest2(z_DEG_subcort_individual(group_CON_subtype==0,i), z_DEG_subcort_individual(group_CON_subtype==3,i));
    p_group_subcort(i,3) = p; t_group_subcort(i,3)   = stats.tstat;
    [~,p,~,stats] = ttest2(z_DEG_subcort_individual(group_CON_subtype==2,i), z_DEG_subcort_individual(group_CON_subtype==3,i));
    p_group_subcort(i,4) = p; t_group_subcort(i,4)   = stats.tstat;
end


p_group_cort = array2table(p_group_cort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
t_group_cort = array2table(t_group_cort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
p_group_subcort = array2table(p_group_subcort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});
t_group_subcort = array2table(t_group_subcort,'VariableNames',{'F_anova','t_con_S1','t_con_S2','t_S1_S2'});

[h, crit_p]=fdr_bh([p_group_cort.F_anova,p_group_subcort.F_anova])% crit_p = 0.0166
sum(h) % 9  4 out of 19 receptors

