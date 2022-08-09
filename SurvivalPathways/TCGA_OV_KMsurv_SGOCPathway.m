%% Import gene expression and survival data
clear
clc

[gene_file, gene_path] = uigetfile({'*.txt'},'Select gene expression file');
[cna_file, cna_path] = uigetfile({'*.xlsx'},'Select 19p13 Loss file');
[S_file, S_path] = uigetfile({'*.xlsx';'*.xls';'*.csv';'*.txt'},'Select survivals file');

gene_data_raw = readtable([gene_path,'/',gene_file],'ReadRowNames',true);
OS_data = readtable([S_path,'/',S_file],'ReadRowNames',true);
cna_data = readtable([cna_path,'/',cna_file],'ReadRowNames',true);

gene_samples = intersect(gene_data_raw.Properties.VariableNames,cna_data.Properties.RowNames);
% Filter NaNs
genesum = sum(gene_data_raw{:,:},2);
rm_genes = isnan(genesum);

%% OPTION 1: Build patient cohorts according to ONLY gene expression LL/LH/HL/HH
genes_to_keep = {'UQCR11','MTHFD2'};
gene_data = gene_data_raw(:,genes_to_keep);
n_cohorts = 2^size(gene_data,2); % 2^n_Genes
gene_cutoffs = prctile(gene_data{:,:},50);
cohort_combi = bsxfun(@gt,gene_data{:,:},gene_cutoffs); 
cohort_grp = bi2de(double(cohort_combi));
gene_data(:,end+1) = table(cohort_grp);
figure; hist(cohort_grp,max(cohort_grp)+1);
%% OPTION 2: Build patient cohorts according to 19p13loss and MTHFD2 gene expression
genes_to_keep = {'MTHFD2'};
gene_data = gene_data_raw(gene_samples,genes_to_keep);

gene_cutoffs = prctile(gene_data{:,:},50);
cohort_combi = bsxfun(@gt,gene_data{:,:},gene_cutoffs); 
foo = [cna_data{gene_samples,1}, double(cohort_combi)];
cohort_grp = bi2de(foo);
gene_data(:,end+1) = table(cohort_grp);
figure; hist(cohort_grp,max(cohort_grp)+1);
% grps: 1 = loss 19p/MTHFD2 low, 3 = loss 19p/MTHFD2 High, 0 = Loss 19p
%% K-M for OS
clear OS PFS
matched_OS = intersect(gene_samples,OS_data.Properties.RowNames);

OS.data = OS_data(matched_OS,{'OverallSurvivalStatus','OverallSurvival_Months_'});
PFS.data = OS_data(matched_OS,{'ProgressionFreeStatus','ProgressFreeSurvival_Months_'});

OS.cens = OS.data{:,'OverallSurvivalStatus'};

OS.time = OS.data{:,'OverallSurvival_Months_'};
OS.group = gene_data{matched_OS,end};
OS.time2 = [OS.time; OS.time]; % Copying OSes of all patients and putting them in 4th group
OS.cens2 = [OS.cens; OS.cens]; % Copying OSes of all patients and putting them in 4th group
OS.group2 = [OS.group; repmat(max(OS.group)+1,size(OS.group))];

PFS.cens = PFS.data{:,'ProgressionFreeStatus'};

PFS.time = PFS.data{:,'ProgressFreeSurvival_Months_'};
PFS.group = gene_data{matched_OS,end};
PFS.time2 = [PFS.time; PFS.time]; % Copying OSes of all patients and putting them in 4th group
PFS.cens2 = [PFS.cens; PFS.cens]; % Copying OSes of all patients and putting them in 4th group
PFS.group2 = [PFS.group; repmat(max(PFS.group)+1,size(PFS.group))];

%% Plot K-M
figure; hold on
for i = unique(OS.group')
    disp(['Cohort #: ',num2str(i)])
    ecdf(OS.time(OS.group==i),'censoring',~OS.cens(OS.group==i),'function','survivor');

end
hold off
title('Overall Survival')
legend('Low UQCR11/Low MTHFD2', 'Low/High', 'High/Low', 'High/High')

figure; hold on
for i = unique(PFS.group')
    disp(['Cohort #: ',num2str(i)])
    ecdf(PFS.time(PFS.group==i),'censoring',~PFS.cens(PFS.group==i),'function','survivor');

end
hold off
title('Progression-Free Survival')
legend('Low UQCR11/Low MTHFD2', 'Low/High', 'High/Low', 'High/High')

%% Plot PFS K-M

[KM_PFS.p, KM_PFS.fh, KM_PFS.stats] = MatSurv(PFS.time,PFS.cens,PFS.group,'GroupsToUse',[0,2],'XStep',12); % Group 0 is low/low U/M; Group 1 is 

%% Plot OS K-M
[KM_OS.p, KM_OS.fh, KM_OS.stats] = MatSurv(OS.time,OS.cens,OS.group,'GroupsToUse',[0,2],'XStep',12);

%% Plot PFS K-M of Low/Low UQCR11/MTHFD2,Low/HighUQCR11/MTHFD2,All
[KM_PFS.p2, KM_PFS.fh2, KM_PFS.stats2] = MatSurv(PFS.time2,PFS.cens2,PFS.group2,'GroupsToUse',[1,3],'XStep',24);

%% Plot OS K-M of Low/Low UQCR11/MTHFD2,Low/HighUQCR11/MTHFD2,All
[KM_OS.p2, KM_OS.fh2, KM_OS.stats2] = MatSurv(OS.time2,OS.cens2,OS.group2,'GroupsToUse',[1,3],'XStep',24);

%% Updated 04-13-2021
% Pathway expression of SGOC related pathways in 19pDEL vs NONDEL
clear x_corr x_19p_corr p_genes_plot x_19p_ttest
gene_list = readtable('MetabGenes.xlsx', 'ReadVariableNames', true);
pathway_list = readtable('UserPathways.xlsx','ReadVariableNames', false);
cna_grp = cna_data{gene_samples,1};

for i = 1:size(pathway_list,1)
    p_genes = gene_list{:,pathway_list{i,1}{1}};
    p_genes_plot{i} = intersect(p_genes, gene_data.Properties.RowNames);
    
    x = gene_data{p_genes_plot{i}, gene_samples}';
    [x_corr{i}.r,x_corr{i}.p] = corr(x, 'type', 'Spearman');
    
    h = clustergram(x_corr{i}.r, 'Standardize', 'None', 'Colormap', redbluecmap,...
        'RowLabels', p_genes_plot{i}, 'ColumnLabels', p_genes_plot{i});
    addTitle(h, pathway_list{i,1}{1})

    for j = 1:size(x,2)
        [~,x_19p_ttest{i}.p(j,1)] = ttest2(x(cna_grp==0,j),x(cna_grp==1,j));
        x_19p_ttest{i}.padj = mafdr(x_19p_ttest{i}.p, 'BHFDR', true);
    end

end





