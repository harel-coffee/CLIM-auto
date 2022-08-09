%% Clustering of CNA/mRNA data based on https://www.nature.com/articles/ng.3398#methods

%% Import Molecular Data
clc
clearvars
[file_cna, path_cna] = uigetfile('*.txt', 'Upload Omics Data, Feature x Sample format');
raw_data = readtable([path_cna,'/',file_cna],'ReadRowNames',1);

%% Reduce CNA by to Cytobands
user_reduce = 1;

if user_reduce
   cytob_all = raw_data{:,'Cytoband'};
   var_names = raw_data.Properties.VariableNames;
   var_names_bool = ones(size(var_names));
   var_names_bool(contains(var_names,'Cytoband')) = 0;
   var_names_bool(contains(var_names,'Gene')) = 0;
   cna_data = raw_data{:,logical(var_names_bool)};
   uniq_cytob = unique(cytob_all);
   new_data = zeros(length(uniq_cytob),size(cna_data,2));
   for i = 1:length(uniq_cytob)
      idx = ismember(cytob_all,uniq_cytob(i));
      new_data(i,:) = mean(cna_data(idx,:),1);
   end
   
   cna_data = array2table(new_data,...
       'RowNames',uniq_cytob,...
       'VariableNames',var_names(logical(var_names_bool)));
%    raw_data = raw_data2;
end

%% Filter "Adjacent Normal" Samples
sample_names = cna_data.Properties.VariableNames;
tcga_samples = contains(sample_names,'TCGA');
tcga_type = cellfun(@(x)x(14:15),sample_names(tcga_samples),'UniformOutput',0);
tcga_type = str2double(tcga_type);
tcga_tumor = tcga_type<10;

keep_idx = true(size(sample_names));
keep_idx(1:length(tcga_tumor)) = tcga_tumor;

cna_data = cna_data(:,keep_idx);
sample_names = cna_data.Properties.VariableNames;

%% Adding Mutation Data
sample_names_CNA = cellfun(@(s)tcgaID(s,15),sample_names,'UniformOutput',0);
cna_data.Properties.VariableNames = sample_names_CNA;

% Import Mutations File
[file_mut, path_mut] = uigetfile('*.txt', 'Upload Mutations Data, Feature x Sample format');
mut_data = readtable([path_mut,'/',file_mut],'ReadRowNames',1);
sample_names_mut = mut_data.Properties.VariableNames;
% Set the mutation score from default 1, to bottom 0.5th percentile, so after transformation 2^.-X, they score high
prctile(reshape(cna_data{:,:},1,numel(cna_data)),0.005);

%% Keep only common samples between CNA and mut datasets
keep_samples = intersect(sample_names_mut,sample_names_CNA);
conv_mut_data = mut_data{:,keep_samples};
conv_mut_data(logical(conv_mut_data)) = conv_mut_score;
keep_data = [cna_data{:,keep_samples};conv_mut_data];
keep_rownames = [cna_data.Properties.RowNames; mut_data.Properties.RowNames];

clust_data = array2table(keep_data,'VariableNames',keep_samples,'RowNames',keep_rownames);

%% Step 1: Randomize initial loading matrices 
% Unsupervised NMF was performed on a gene-by-sample matrix X first with 
% 20 randomly initialized instances of NMF using the MATLAB 
% (MathWorks, R2013a) multiplicative update NMF solver for 10 steps.

X = -clust_data{:,:}; % Put more weightage on deletions
X = 2.^X; % Convert to non-negative values
k_max = ceil(sqrt(size(clust_data,2)));
P = size(clust_data,2);
N = size(clust_data,1);
n_top_G = 100;
step1a_opt = statset('MaxIter',1000,'Display','final','UseParallel',1);
step1b_opt = statset('MaxIter',1000,'Display','final','TolFun',1e-6,'TolX',1e-6);
%% Run for multiple cluster tries 
for i = 20:2*k_max
    [W_1{i},H_1{i},ranked_1{i},D_1a(i),D_1b(i)] = ...
        rankedGnmf(X,i,n_top_G,step1a_opt);
    disp(['kth TRY # ', num2str(i), ' OF ', num2str(2*k_max)])
end

figure; plot(D_1a,'ro-')
hold on; plot(D_1b,'o-'); hold off
save_filename = strrep(file_cna,'_focal_data.txt','');
save([save_filename,'_Step1.mat'])

%% Step 2a: Bootstrapping 
clear W_bs H_bs ranked_Gbs ranked_tmp consensus_G
% Gene-by-factor consensus matrix populated by repeating Steps 1 & 2
% N_bootsrap times for ~80% of data and oversampling 5-fold

k = 22;

N_bootstrap = 200;
O_bootstrap = 5;
f_bootstrap = 0.8;

W_bs = cell(1,N_bootstrap);
H_bs = W_bs;
ranked_Gbs = W_bs;

for b = 1:N_bootstrap
    disp(['BOOTSTRAP # ', num2str(b), ' OF ', num2str(N_bootstrap)])
    samples_bootstrap = [];
    for i = 1:O_bootstrap
        samples_bootstrap = [samples_bootstrap, randi([1 P],1,floor(0.8*P))]; 
    end
    
    X_bootstrap = X(:,samples_bootstrap);
    
    [W_bs{b},H_bs{b},ranked_Gbs{b}] = ...
        rankedGnmf(X_bootstrap,k,n_top_G,step1a_opt);
end

save([save_filename,'_Step2.mat'])

% Step 2b: Consensus clustering of Bootstrapped results
consensus_G = zeros(N,N);
for b = 1:N_bootstrap
    ranked_tmp = ranked_Gbs{b};
    ranked_tmp = ranked_tmp(1:30,:);
    for i = 1:k
        clear linear_idx
        found_G_idx = nchoosek(ranked_tmp(:,i),2);
        linear_idx = sub2ind(size(consensus_G),found_G_idx(:,1),found_G_idx(:,2));
        linear_idx = [linear_idx...
            ;sub2ind(size(consensus_G),found_G_idx(:,2),found_G_idx(:,1))];
        consensus_G(linear_idx) = consensus_G(linear_idx) + 1;
    end
end

clustergram(consensus_G,'Standardize',2,'Cluster',2,'Colormap',redbluecmap)
save([save_filename,'_Step2.mat'])

% Step 4: CNA Features from consensus matrix G
% Rows are Genes/Cytobands, Columns are factors
% consensus_G(i,j) represents number of times gene i, occured as a top-500
% feature of factor j.
clear choose_c
Y = pdist(consensus_G,'correlation');
Z = linkage(Y,'complete');
n_cutoff = 200;
cutoff_c = linspace(1.15,1.16,n_cutoff);
for i = 1:n_cutoff
consensus_G_clust = cluster(Z, 'cutoff',cutoff_c(i));
c_k = length(unique(consensus_G_clust))
    if c_k<k+5&&c_k>k-5
       choose_c = cutoff_c(i);
       
    end
end
% Step 4b:
close all
consensus_G_clust = cluster(Z, 'cutoff',choose_c);
figure;hist(consensus_G_clust,k)
title('Consensus Gene Cluster Histogram')
G0 = 0.01*ones(N,k);
for i = 1:k
    idx = find(consensus_G_clust==i);
    G0(idx,i) = 1;
end

% Final NNMF with seeded G0
[W_final,H_final,D_final] = nnmf(X,k,'w0',G0,'algo','mult','opt',step1a_opt);
fooh = reshape(H_final,numel(H_final),1);
foow = reshape(W_final,numel(W_final),1);
figure; hist(fooh,100);
title('H values histogram')
figure; hist(foow,100);
title('W values histogram')

%% Plotting NNMF results
W_final_mean = mean(W_final,2);
[~,idx_top25] = maxk(W_final_mean,25);
W_top25 = clust_data.Properties.RowNames(idx_top25);

[W_meansort,W_sort_idx] = sortrows(W_final_mean,'descend');
Wz_sort = zscore(W_final,1);
Wz_sort = Wz_sort(W_sort_idx,:);
W_sortnames = clust_data.Properties.RowNames(W_sort_idx);

%% Save files
W_avgscore = table(W_sortnames,W_meansort,'VariableNames',{'Feature','Wscore'});