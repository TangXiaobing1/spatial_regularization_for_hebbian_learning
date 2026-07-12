S = load('microns_example_all.mat');
% make sure this file is in your path

k =138;
tree_k = S.out_trees{k};
spine_indices = find(tree_k.R == 4);  % spine index
len = len_tree(tree_k);
spine_lengths = len(spine_indices);


% remove spine with lenth > 4μm
long_spine_indices = spine_indices(spine_lengths > 4);
if ~isempty(long_spine_indices)
    tree_k_corrected = delete_tree(tree_k, long_spine_indices);
    disp(['Removed ', num2str(length(long_spine_indices)), ' spines longer than 4 microns.']);
else
    disp('No spines longer than 4 microns found.');
end

spine_indices = find(tree_k_corrected.R == 4);  % spine index
len = len_tree(tree_k_corrected);
spine_lengths = len(spine_indices);
spine_diameters = tree_k_corrected.D(spine_indices);% Computespinediameter


% spine classification
% Based on Harris & Stevens (1989)
mushroom = spine_diameters > 0.6 & spine_lengths < 0.8;
filiform = spine_diameters < 0.3 & spine_lengths > 1;
stubby = spine_diameters > 0.4 & spine_lengths < 0.5;

tree_k_mushroom = delete_tree(tree_k_corrected, spine_indices(mushroom==0));
tree_k_filiform = delete_tree(tree_k_corrected, spine_indices(filiform==0));
tree_k_stubby = delete_tree(tree_k_corrected, spine_indices(stubby==0));

% Visualization

histogram(tree_k_corrected.D(spine_indices));
tree_k_dendrite_arbor = delete_tree(tree_k,spine_indices);
tree_k_dendrite = repair_tree(tree_k_dendrite_arbor);
plot_tree(tree_k_dendrite);
title(['Tree ',num2str(k),' without spine']);
plot_tree(tree_k);
title('Tree',num2str(k));

[p_k_s,p_sd_k_s,leaf_N_k_s,leaf_num_k_s,error_N_k_s,b_k_s]=perfection_index(tree_k_corrected);
error_N_k_s=[]; % Error bars are not interpreted here
plot_subtree_distribution(leaf_N_k_s, leaf_num_k_s, error_N_k_s, b_k_s, p_k_s, p_sd_k_s);
title(['Subtree Distribution of Tree ',num2str(k),' with spine']);
[p_k,p_sd_k,leaf_N_k,leaf_num_k,error_N_k,b_k]=perfection_index(tree_k_dendrite);
error_N_k=[]; % Error bars are not interpreted here
plot_subtree_distribution(leaf_N_k, leaf_num_k, error_N_k, b_k, p_k, p_sd_k);
title(['Subtree Distribution of Tree ',num2str(k),' without spine']);
[p_k_m,p_sd_k_m,leaf_N_k_m,leaf_num_k_m,error_N_k_m,b_k_m]=perfection_index(tree_k_mushroom);
error_N_k_m=[]; % Error bars are not interpreted here
plot_subtree_distribution(leaf_N_k_m, leaf_num_k_m, error_N_k_m, b_k_m, p_k_m, p_sd_k_m);
title(['Subtree Distribution of Tree ', num2str(k), ' with mushroom type spine']);
[p_k_f,p_sd_k_f,leaf_N_k_f,leaf_num_k_f,error_N_k_f,b_k_f]=perfection_index(tree_k_filiform);
error_N_k_f=[]; % Error bars are not interpreted here
plot_subtree_distribution(leaf_N_k_f, leaf_num_k_f, error_N_k_f, b_k_f, p_k_f, p_sd_k_f);
title(['Subtree Distribution of Tree ', num2str(k), ' with filiform type spine']);
[p_k_st,p_sd_k_st,leaf_N_k_st,leaf_num_k_st,error_N_k_st,b_k_st]=perfection_index(tree_k_stubby);
error_N_k_st=[]; % Error bars are not interpreted here
plot_subtree_distribution(leaf_N_k_st, leaf_num_k_st, error_N_k_st, b_k_st, p_k_st, p_sd_k_st);
title(['Subtree Distribution of Tree ', num2str(k), ' with stubby type spine']);
