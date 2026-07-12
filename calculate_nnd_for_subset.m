function mean_nnd = calculate_nnd_for_subset(tree, subset_spine_indices, biological_sect)
% Input:
% biological_sect = find_biological_branches(tree);
% Output:

if isempty(subset_spine_indices) || length(subset_spine_indices) < 2
    mean_nnd = nan;
    return;
end

direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
parent_paths_matrix = ipar_tree(tree);
subset_attachment_points = direct_parent_indices(subset_spine_indices);

num_branches = size(biological_sect, 1);

all_isi_for_subset = [];
for i = 1:num_branches
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    spines_on_branch_mask = ismember(subset_attachment_points, branch_nodes);
    path_lengths_on_branch = all_path_lengths_from_root(subset_attachment_points(spines_on_branch_mask));
    
    if length(path_lengths_on_branch) >= 2
        sorted_lengths = sort(path_lengths_on_branch);
        branch_isi = diff(sorted_lengths);
        all_isi_for_subset = [all_isi_for_subset; branch_isi];
    end
end

if isempty(all_isi_for_subset)
    mean_nnd = nan;
else
    mean_nnd = mean(all_isi_for_subset);
end

end