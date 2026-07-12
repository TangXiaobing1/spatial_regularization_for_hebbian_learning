function all_isi_values = calculate_isi(tree, spine_type)
% Stability and dynamics of dendritic spines in macaque prefrontal cortex
% Input:
% Output:

if nargin < 2
    spine_type = 4;
end

direct_parent_indices = idpar_tree(tree);
all_path_lengths = Pvec_tree(tree, len_tree(tree));
spine_indices = tree.R == spine_type;
spine_attachment_points = direct_parent_indices(spine_indices);

biological_sect = find_biological_branches(tree, spine_type);
num_branches = size(biological_sect, 1);
parent_paths_matrix = ipar_tree(tree);

all_isi_values = [];

for i = 1:num_branches
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    spines_on_this_branch_mask = ismember(spine_attachment_points, branch_nodes);
    
    path_lengths_of_spines_on_branch = all_path_lengths(spine_attachment_points(spines_on_this_branch_mask));
    
    if length(path_lengths_of_spines_on_branch) >= 2
        sorted_lengths = sort(path_lengths_of_spines_on_branch);
        
        branch_isi = diff(sorted_lengths);
        
        all_isi_values = [all_isi_values; branch_isi];
    end
end

disp(['计算完成。共找到 ', num2str(length(all_isi_values)), ' 个 Inter-Spine Intervals (ISIs)。']);

end