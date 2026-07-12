function spines_sorted_by_branch = get_spines_sorted_by_branch(tree, all_spine_ids, spine_type_for_branches)
% Input:
% Output:

if nargin < 3
    spine_type_for_branches = 4;
end

biological_sect = find_biological_branches(tree, spine_type_for_branches);

parent_paths_matrix = ipar_tree(tree);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
all_attachment_points = direct_parent_indices(all_spine_ids);

spines_sorted_by_branch = cell(size(biological_sect, 1), 1);

for i = 1:size(biological_sect, 1)
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    spines_on_branch_mask = ismember(all_attachment_points, branch_nodes);
    spine_ids_on_branch = all_spine_ids(spines_on_branch_mask);
    
    if ~isempty(spine_ids_on_branch)
        attachment_points_on_branch = all_attachment_points(spines_on_branch_mask);
        path_lengths_on_branch = all_path_lengths_from_root(attachment_points_on_branch);
        
        [~, sort_idx] = sort(path_lengths_on_branch);
        
        spines_sorted_by_branch{i} = spine_ids_on_branch(sort_idx);
    end
end

end