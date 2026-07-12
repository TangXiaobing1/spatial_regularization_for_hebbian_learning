function final_biological_sect = find_biological_branches(tree, spine_type, min_branch_length)

% Input:
% Output:


if nargin < 2, spine_type = 4; end
if nargin < 3, min_branch_length = 5; end

parent_paths_matrix = ipar_tree(tree);
all_segment_lengths = len_tree(tree);

[initial_sects] = find_all_potential_branches(tree, spine_type, parent_paths_matrix);
if isempty(initial_sects), final_biological_sect = []; return; end

branch_lengths = zeros(size(initial_sects, 1), 1);
for i = 1:size(initial_sects, 1)
    branch_lengths(i) = get_path_length(initial_sects(i,:), parent_paths_matrix, all_segment_lengths);
end

long_branches_mask = branch_lengths >= min_branch_length;
long_sects = initial_sects(long_branches_mask, :);

while true
    unique_nodes_in_long_branches = unique(long_sects(:));
    
    node_degrees = histcounts(long_sects(:), [unique_nodes_in_long_branches; max(unique_nodes_in_long_branches)+1]);
    
    nodes_to_merge = unique_nodes_in_long_branches(node_degrees == 2);
    
    if isempty(nodes_to_merge)
        break;
    end
    
    merged_in_this_iteration = false;
    for i = 1:length(nodes_to_merge)
        merge_node = nodes_to_merge(i);
        
        [row1, ~] = find(long_sects == merge_node);
        if length(row1) ~= 2, continue; end
        
        branch1 = long_sects(row1(1), :);
        branch2 = long_sects(row1(2), :);
        
        new_start_node = setdiff(branch1, merge_node);
        new_end_node = setdiff(branch2, merge_node);
        
        if isempty(new_start_node) || isempty(new_end_node), continue; end

        merged_branch = [new_start_node, new_end_node];
        
        long_sects(row1, :) = [];
        
        long_sects(end+1, :) = merged_branch;
        
        merged_in_this_iteration = true;
        break;
    end
    
    if ~merged_in_this_iteration
        break;
    end
end

final_biological_sect = long_sects;
disp(['分析完成。共找到 ', num2str(size(final_biological_sect, 1)), ' 条合并后的生物学分支路径。']);
end


% --- Helper function ---

function sects = find_all_potential_branches(tree, spine_type, parent_paths_matrix)
    is_true_branch_point = find_true_branch_points(tree, spine_type);
    is_true_termination_point = find_non_spine_termination_points(tree, spine_type);
    is_true_boundary = is_true_branch_point | is_true_termination_point;
    true_boundary_indices = find(is_true_boundary);
    direct_parent_indices = idpar_tree(tree);
    path_start_nodes = unique([1; find(is_true_branch_point)]);
    sects = [];
    for i = 1:length(path_start_nodes)
        current_start_node = path_start_nodes(i);
        child_indices = find(direct_parent_indices == current_start_node);
        non_spine_children = child_indices(tree.R(child_indices) ~= spine_type);
        for j = 1:length(non_spine_children)
            current_node = non_spine_children(j);
            while true
                if ismember(current_node, true_boundary_indices)
                    path_end_point = current_node; break;
                end
                next_child_indices = find(direct_parent_indices == current_node);
                next_non_spine_child = next_child_indices(tree.R(next_child_indices) ~= spine_type);
                if length(next_non_spine_child) == 1
                    current_node = next_non_spine_child;
                else
                    path_end_point = current_node; break;
                end
            end
            sects = [sects; current_start_node, path_end_point];
        end
    end
end

function len = get_path_length(sect, parent_paths_matrix, all_segment_lengths)
    start_node = sect(1); end_node = sect(2);
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node, 1);
    if isempty(start_idx_in_path)
        len = inf; return;
    end
    path_nodes = fliplr(path_to_root(1:start_idx_in_path));
    if length(path_nodes) < 2, len = 0; else
        len = sum(all_segment_lengths(path_nodes(2:end)));
    end
end