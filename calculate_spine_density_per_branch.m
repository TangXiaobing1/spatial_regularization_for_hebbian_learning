function branch_info = calculate_spine_density_per_branch(tree, spine_type)
% Input:
% Output:
% dependencies:function is_true_branch_point_vec = find_true_branch_points(tree, spine_type)

if ~isfield(tree, 'R') || ~isfield(tree, 'X')
    error('输入的 tree_k 结构体不完整，缺少 R 或坐标字段。');
end
if nargin < 2
    spine_type = 4;
end

[sect, ~] = dissect_tree(tree);
num_branches = size(sect, 1);

direct_parent_indices = idpar_tree(tree);
all_segment_lengths = len_tree(tree);
parent_paths_matrix = ipar_tree(tree);
is_branch_point_vec = find_true_branch_points(tree); 

branch_info = struct('branch_id', {},'branch_order', {}, ...
    'start_node', {}, 'end_node', {}, ...
                     'nodes_on_branch', {}, 'length', {}, ...
                     'spine_count', {}, 'density', {});

for i = 1:num_branches
    start_node = sect(i, 1);
    end_node = sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0); 
    
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path)
        warning(['分支 ', num2str(i), ' 的起始节点 ', num2str(start_node), ' 不在末端节点 ', num2str(end_node), ' 的路径上，跳过此分支。']);
        continue;
    end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    if length(branch_nodes) < 2
        branch_length = 0;
    else
        branch_segment_indices = branch_nodes(2:end);
        branch_length = sum(all_segment_lengths(branch_segment_indices));
    end
    
    spine_count_on_branch = 0;
    for j = 1:length(branch_nodes)
        current_branch_node = branch_nodes(j);
        attached_spines = find(direct_parent_indices == current_branch_node & tree.R == spine_type);
        spine_count_on_branch = spine_count_on_branch + length(attached_spines);
    end
    
    if branch_length > 0
        density_on_branch = spine_count_on_branch / branch_length;
    else
        density_on_branch = 0; 
    end

    branch_points_on_path_count = sum(is_branch_point_vec(path_to_root));
    
    branch_info(i).branch_id = i;
    branch_info(i).branch_order = branch_points_on_path_count;
    branch_info(i).start_node = start_node;
    branch_info(i).end_node = end_node;
    branch_info(i).nodes_on_branch = branch_nodes;
    branch_info(i).length = branch_length;
    branch_info(i).spine_count = spine_count_on_branch;
    branch_info(i).density = density_on_branch;
    
end


end