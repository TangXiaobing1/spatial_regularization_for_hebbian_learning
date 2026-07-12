function is_true_branch_point_vec = find_true_branch_points(tree, spine_type)
% some kind like function dissect_tree in TREES toolbox, but optimized for neurons with spines
% Input:
% Output:

if ~isfield(tree, 'R') || ~isfield(tree, 'X')
    error('输入的 tree_k 结构体不完整，缺少 R 或坐标字段。');
end
if nargin < 2
    spine_type = 4;
end

num_nodes = size(tree.X, 1);
direct_parent_indices = idpar_tree(tree);

is_true_branch_point_vec = false(num_nodes, 1);

for i = 1:num_nodes
    child_indices = find(direct_parent_indices == i);
    
    if isempty(child_indices)
        continue;
    end
    
    child_types = tree.R(child_indices);
    
    non_spine_child_count = sum(child_types ~= spine_type);
    
    if non_spine_child_count >= 2
        is_true_branch_point_vec(i) = true;
    end
end

end
