function is_true_termination_point_vec = find_non_spine_termination_points(tree, spine_type)
% Input:
% Output:


if nargin < 2
    spine_type = 4;
end

num_nodes = size(tree.X, 1);
direct_parent_indices = idpar_tree(tree);
is_true_termination_point_vec = false(num_nodes, 1);

for i = 1:num_nodes
    if tree.R(i) == spine_type
        continue;
    end
    
    child_indices = find(direct_parent_indices == i);
    
    if isempty(child_indices)
        is_true_termination_point_vec(i) = true;
        continue;
    end
    
    child_types = tree.R(child_indices);
    non_spine_child_count = sum(child_types ~= spine_type);
    
    if non_spine_child_count == 0
        is_true_termination_point_vec(i) = true;
    end
end

end