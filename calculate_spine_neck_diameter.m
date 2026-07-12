function neck_diameters = calculate_spine_neck_diameter(tree, spine_type)
% Input:
% Output:
% Dependencies：TREES toolbox

if nargin < 2
    spine_type = 4;
end

spine_head_indices = find(tree.R == spine_type);
num_spines = length(spine_head_indices);

if num_spines == 0
    disp('未在树中找到任何树突棘。');
    neck_diameters = [];
    return;
end

direct_parent_indices = idpar_tree(tree);

neck_diameters = zeros(num_spines, 1);

for i = 1:num_spines
    current_spine_head = spine_head_indices(i);
    current_node = current_spine_head;
    
    current_spine_path = [];
    
    while true
        if tree.R(current_node) ~= spine_type
            break;
        end
        
        current_spine_path = [current_spine_path, current_node]; % #ok<AGROW>
        
        parent_node = direct_parent_indices(current_node);
        
        if parent_node == current_node || parent_node == 0
            break;
        end
        
        current_node = parent_node;
    end
    
    if ~isempty(current_spine_path)
        current_spine_path = [current_spine_path, current_node];
    end

    if isempty(current_spine_path)
        neck_diameters(i) = tree.D(current_spine_head);
    else
        path_diameters = tree.D(current_spine_path);
        neck_diameters(i) = min(path_diameters);
    end
end

end