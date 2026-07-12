function [spine_head_volumes, spine_head_volumes_zscore, spine_base_ids_in_order] = extract_spine_head_volume_feature(tree, spine_type)
% Input:
% Output:

if nargin < 2
    spine_type = 4;
end

is_terminal_vec = T_tree(tree);
is_spine_vec = (tree.R == spine_type);
spine_head_indices = find(is_terminal_vec & is_spine_vec);

num_spines = length(spine_head_indices);
if num_spines == 0
    spine_head_volumes = [];
    spine_head_volumes_zscore = [];
    spine_base_ids_in_order = [];
    return;
end

spine_base_ids_in_order = zeros(num_spines, 1);
direct_parent_indices = idpar_tree(tree);
all_volumes_vec = vol_tree(tree);
spine_head_volumes = zeros(num_spines, 1);

for i = 1:num_spines
    current_spine_head = spine_head_indices(i);
    current_node = current_spine_head;
    
    current_spine_path_nodes = [];
    attachment_point_on_backbone = -1;
    while true
        if tree.R(current_node) ~= spine_type
            attachment_point_on_backbone = current_node;
            break;
        end
        current_spine_path_nodes = [current_node, current_spine_path_nodes];
        parent_node = direct_parent_indices(current_node);
        if parent_node == current_node || parent_node == 0
            break;
        end
        current_node = parent_node;
    end
    
    if ~isempty(current_spine_path_nodes)
        spine_base_ids_in_order(i) = current_spine_path_nodes(1);
    else
        spine_base_ids_in_order(i) = current_spine_head;
    end
    
    spine_head_volumes(i) = sum(all_volumes_vec(current_spine_path_nodes));
end

%% 3. Z-scorestandardized
if num_spines > 1 && std(spine_head_volumes) > 0
    spine_head_volumes_zscore = zscore(spine_head_volumes);
else
    spine_head_volumes_zscore = zeros(size(spine_head_volumes));
end

end
