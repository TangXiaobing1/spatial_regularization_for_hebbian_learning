function [head_volumes_terminal, spine_head_ids, spine_volumes_total, dist_from_root] = calculate_spine_head_volumes(tree, spine_type)
% Input:
% Output:

if nargin < 2
    spine_type = 4;
end

all_vols = vol_tree(tree);
all_dists = Pvec_tree(tree);
is_terminal = T_tree(tree);
direct_parent = idpar_tree(tree);

% --- 2. Identify Spine Heads ---
is_spine = (tree.R == spine_type);

spine_head_ids = find(is_terminal & is_spine);

num_spines = length(spine_head_ids);
if num_spines == 0
    head_volumes_terminal = []; 
    spine_head_ids = []; 
    spine_volumes_total = []; 
    dist_from_root = [];
    warning('未找到任何棘突 (Region Type = %d)', spine_type);
    return;
end

spine_head_ids = spine_head_ids(:);

head_volumes_terminal = zeros(num_spines, 1);
spine_volumes_total = zeros(num_spines, 1);
dist_from_root = zeros(num_spines, 1);

for i = 1:num_spines
    idx_head = spine_head_ids(i);
    
    dist_from_root(i) = all_dists(idx_head);
    
    head_volumes_terminal(i) = all_vols(idx_head);
    
    current_node = idx_head;
    current_spine_volume = 0;
    
    while true
        current_spine_volume = current_spine_volume + all_vols(current_node);
        
        parent_node = direct_parent(current_node);
        
        if parent_node == 0 || parent_node == current_node || tree.R(parent_node) ~= spine_type
            break;
        end
        
        current_node = parent_node;
    end
    
    spine_volumes_total(i) = current_spine_volume;
end

end