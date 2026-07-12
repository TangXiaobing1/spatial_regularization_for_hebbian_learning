function [profile_data] = analyze_mexican_hat_profile_v2(tree, varargin)
% Input:

    % --- 0. Parse parameters ---
    p = inputParser;
    addRequired(p, 'tree');
    addParameter(p, 'LargeThreshold', 0.2);
    addParameter(p, 'DiameterThreshold', 2);
    parse(p, tree, varargin{:});
    
    large_threshold_ratio = p.Results.LargeThreshold;
    dia_thresh = p.Results.DiameterThreshold;

    % --- 1. Identify spine structure ---
    fprintf('正在识别 Spine Head 和结构...\n');
    [spine_group_ids, head_node_indices] = identify_spine_heads_by_diameter(tree, 'DiameterThreshold', dia_thresh);
    
    node_vols = vol_tree(tree);
    
    % --- 2. Aggregate data ---
    unique_spine_ids = unique(spine_group_ids);
    unique_spine_ids(unique_spine_ids == 0) = []; 
    
    num_spines = length(unique_spine_ids);
    spine_vols = zeros(num_spines, 1);
    spine_anchors = zeros(num_spines, 1);
    dendrite_anchors = zeros(num_spines, 1);
    spine_anchor_locs = zeros(num_spines, 1);
    
    parents = idpar_tree(tree);
    path_lengths_to_root = Pvec_tree(tree); 
    
    fprintf('正在聚合 %d 个 Spine 的体积和位置...\n', num_spines);
    
    for k = 1:num_spines
        s_id = unique_spine_ids(k);
        
        % 2.1 Computevolume
        is_head_in_spine = ismember(head_node_indices, find(spine_group_ids == s_id));
        current_head_indices = head_node_indices(is_head_in_spine);
        spine_vols(k) = sum(node_vols(current_head_indices));
        
        % 2.2 Find Spine Anchor
        all_nodes_in_spine = find(spine_group_ids == s_id);
        anchor_node = [];
        
        for i = 1:length(all_nodes_in_spine)
            curr = all_nodes_in_spine(i);
            p = parents(curr);
            if p > 0 && spine_group_ids(p) ~= s_id
                anchor_node = curr;
                break;
            end
        end
        
        if isempty(anchor_node)
            anchor_node = all_nodes_in_spine(1); 
        end
        spine_anchors(k) = anchor_node;
        
        dendrite_parent = parents(anchor_node);
        
        dendrite_anchors(k) = dendrite_parent;
        
        if dendrite_parent > 0 && dendrite_parent <= length(path_lengths_to_root)
            spine_anchor_locs(k) = path_lengths_to_root(dendrite_parent);
        else
            spine_anchor_locs(k) = path_lengths_to_root(anchor_node);
        end
    end
    
    % --- 3. Data preprocessing ---
    log_vols = log10(spine_vols + eps);
    z_vols = (log_vols - mean(log_vols)) / std(log_vols);
    
    sorted_v = sort(z_vols, 'descend');
    cut_idx = ceil(length(z_vols) * large_threshold_ratio);
    threshold_val = sorted_v(cut_idx);
    is_large_global = (z_vols >= threshold_val);
    
    
    biological_sect = find_biological_branches(tree); 
    parent_paths = ipar_tree(tree);
    
    all_dists = [];
    all_neigh_vols = [];
    all_shuffled_vols = [];
    
    fprintf('正在进行分支邻居分析...\n');
    fprintf('共找到 %d 条拓扑分支。\n', size(biological_sect, 1));
    
    for b = 1:size(biological_sect, 1)
        path_nodes = parent_paths(biological_sect(b, 2), :);
        branch_nodes = path_nodes(path_nodes > 0);
        
        mask = ismember(dendrite_anchors, branch_nodes);
        
        if sum(mask) < 2, continue; end 
        
        local_pos = spine_anchor_locs(mask);
        local_vols = z_vols(mask);
        local_is_large = is_large_global(mask);
        
        local_vols_shuffled = local_vols(randperm(length(local_vols)));
        
        large_indices = find(local_is_large);
        
        for i = 1:length(large_indices)
            center_idx = large_indices(i);
            center_pos = local_pos(center_idx);
            
            dists = abs(local_pos - center_pos);
            valid_neighbors = (dists > 0);
            
            d_vec = dists(valid_neighbors);
            v_vec = local_vols(valid_neighbors);
            s_vec = local_vols_shuffled(valid_neighbors);
            
            all_dists = [all_dists; d_vec];
            all_neigh_vols = [all_neigh_vols; v_vec];
            all_shuffled_vols = [all_shuffled_vols; s_vec];
        end
    end
    
    profile_data.distances = all_dists;
    profile_data.neighbor_vols = all_neigh_vols;
    profile_data.shuffled_vols = all_shuffled_vols;
    
    fprintf('分析完成。共收集 %d 组邻居关系。\n', length(all_dists));
end