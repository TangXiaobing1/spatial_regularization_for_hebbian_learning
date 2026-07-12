function [profile_data] = analyze_mexican_hat_profile(tree, spine_ids, head_volumes, large_threshold_ratio)
% Input:

    if nargin < 4, large_threshold_ratio = 0.2; end

    % --- 1. Data preprocessing (Log + Z-score) ---
    log_vols = log10(head_volumes + eps);
    z_vols = (log_vols - mean(log_vols)) / std(log_vols);
    
    % Definelarge-spine threshold
    sorted_v = sort(z_vols, 'descend');
    cut_idx = ceil(length(z_vols) * large_threshold_ratio);
    threshold_val = sorted_v(cut_idx);
    is_large_global = (z_vols >= threshold_val);

    % --- 2. Topology mapping (Group by branch) ---
    direct_parent = idpar_tree(tree);
    all_path_lengths = Pvec_tree(tree);
    spine_bases = zeros(size(spine_ids));
    
    % Trace back to the attachment point
    for i = 1:length(spine_ids)
        curr = spine_ids(i);
        while tree.R(direct_parent(curr)) == 4 
            curr = direct_parent(curr);
            if curr <= 1, break; end 
        end
        spine_bases(i) = direct_parent(curr);
    end
    
    biological_sect = find_biological_branches(tree, 4);
    parent_paths = ipar_tree(tree);
    
    all_dists = [];
    all_neigh_vols = [];
    all_shuffled_vols = [];
    
    for b = 1:size(biological_sect, 1)
        % Getspines on the branch
        path_nodes = parent_paths(biological_sect(b, 2), :);
        branch_nodes = path_nodes(path_nodes > 0);
        mask = ismember(spine_bases, branch_nodes);
        
        if sum(mask) < 3, continue; end
        
        % Getdata for the current branch
        local_pos = all_path_lengths(spine_bases(mask));
        local_vols = z_vols(mask);
        local_is_large = is_large_global(mask);
        
        local_vols_shuffled = local_vols(randperm(length(local_vols)));
        
        large_indices = find(local_is_large);
        
        for i = 1:length(large_indices)
            center_idx = large_indices(i);
            center_pos = local_pos(center_idx);
            
            % Computedistances to all other spines
            dists = abs(local_pos - center_pos);
            
            valid_neighbors = (dists > 0);
            
            % 1. distance
            d_vec = dists(valid_neighbors);
            v_vec = local_vols(valid_neighbors);
            % 3. shuffled volumes
            s_vec = local_vols_shuffled(valid_neighbors);
            
            % Add to the output pool
            all_dists = [all_dists; d_vec];
            all_neigh_vols = [all_neigh_vols; v_vec];
            all_shuffled_vols = [all_shuffled_vols; s_vec];
        end
    end
    
    profile_data.distances = all_dists;
    profile_data.neighbor_vols = all_neigh_vols;
    profile_data.shuffled_vols = all_shuffled_vols;
end