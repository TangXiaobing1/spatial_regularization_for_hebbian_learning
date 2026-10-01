function results = within_branch_permutation_test_v2(tree, target_percentile, neighbor_percentile, n_permutations, distance_range, P_obs_external)
% Input:
% distance_range: [min, max] distancerange，default[2, 14]

    %% Default parameters
    if nargin < 2 || isempty(target_percentile)
        target_percentile = [80, 100];  % defaultTop 20%
    end
    if nargin < 3 || isempty(neighbor_percentile)
        neighbor_percentile = [90, 100];  % defaultTop 10%
    end
    if nargin < 4 || isempty(n_permutations)
        n_permutations = 1000;
    end
    if nargin < 5 || isempty(distance_range)
        distance_range = [2, 14];  % default2-14 μm
    end
    if nargin < 6
        P_obs_external = [];
    end
    
    [spine_ids, head_indices] = identify_spine_heads_by_diameter(tree);
    
    if isempty(head_indices) || isempty(spine_ids)
        error('未能识别到spine head');
    end
    
    node_volumes = vol_tree(tree);
    head_node_volumes = node_volumes(head_indices);
    
    head_spine_ids = spine_ids(head_indices);
    num_spines_all = max(spine_ids);
    
    spine_volumes = accumarray(...
        head_spine_ids, ...
        head_node_volumes, ...
        [num_spines_all, 1], ...
        @sum, ...
        0 ...
    );
    
    valid_spine_idx = spine_volumes > 0;
    spine_volumes = spine_volumes(valid_spine_idx);
    n_spines = length(spine_volumes);
    
    if n_spines < 10
        error('有效spine数量不足（<10）');
    end
    
    % 1.2 Getspineposition
    parents = idpar_tree(tree);
    path_lengths_to_root = Pvec_tree(tree);
    
    spine_positions = zeros(n_spines, 1);
    valid_spine_ids = find(valid_spine_idx);
    
    for s = 1:n_spines
        spine_id = valid_spine_ids(s);
        all_nodes_in_spine = find(spine_ids == spine_id);
        
        anchor_node = [];
        for i = 1:length(all_nodes_in_spine)
            curr = all_nodes_in_spine(i);
            p = parents(curr);
            if p > 0 && spine_ids(p) ~= spine_id
                anchor_node = curr;
                break;
            end
        end
        
        if isempty(anchor_node)
            anchor_node = all_nodes_in_spine(1);
        end
        
        dendrite_parent = parents(anchor_node);
        if dendrite_parent > 0 && dendrite_parent <= length(path_lengths_to_root)
            spine_positions(s) = path_lengths_to_root(dendrite_parent);
        else
            spine_positions(s) = path_lengths_to_root(anchor_node);
        end
    end
    
    target_threshold = prctile(spine_volumes, target_percentile(1));
    neighbor_threshold = prctile(spine_volumes, neighbor_percentile(1));
    
    is_target = spine_volumes >= target_threshold;
    is_neighbor_target = spine_volumes >= neighbor_threshold;
    
    n_target = sum(is_target);
    
    if n_target < 2
        error('目标突触数量不足（<2）');
    end
    
    % 1.4 ComputeP_obs
    if ~isempty(P_obs_external)
        P_obs = P_obs_external;
        P_obs_computed = compute_proportion_in_range_v2(spine_positions, is_target, is_neighbor_target, distance_range);
    else
        P_obs = compute_proportion_in_range_v2(spine_positions, is_target, is_neighbor_target, distance_range);
        P_obs_computed = P_obs;
    end
    
    P_null_distribution = zeros(n_permutations, 1);
    
    biological_sect = find_biological_branches(tree);
    parent_paths = ipar_tree(tree);
    
    for perm = 1:n_permutations
        volumes_shuffled = shuffle_volumes_within_branches_v2(...
            spine_volumes, spine_positions, biological_sect, parent_paths, valid_spine_ids);
        
        target_threshold_shuffled = prctile(volumes_shuffled, target_percentile(1));
        neighbor_threshold_shuffled = prctile(volumes_shuffled, neighbor_percentile(1));
        
        is_target_shuffled = volumes_shuffled >= target_threshold_shuffled;
        is_neighbor_shuffled = volumes_shuffled >= neighbor_threshold_shuffled;
        
        P_null_distribution(perm) = compute_proportion_in_range_v2(...
            spine_positions, is_target_shuffled, is_neighbor_shuffled, distance_range);
    end
    
    eps_val = 1e-10;
    p_value = (sum(P_null_distribution <= P_obs) + eps_val) / (n_permutations + eps_val);
    
    results = struct();
    results.P_obs = P_obs;
    results.P_obs_computed = P_obs_computed;
    results.P_null_distribution = P_null_distribution;
    results.p_value = p_value;
    results.target_percentile = target_percentile;
    results.neighbor_percentile = neighbor_percentile;
    results.distance_range = distance_range;
    results.target_threshold = target_threshold;
    results.neighbor_threshold = neighbor_threshold;
    results.n_spines = n_spines;
    results.n_target = n_target;
    results.n_permutations = n_permutations;
end

function P = compute_proportion_in_range_v2(positions, is_target, is_neighbor_target, distance_range)
    min_dist = distance_range(1);
    max_dist = distance_range(2);
    
    target_indices = find(is_target);
    n_target = length(target_indices);
    
    proportions = zeros(n_target, 1);
    
    for i = 1:n_target
        center_idx = target_indices(i);
        center_pos = positions(center_idx);
        
        % Computedistance
        dists = abs(positions - center_pos);
        
        neighbor_idx = (dists >= min_dist) & (dists <= max_dist) & (dists > 0);
        
        if sum(neighbor_idx) >= 1
            proportions(i) = mean(is_neighbor_target(neighbor_idx));
        else
            proportions(i) = NaN;
        end
    end
    
    if sum(~isnan(proportions)) >= 1
        P = nanmean(proportions);
    else
        P = NaN;
    end
end

function volumes_shuffled = shuffle_volumes_within_branches_v2(...
    spine_volumes, spine_positions, biological_sect, parent_paths, valid_spine_ids)
    
    n_spines = length(spine_volumes);
    volumes_shuffled = zeros(n_spines, 1);
    shuffled_flags = false(n_spines, 1);
    
    for b = 1:size(biological_sect, 1)
        path_nodes = parent_paths(biological_sect(b, 2), :);
        branch_nodes = path_nodes(path_nodes > 0);
        
        mask = false(n_spines, 1);
        for s = 1:n_spines
            spine_node = valid_spine_ids(s);
            if ismember(spine_node, branch_nodes)
                mask(s) = true;
            end
        end
        
        if sum(mask) >= 2
            local_vols = spine_volumes(mask);
            local_vols_shuffled = local_vols(randperm(length(local_vols)));
            volumes_shuffled(mask) = local_vols_shuffled;
            shuffled_flags(mask) = true;
        end
    end
    
    unshuffled = ~shuffled_flags;
    if sum(unshuffled) > 0
        volumes_shuffled(unshuffled) = spine_volumes(unshuffled);
    end
end
