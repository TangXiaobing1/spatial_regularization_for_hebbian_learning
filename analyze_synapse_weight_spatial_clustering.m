function [results] = analyze_synapse_weight_spatial_clustering(tree, spine_ids, head_volumes, scan_windows, do_plot)
% ANALYZE_ADAPTIVE_SPATIAL_CLUSTERING_V3:

    if nargin < 5, do_plot = false; end
    if nargin < 4 || isempty(scan_windows), scan_windows = [2:1:30, 32:2:50]; end

    % --- 1. Data preprocessing & Topology analysis ---
    log_vols = zscore(log10(head_volumes + eps)); 
    step_size = 1.0; 
    
    direct_parent = idpar_tree(tree);
    all_path_lengths = Pvec_tree(tree);
    spine_bases = zeros(size(spine_ids));
    
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
    branches = struct('dists', {}, 'vols', {});
    valid_branch_count = 0;
    
    for b = 1:size(biological_sect, 1)
        path_nodes = parent_paths(biological_sect(b, 2), :);
        branch_nodes = path_nodes(path_nodes > 0);
        mask = ismember(spine_bases, branch_nodes);
        if sum(mask) < 2, continue; end
        
        b_dists = all_path_lengths(spine_bases(mask));
        [sorted_dists, sort_idx] = sort(b_dists);
        current_branch_vols = log_vols(mask);
        sorted_vols = current_branch_vols(sort_idx);
        
        valid_branch_count = valid_branch_count + 1;
        branches(valid_branch_count).dists = sorted_dists;
        branches(valid_branch_count).vols = sorted_vols;
    end
    
    num_wins = length(scan_windows);
    mean_obs = zeros(num_wins, 1);
    mean_branch = zeros(num_wins, 1);
    mean_global = zeros(num_wins, 1); % new Global
    
    cached_data = cell(num_wins, 1);
    
    for w_idx = 1:num_wins
        win_size = scan_windows(w_idx);
        d_obs_list = [];
        d_bra_list = [];
        d_glo_list = [];
        
        for b = 1:valid_branch_count
            b_dists = branches(b).dists;
            b_vols = branches(b).vols;
            
            if (max(b_dists) - min(b_dists)) < win_size, continue; end
            starts = min(b_dists) : step_size : (max(b_dists) - win_size);
            
            for s = starts
                in_win = (b_dists >= s & b_dists <= (s + win_size));
                k = sum(in_win);
                
                if k >= 2
                    win_vols = b_vols(in_win);
                    
                    % 1. Local (Observed)
                    d_obs_list(end+1) = mean(pdist(win_vols, 'euclidean'));
                    
                    % 2. Branch Random (Density Matched)
                    rand_bra = randsample(length(b_vols), k);
                    d_bra_list(end+1) = mean(pdist(b_vols(rand_bra), 'euclidean'));
                    
                    % 3. Global Random (Density Matched) - new
                    rand_glo = randsample(length(log_vols), k);
                    d_glo_list(end+1) = mean(pdist(log_vols(rand_glo), 'euclidean'));
                end
            end
        end
        
        if isempty(d_obs_list)
            mean_obs(w_idx) = NaN; mean_branch(w_idx) = NaN; mean_global(w_idx) = NaN;
        else
            mean_obs(w_idx) = mean(d_obs_list);
            mean_branch(w_idx) = mean(d_bra_list);
            mean_global(w_idx) = mean(d_glo_list);
        end
        
        cached_data{w_idx}.obs = d_obs_list;
        cached_data{w_idx}.bra = d_bra_list;
    end
    
    delta = mean_branch - mean_obs;
    
    max_delta = max(delta);
    
    if max_delta <= 0
        detect_idx = 1;
        effective_radius = NaN;
    else
        threshold = max_delta * 0.10; 
        candidates = find(delta < threshold);
        
        if isempty(candidates)
            detect_idx = num_wins;
        else
            detect_idx = candidates(1);
        end
        effective_radius = scan_windows(detect_idx);
    end
    
    if isnan(effective_radius)
        test_idx = [];
    else
        [~, max_diff_idx] = max(delta(1:detect_idx));
        test_idx = max_diff_idx;
    end

    if isempty(test_idx)
        results.ElbowRadius = NaN;
        results.IsSignificant = false;
        results.PValue = NaN;
        results.EffectSize = NaN;
    else
        target = cached_data{test_idx};
        obs = target.obs;
        bra = target.bra;
        
        if length(obs) > 10
            [h, p] = ttest(obs, bra, 'Tail', 'left');
            diffs = bra - obs;
            es = mean(diffs) / std(diffs);
            
            results.ElbowRadius = effective_radius;
            results.PeakEffectRadius = scan_windows(test_idx);
            results.IsSignificant = (h == 1);
            results.PValue = p;
            results.EffectSize = es;
        else
            results.ElbowRadius = NaN; results.IsSignificant = false; results.PValue = NaN; results.EffectSize = NaN;
        end
    end
    
    if do_plot
        figure('Color', 'w', 'Name', 'Spatial Clustering V3'); hold on;
        
        plot(scan_windows, mean_global, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Global Random');
        plot(scan_windows, mean_branch, 'g--^', 'LineWidth', 1.5, 'DisplayName', 'Branch Random');
        plot(scan_windows, mean_obs, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Spatial Window');
        
        if ~isnan(effective_radius)
            xline(effective_radius, '--k', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Convergence: %.1f \\mum', effective_radius));
            
            plot(results.PeakEffectRadius, mean_obs(test_idx), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
                'DisplayName', sprintf('Peak Effect (d=%.2f)', results.EffectSize));
        end
        
        title({sprintf('Weight Clustering Analysis'), ...
               sprintf('Sig: %d | P=%.1e | Effect Size=%.2f', results.IsSignificant, results.PValue, results.EffectSize)}, ...
               'FontSize', 12);
        xlabel('Window Size (\mum)');
        ylabel('Avg Pairwise Distance (Z-score Log Vol)');
        legend('Location', 'best');
        grid on; hold off;
    end
end