function [results_struct, elbow_radius] = analyze_spatial_window_similarity_1d(tree, all_spine_ids, all_features, window_sizes, spine_type)
% ANALYZE_SPATIAL_WINDOW_SIMILARITY_1D: Spine similarity within physical-distance windows (1D)
% Input:
% all_spine_ids: [N x 1] spine basenodeID
% Output:

if nargin < 5, spine_type = 4; end
if nargin < 4, window_sizes = [2:0.5:50]; end

%% 1. Topology analysis ---
fprintf('正在解析分支拓扑结构...\n');

biological_sect = find_biological_branches(tree, spine_type);
direct_parent = idpar_tree(tree);
all_path_lengths = Pvec_tree(tree);
parent_paths = ipar_tree(tree);

spine_bases = zeros(size(all_spine_ids));
for i = 1:length(all_spine_ids)
    curr = all_spine_ids(i);
    while tree.R(direct_parent(curr)) == spine_type
        curr = direct_parent(curr);
        if curr == 0 || curr == 1, break; end
    end
    spine_bases(i) = direct_parent(curr);
end

% Assign spines to branches
branch_data = struct('spine_ids', {}, 'dists', {});
valid_branch_count = 0;

for i = 1:size(biological_sect, 1)
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx = find(path_to_root == start_node, 1);
    if isempty(start_idx), continue; end
    branch_nodes = path_to_root(1:start_idx);
    
    mask = ismember(spine_bases, branch_nodes);
    curr_spines = all_spine_ids(mask);
    curr_bases = spine_bases(mask);
    
    if ~isempty(curr_spines)
        dists = all_path_lengths(curr_bases);
        [sorted_dists, sort_idx] = sort(dists);
        valid_branch_count = valid_branch_count + 1;
        branch_data(valid_branch_count).spine_ids = curr_spines(sort_idx);
        branch_data(valid_branch_count).dists = sorted_dists;
    end
end

fprintf('分析完成。共找到 %d 条合并后的生物学分支路径。\n', valid_branch_count);

%% 2. Sliding-window analysis ---
num_windows = length(window_sizes);
results_global = zeros(num_windows, 1);
results_branch = zeros(num_windows, 1);
results_window = zeros(num_windows, 1);

step_size = 1.0;

for w_idx = 1:num_windows
    L = window_sizes(w_idx);
    
    d_obs_all = []; d_glo_all = []; d_bra_all = [];
    
    for b = 1:valid_branch_count
        b_spines = branch_data(b).spine_ids;
        b_dists = branch_data(b).dists;
        if isempty(b_dists), continue; end
        
        min_d = min(b_dists);
        max_d = max(b_dists);
        if (max_d - min_d) < 1, continue; end
        
        starts = min_d : step_size : (max_d - L);
        
        for s = starts
            in_window_idx = find(b_dists >= s & b_dists <= (s + L));
            k = length(in_window_idx);
            
            if k >= 2
                window_spine_ids = b_spines(in_window_idx);
                d_obs_all(end+1) = calc_mean_dist_1d(window_spine_ids, all_spine_ids, all_features);
                
                % Global Random (Density Matched)
                rand_global = randsample(all_spine_ids, k);
                d_glo_all(end+1) = calc_mean_dist_1d(rand_global, all_spine_ids, all_features);
                
                % Branch Random (Density Matched)
                rand_branch = b_spines(randsample(length(b_spines), k));
                d_bra_all(end+1) = calc_mean_dist_1d(rand_branch, all_spine_ids, all_features);
            end
        end
    end
    
    results_global(w_idx) = mean(d_glo_all);
    results_branch(w_idx) = mean(d_bra_all);
    results_window(w_idx) = mean(d_obs_all);
end

results_struct.window_sizes = window_sizes;
results_struct.mean_global = results_global;
results_struct.mean_branch = results_branch;
results_struct.mean_observed = results_window;

[elbow_radius, elbow_y, ~] = find_elbow_point(window_sizes, results_window);
fprintf('自动检测到的有效作用半径 (Elbow Radius): %.2f um\n', elbow_radius);

% Check the current figure visibility
fig_visible = get(0, 'DefaultFigureVisible');
if strcmp(fig_visible, 'on')
    figure('Name', 'Spatial Window Analysis (1D)', 'Color', 'w', 'Position', [100 100 800 600]);
    hold on;
    
    plot(window_sizes, results_global, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'Global Random (Density Matched)');
    plot(window_sizes, results_branch, 'g--^', 'LineWidth', 1.5, 'DisplayName', 'Branch Random (Density Matched)');
    plot(window_sizes, results_window, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Spatial Window (Observed)');
    
    xline(elbow_radius, '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Effective Radius: %.1f \\mum', elbow_radius));
    plot(elbow_radius, elbow_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    
    xlabel('Spatial Window Size (\mum)', 'FontSize', 12);
    ylabel('Average Euclid Pairwise Distance', 'FontSize', 12);
    title({'Morphological Similarity vs. Physical Distance (1D)', ...
           sprintf('Detected Interaction Range: ~%.1f \\mum', elbow_radius)}, 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    hold off;
end

end

function d = calc_mean_dist_1d(ids, all_ids, all_feats)
    [~, rows] = ismember(ids, all_ids);
    feats = all_feats(rows, :);
    if size(feats, 1) < 2
        d = NaN;
    else
        % In 1D, pdist is the absolute difference
        d = mean(pdist(feats, 'euclidean'));
    end
end
