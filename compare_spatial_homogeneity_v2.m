function [results] = compare_spatial_homogeneity_v2(tree, all_spine_ids, all_normalized_features, cluster_size, num_permutations)
% Input:
% Output:

if nargin < 5, num_permutations = 5000; end
if nargin < 4, error('必须提供 cluster_size。'); end

fprintf('正在进行拓扑分析 (Cluster Size = %d)...\n', cluster_size);

spine_type_for_branches = tree.R(all_spine_ids(1));
biological_sect = find_biological_branches(tree, spine_type_for_branches);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
parent_paths_matrix = ipar_tree(tree);
all_attachment_points = direct_parent_indices(all_spine_ids);

spines_sorted_by_branch = cell(size(biological_sect, 1), 1);
for i = 1:size(biological_sect, 1)
    path_to_root = parent_paths_matrix(biological_sect(i, 2), :); 
    path_to_root = path_to_root(path_to_root > 0);
    start_idx = find(path_to_root == biological_sect(i, 1));
    if isempty(start_idx), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx));
    
    mask = ismember(all_attachment_points, branch_nodes);
    spine_ids_on_branch = all_spine_ids(mask);
    attachment_pts = all_attachment_points(mask);
    
    if ~isempty(spine_ids_on_branch)
        dists = all_path_lengths_from_root(attachment_pts);
        [~, sort_idx] = sort(dists);
        spines_sorted_by_branch{i} = spine_ids_on_branch(sort_idx);
    end
end

valid_branch_indices = find(cellfun(@length, spines_sorted_by_branch) >= cluster_size);
if isempty(valid_branch_indices)
    error('没有分支包含足够的棘突 (需 >= %d)。', cluster_size);
end
fprintf('分析完成。共有 %d 个有效分支可用于抽样。\n', length(valid_branch_indices));

disp(['开始执行 ', num2str(num_permutations), ' 次置换模拟...']);

null_dist_global   = zeros(num_permutations, 1);
null_dist_adjacent = zeros(num_permutations, 1);
null_dist_branch   = zeros(num_permutations, 1);

for i = 1:num_permutations
    global_ids = randsample(all_spine_ids, cluster_size);
    null_dist_global(i) = get_dist(global_ids, all_spine_ids, all_normalized_features);
    
    random_idx_ptr = randi(length(valid_branch_indices)); 
    chosen_branch_idx = valid_branch_indices(random_idx_ptr);
    spines_on_this_branch = spines_sorted_by_branch{chosen_branch_idx};
    num_spines_on_branch = length(spines_on_this_branch);
    
    start_pos = randi(num_spines_on_branch - cluster_size + 1);
    adjacent_ids = spines_on_this_branch(start_pos : start_pos + cluster_size - 1);
    null_dist_adjacent(i) = get_dist(adjacent_ids, all_spine_ids, all_normalized_features);
    
    branch_random_indices = randsample(num_spines_on_branch, cluster_size);
    branch_random_ids = spines_on_this_branch(branch_random_indices);
    null_dist_branch(i) = get_dist(branch_random_ids, all_spine_ids, all_normalized_features);
    
    if mod(i, 1000) == 0, fprintf('进度: %d / %d\n', i, num_permutations); end
end

results.dist_global   = null_dist_global;
results.dist_adjacent = null_dist_adjacent;
results.dist_branch   = null_dist_branch;

results.mean_global   = mean(null_dist_global);
results.mean_adjacent = mean(null_dist_adjacent);
results.mean_branch   = mean(null_dist_branch);

% T-Tests
[~, results.p_global_vs_branch] = ttest2(null_dist_global, null_dist_branch);
[~, results.p_branch_vs_adjacent] = ttest2(null_dist_branch, null_dist_adjacent);

fprintf('\n--- 结果摘要 (Cluster Size = %d) ---\n', cluster_size);
fprintf('1. 全局随机均值 (Global):   %.4f\n', results.mean_global);
fprintf('2. 分支内随机均值 (Branch): %.4f (分支背景效应)\n', results.mean_branch);
fprintf('3. 空间邻近均值 (Adjacent): %.4f (局部相互作用)\n', results.mean_adjacent);
fprintf('\n');
fprintf('P-Value [Global vs Branch]:   %.4e (是否受分支约束?)\n', results.p_global_vs_branch);
fprintf('P-Value [Branch vs Adjacent]: %.4e (是否存在邻近效应?)\n', results.p_branch_vs_adjacent);

figure('Name', '三模型形态同质性对比', 'Color', 'w', 'Position', [100 100 900 600]);
hold on;

histogram(null_dist_global, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Global Random');
histogram(null_dist_branch, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.7 0.3], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Within-Branch Random');
histogram(null_dist_adjacent, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.4 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Spatially Adjacent');

yl = ylim;
line([results.mean_global results.mean_global], yl, 'Color', [0.8 0 0], 'LineWidth', 2, 'LineStyle', '--');
line([results.mean_branch results.mean_branch], yl, 'Color', [0 0.6 0], 'LineWidth', 2, 'LineStyle', '--');
line([results.mean_adjacent results.mean_adjacent], yl, 'Color', [0 0 0.8], 'LineWidth', 2, 'LineStyle', '--');

xlabel('Average Pairwise Euclidean Distance');
ylabel('Probability Density');
title({['Cluster Size = ' num2str(cluster_size)], ...
       ['Global(' num2str(results.mean_global, '%.2f') ') > Branch(' num2str(results.mean_branch, '%.2f') ') > Adj(' num2str(results.mean_adjacent, '%.2f') ')']});
legend('Location', 'best');
grid on; hold off;

end

function avg_dist = get_dist(ids, all_ids, all_feats)
    [~, rows] = ismember(ids, all_ids);
    feats = all_feats(rows, :);
    dists = pdist(feats, 'euclidean');
    avg_dist = mean(dists);
end