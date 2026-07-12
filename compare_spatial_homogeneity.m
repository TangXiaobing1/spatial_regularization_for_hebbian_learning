function [results] = compare_spatial_homogeneity(tree, all_spine_ids, all_normalized_features, cluster_size, num_permutations)
% Input:
% Output:

if nargin < 5, num_permutations = 5000; end
if nargin < 4, error('必须提供要测试的 cluster_size。'); end

disp('正在进行一次性拓扑分析...');
spine_type_for_branches = tree.R(all_spine_ids(1));
biological_sect = find_biological_branches(tree, spine_type_for_branches);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
parent_paths_matrix = ipar_tree(tree);
all_attachment_points = direct_parent_indices(all_spine_ids);

spines_sorted_by_branch = cell(size(biological_sect, 1), 1);
for i = 1:size(biological_sect, 1)
    start_node = biological_sect(i, 1); end_node = biological_sect(i, 2);
    path_to_root = parent_paths_matrix(end_node, :); path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node); if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    spines_on_branch_mask = ismember(all_attachment_points, branch_nodes);
    spine_ids_on_branch = all_spine_ids(spines_on_branch_mask);
    attachment_points_on_branch = all_attachment_points(spines_on_branch_mask);
    
    if ~isempty(spine_ids_on_branch)
        path_lengths_on_branch = all_path_lengths_from_root(attachment_points_on_branch);
        [~, sort_idx] = sort(path_lengths_on_branch);
        spines_sorted_by_branch{i} = spine_ids_on_branch(sort_idx);
    end
end
disp('拓扑分析与棘排序完成。');

disp(['正在进行 ', num2str(num_permutations), ' 次置换检验，测试簇大小为 ', num2str(cluster_size), '...']);
null_dist_global = zeros(num_permutations, 1);
null_dist_adjacent = zeros(num_permutations, 1);

valid_branches = find(cellfun(@length, spines_sorted_by_branch) >= cluster_size);
if isempty(valid_branches)
    error('没有任何分支包含足够数量的棘来创建大小为 %d 的邻近簇。', cluster_size);
end

for i = 1:num_permutations
    global_random_ids = randsample(all_spine_ids, cluster_size);
    [~, feature_rows] = ismember(global_random_ids, all_spine_ids);
    null_dist_global(i) = calculate_intra_cluster_similarity(all_normalized_features(feature_rows, :));
    
    random_branch_idx = randsample(valid_branches, 1);
    spines_on_chosen_branch = spines_sorted_by_branch{random_branch_idx};
    start_pos = randi(length(spines_on_chosen_branch) - cluster_size + 1);
    adjacent_random_ids = spines_on_chosen_branch(start_pos : start_pos + cluster_size - 1);
    [~, feature_rows] = ismember(adjacent_random_ids, all_spine_ids);
    null_dist_adjacent(i) = calculate_intra_cluster_similarity(all_normalized_features(feature_rows, :));
    
    if mod(i, 500) == 0, fprintf('已完成 %d / %d 次模拟...\n', i, num_permutations); end
end
disp('检验完成。');

results.null_dist_global_random = null_dist_global;
results.null_dist_adjacent_random = null_dist_adjacent;
results.mean_global = mean(null_dist_global);
results.mean_adjacent = mean(null_dist_adjacent);

[~, p_ttest] = ttest2(null_dist_global, null_dist_adjacent);
results.p_value_ttest = p_ttest;

fprintf('\n--- 结果摘要 ---\n');
fprintf('全局随机簇的平均形态距离: %.4f\n', results.mean_global);
fprintf('空间邻近簇的平均形态距离: %.4f\n', results.mean_adjacent);
fprintf('比较两个分布的双样本t检验 p-value: %.4f\n', results.p_value_ttest);


figure('Name', '两种零假设模型的形态距离分布对比');
hold on;
h1 = histogram(results.null_dist_global_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.7);
h2 = histogram(results.null_dist_adjacent_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.996 0.311 0.893], 'FaceAlpha', 0.7);

line([results.mean_global, results.mean_global], ylim, 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'LineStyle', '--');
line([results.mean_adjacent, results.mean_adjacent], ylim, 'Color', [0.878, 0.275, 0.730], 'LineWidth', 2, 'LineStyle', '--');

title({'零假设模型形态距离分布对比', sprintf('簇大小 = %d, t-test p = %.4f', cluster_size, results.p_value_ttest)});
xlabel('Average pairwise Euclidean distance in clusters');
ylabel('Probability density');
legend([h1, h2], {'Global Random', 'Spatially Adjacent'}, 'Location', 'northeast');
grid on;
hold off;

end