function [results] = test_cluster_homogeneity(tree, all_spine_ids, all_normalized_features, hotspot_info, num_permutations)
% Input:
% Output:

if nargin < 5, num_permutations = 5000; end

observed_similarities = zeros(length(hotspot_info), 1);
for i = 1:length(hotspot_info)
    cluster_ids = hotspot_info(i).spine_ids_in_hotspot;
    [~, feature_rows] = ismember(cluster_ids, all_spine_ids);
    cluster_features = all_normalized_features(feature_rows, :);
    observed_similarities(i) = calculate_intra_cluster_similarity(cluster_features);
end

results.observed_individual_similarities = observed_similarities; 
results.observed_mean_homogeneity = mean(observed_similarities);
fprintf('真实空间簇的平均形态距离为: %.4f\n', results.observed_mean_homogeneity);

disp('正在进行一次性拓扑分析...');
spine_type_for_branches = tree.R(all_spine_ids(1));
biological_sect = find_biological_branches(tree, spine_type_for_branches);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
parent_paths_matrix = ipar_tree(tree);
all_attachment_points = direct_parent_indices(all_spine_ids);
cluster_sizes = [hotspot_info.spine_count];

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


disp(['正在进行 ', num2str(num_permutations), ' 次置换检验...']);
null_dist_global = zeros(num_permutations, 1);
null_dist_adjacent = zeros(num_permutations, 1);

for i = 1:num_permutations
    similarities_global_random = [];
    similarities_adjacent_random = [];
    
    for j = 1:length(cluster_sizes)
        current_cluster_size = cluster_sizes(j);
        
        global_random_ids = randsample(all_spine_ids, current_cluster_size);
        [~, feature_rows] = ismember(global_random_ids, all_spine_ids);
        sim_global = calculate_intra_cluster_similarity(all_normalized_features(feature_rows, :));
        similarities_global_random = [similarities_global_random, sim_global];
        
        valid_branches = find(cellfun(@length, spines_sorted_by_branch) >= current_cluster_size);
        if isempty(valid_branches), continue; end
        random_branch_idx = randsample(valid_branches, 1);
        
        spines_on_chosen_branch = spines_sorted_by_branch{random_branch_idx};
        start_pos = randi(length(spines_on_chosen_branch) - current_cluster_size + 1);
        
        adjacent_random_ids = spines_on_chosen_branch(start_pos : start_pos + current_cluster_size - 1);
        [~, feature_rows] = ismember(adjacent_random_ids, all_spine_ids);
        sim_adjacent = calculate_intra_cluster_similarity(all_normalized_features(feature_rows, :));
        similarities_adjacent_random = [similarities_adjacent_random, sim_adjacent];
    end
    
    null_dist_global(i) = mean(similarities_global_random);
    if ~isempty(similarities_adjacent_random)
        null_dist_adjacent(i) = mean(similarities_adjacent_random);
    else
        null_dist_adjacent(i) = nan;
    end
    
    if mod(i, 500) == 0, fprintf('已完成 %d / %d 次模拟...\n', i, num_permutations); end
end
disp('检验完成。');

results.null_dist_global_random = null_dist_global;
results.null_dist_adjacent_random = rmmissing(null_dist_adjacent);
results.p_value_global_random = sum(results.null_dist_global_random <= results.observed_mean_homogeneity) / num_permutations;
results.p_value_adjacent_random = sum(results.null_dist_adjacent_random <= results.observed_mean_homogeneity) / length(results.null_dist_adjacent_random);

disp('正在生成四张分析图表...');

figure('Name', '检验1: vs. 全局随机');
histogram(results.null_dist_global_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
h_obs1 = line([results.observed_mean_homogeneity, results.observed_mean_homogeneity], ylim, 'Color', 'r', 'LineWidth', 2.5);
title({'检验1: 真实热点 vs. 全局随机簇', sprintf('p = %.4f', results.p_value_global_random)});
xlabel('平均簇内形态距离');
ylabel('概率密度');
legend(h_obs1, {'真实热点平均值'}); % Create a legend entry only for the red line
grid on;

figure('Name', '检验2: vs. 空间邻近簇');
histogram(results.null_dist_adjacent_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on;
h_obs2 = line([results.observed_mean_homogeneity, results.observed_mean_homogeneity], ylim, 'Color', 'r', 'LineWidth', 2.5);
title({'检验2: 真实热点 vs. 空间邻近簇', sprintf('p = %.4f', results.p_value_adjacent_random)});
xlabel('平均簇内形态距离');
ylabel('概率密度');
legend(h_obs2, {'hotspot region'});
grid on;

figure('Name', '两种零假设分布与真实平均值对比');
hold on;
h1 = histogram(results.null_dist_global_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.7);
h2 = histogram(results.null_dist_adjacent_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.7);
h_obs3 = line([results.observed_mean_homogeneity, results.observed_mean_homogeneity], ylim, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '--');
title('两种随机模型的分布与真实平均值对比');
xlabel('平均簇内形态距离');
ylabel('概率密度');
legend([h1, h2, h_obs3], {'global random spines', 'adjacent spines', 'hotspot region'});
grid on;

figure('Name', '两种零假设分布与每个真实热点值对比');
hold on;
h4 = histogram(results.null_dist_global_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6);
h5 = histogram(results.null_dist_adjacent_random, 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.6);
num_real_clusters = length(results.observed_individual_similarities);
h_scatter = plot(results.observed_individual_similarities, zeros(num_real_clusters, 1), 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
title('两种随机模型的分布与每个真实热点值对比');
xlabel('簇内形态距离');
ylabel('概率密度');
legend([h4, h5, h_scatter], {'全局随机簇', '空间邻近簇', '每一个真实热点的值'});
grid on;

end