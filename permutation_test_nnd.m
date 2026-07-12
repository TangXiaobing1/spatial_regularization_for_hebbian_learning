function [nnd_distribution, observed_nnd, p_value] = permutation_test_nnd(tree, target_spine_indices, all_spine_indices,spine_type_for_branches, num_permutations)
% Input:
% Output:

if nargin < 5
    num_permutations = 10000;
end

% spine_type_for_branches = tree.R(all_spine_indices(1));
if nargin < 4
    num_permutations = 10000;
    spine_type_for_branches = 4;
end

biological_sect = find_biological_branches(tree, spine_type_for_branches);

observed_nnd = calculate_nnd_for_subset(tree, target_spine_indices, biological_sect);
fprintf('观测到的NND为: %.4f µm\n', observed_nnd);

disp(['正在进行 ', num2str(num_permutations), ' 次置换检验...']);
nnd_distribution = zeros(num_permutations, 1);
M = length(target_spine_indices);

for i = 1:num_permutations
    random_subset_indices = randsample(all_spine_indices, M);
    
    nnd_distribution(i) = calculate_nnd_for_subset(tree, random_subset_indices, biological_sect);
    
    if mod(i, 1000) == 0
        fprintf('已完成 %d / %d 次模拟...\n', i, num_permutations);
    end
end
disp('置换检验完成。');

p_value = sum(nnd_distribution <= observed_nnd) / num_permutations;
fprintf('统计检验p值为: %.4f\n', p_value);

figure;
histogram(nnd_distribution, 100, 'Normalization', 'pdf', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
line([observed_nnd, observed_nnd], ylim, 'Color', 'r', 'LineWidth', 2);
xlabel('Average nearest neighbor distance (NND, µm)');
ylabel('Probability Density');
title({'NND的置换检验结果', '红色竖线 = 真实观测值 | 灰色分布 = 随机模型'});
legend({'Random distribution (null hypothesis)', 'observed'});
grid on;

all_data_for_plot = [nnd_distribution; observed_nnd];
min_val = min(all_data_for_plot);
max_val = max(all_data_for_plot);
data_range = max_val - min_val;
if data_range == 0
    margin = 1;
else
    margin = 0.1 * data_range;
end

xlim([min_val - margin, max_val + margin]);



end
