function [U, S, V, volume_matrix_normalized, bin_edges] = analyze_spine_volume_svd(tree, spine_type, bin_length, min_branch_length)
% Input:
% Output:

%% Parameter setup
if nargin < 2, spine_type = 4; end
if nargin < 3, bin_length = 2; end
if nargin < 4, min_branch_length = 5; end

fprintf('=== 步骤0: 获取Spine数据 ===\n');

[head_volumes_terminal, spine_head_ids, ~, dist_from_root] = ...
    calculate_spine_head_volumes(tree, spine_type);

biological_sect = find_biological_branches(tree, spine_type, min_branch_length);

num_branches = size(biological_sect, 1);
fprintf('共 %d 条生物学分支，%d 个spine heads\n', num_branches, length(spine_head_ids));

fprintf('=== 步骤1: 绝对对齐 ===\n');

parent_paths_matrix = ipar_tree(tree);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));

all_attachment_points = direct_parent_indices(spine_head_ids);

spines_per_branch = cell(num_branches, 1);
local_distances = cell(num_branches, 1);

for i = 1:num_branches
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    spines_on_branch_mask = ismember(all_attachment_points, branch_nodes);
    spine_ids_on_branch = spine_head_ids(spines_on_branch_mask);
    
    if ~isempty(spine_ids_on_branch)
        attachment_points_on_branch = all_attachment_points(spines_on_branch_mask);
        path_lengths_on_branch = all_path_lengths_from_root(attachment_points_on_branch);
        
        start_dist = all_path_lengths_from_root(start_node);
        local_dist = path_lengths_on_branch - start_dist;
        
        % Sort
        [local_dist, sort_idx] = sort(local_dist);
        spines_per_branch{i} = spine_ids_on_branch(sort_idx);
        local_distances{i} = local_dist;
    end
end

fprintf('=== 步骤2: 空间求和 (Bin长度=%.1f μm) ===\n', bin_length);

max_local_dist = 0;
for i = 1:num_branches
    if ~isempty(local_distances{i})
        curr_max = max(local_distances{i});
        if curr_max > max_local_dist
            max_local_dist = curr_max;
        end
    end
end

% Create bin edges
num_bins = ceil(max_local_dist / bin_length);
bin_edges = 0:bin_length:num_bins * bin_length;

% Initializematrix（rows = branches, columns = bins）
volume_matrix = zeros(num_branches, num_bins);

for i = 1:num_branches
    if isempty(spines_per_branch{i}), continue; end
    
    branch_spines = spines_per_branch{i};
    branch_dists = local_distances{i};
    
    for j = 1:length(branch_spines)
        spine_idx = find(spine_head_ids == branch_spines(j));
        if isempty(spine_idx), continue; end
        volume = head_volumes_terminal(spine_idx);
        
        bin_idx = floor(branch_dists(j) / bin_length) + 1;
        if bin_idx > num_bins, bin_idx = num_bins; end
        
        volume_matrix(i, bin_idx) = volume_matrix(i, bin_idx) + volume;
    end
end

fprintf('体积矩阵大小: %d x %d (分支 x Bins)\n', size(volume_matrix, 1), size(volume_matrix, 2));

fprintf('=== 步骤3: 消除贫富差距 (行归一化) ===\n');

row_sums = sum(volume_matrix, 2);
volume_matrix_normalized = zeros(size(volume_matrix));

for i = 1:num_branches
    if row_sums(i) > 0
        volume_matrix_normalized(i, :) = volume_matrix(i, :) / row_sums(i);
    end
end

fprintf('归一化完成。有效分支数: %d\n', sum(row_sums > 0));

fprintf('=== 步骤4: 终极审判 (SVD) ===\n');

[U, S, V] = svd(volume_matrix_normalized, 'econ');

fprintf('SVD分解完成。\n');
fprintf('  - U 大小: %d x %d (分支特征)\n', size(U, 1), size(U, 2));
fprintf('  - S 大小: %d x %d (奇异值对角矩阵)\n', size(S, 1), size(S, 2));
fprintf('  - V 大小: %d x %d (空间Bin特征)\n', size(V, 1), size(V, 2));

singular_values = diag(S);
total_variance = sum(singular_values.^2);
explained_variance = (singular_values.^2 / total_variance) * 100;
fprintf('\n前5个奇异值解释方差比:\n');
for i = 1:min(5, length(singular_values))
    fprintf('  SV%d: %.2f%%\n', i, explained_variance(i));
end

fprintf('\n分析完成！\n');
end
