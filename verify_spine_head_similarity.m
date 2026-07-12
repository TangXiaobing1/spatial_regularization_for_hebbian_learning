function verify_spine_head_similarity(tree, spine_type, max_cluster_size)
% Input:

if nargin < 3, max_cluster_size = 30; end
if nargin < 2, spine_type = 4; end

fprintf('正在提取 Head Volume 并进行对数变换...\n');

is_terminal = T_tree(tree);
is_spine = (tree.R == spine_type);
spine_head_ids = find(is_terminal & is_spine);

if isempty(spine_head_ids)
    error('未找到任何棘突 (Region Type = %d)', spine_type);
end

all_vols = vol_tree(tree);
raw_head_volumes = all_vols(spine_head_ids);

log_head_volumes = log10(raw_head_volumes + eps);

norm_log_volumes = zscore(log_head_volumes);

fprintf('共分析 %d 个棘突。\n', length(spine_head_ids));

fprintf('正在解析分支拓扑结构 (含回溯逻辑)...\n');

biological_sect = find_biological_branches(tree, spine_type);
direct_parent = idpar_tree(tree);
all_path_lengths = Pvec_tree(tree);
parent_paths = ipar_tree(tree);

spine_bases = zeros(size(spine_head_ids));
for i = 1:length(spine_head_ids)
    curr = spine_head_ids(i);
    while tree.R(direct_parent(curr)) == spine_type
        curr = direct_parent(curr);
        if curr == 0 || curr == 1, break; end % Prevent an infinite loop
    end
    spine_bases(i) = direct_parent(curr);
end

% 3. Assign spines to branches
all_attachment_points = spine_bases;
spines_on_branches = {}; 
valid_branch_count = 0;

for i = 1:size(biological_sect, 1)
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx = find(path_to_root == start_node, 1);
    if isempty(start_idx), continue; end
    branch_nodes = path_to_root(1:start_idx);
    
    mask = ismember(all_attachment_points, branch_nodes);
    current_branch_spines = spine_head_ids(mask);
    current_attachments = all_attachment_points(mask);
    
    if ~isempty(current_branch_spines)
        dists = all_path_lengths(current_attachments);
        [~, sort_idx] = sort(dists);
        spines_on_branches{valid_branch_count+1} = current_branch_spines(sort_idx);
        valid_branch_count = valid_branch_count + 1;
    end
end
cluster_sizes = 2:max_cluster_size;
mean_global = zeros(length(cluster_sizes), 1);
mean_branch = zeros(length(cluster_sizes), 1);
mean_adjacent = zeros(length(cluster_sizes), 1);

num_permutations = 1000;

fprintf('开始三模型对比分析 (Max Size = %d)...\n', max_cluster_size);
h_wait = waitbar(0, '正在进行置换检验...');

for k = 1:length(cluster_sizes)
    c_size = cluster_sizes(k);
    waitbar(k/length(cluster_sizes), h_wait, sprintf('分析 Cluster Size = %d...', c_size));
    
    valid_branches_idx = find(cellfun(@length, spines_on_branches) >= c_size);
    
    if isempty(valid_branches_idx)
        mean_global(k:end) = NaN; mean_branch(k:end) = NaN; mean_adjacent(k:end) = NaN;
        break; % No sufficiently long branches remain，Stop
    end
    
    d_global = zeros(num_permutations, 1);
    d_branch = zeros(num_permutations, 1);
    d_adj = zeros(num_permutations, 1);
    
    for p = 1:num_permutations
        % 1. Global Random
        g_ids = randsample(spine_head_ids, c_size);
        d_global(p) = calc_dist(g_ids, spine_head_ids, norm_log_volumes);
        
        rand_ptr = randi(length(valid_branches_idx));
        b_idx = valid_branches_idx(rand_ptr);
        target_branch_spines = spines_on_branches{b_idx};
        n_spines = length(target_branch_spines);
        
        b_rand_ids = target_branch_spines(randsample(n_spines, c_size));
        d_branch(p) = calc_dist(b_rand_ids, spine_head_ids, norm_log_volumes);
        
        start_pos = randi(n_spines - c_size + 1);
        adj_ids = target_branch_spines(start_pos : start_pos + c_size - 1);
        d_adj(p) = calc_dist(adj_ids, spine_head_ids, norm_log_volumes);
    end
    
    mean_global(k) = mean(d_global);
    mean_branch(k) = mean(d_branch);
    mean_adjacent(k) = mean(d_adj);
end
close(h_wait);

%% 4. Visualization ---
figure('Name', 'Head Volume Similarity Analysis', 'Color', 'w', 'Position', [100 100 900 600]);
hold on;

plot(cluster_sizes, mean_global, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Global Random');
plot(cluster_sizes, mean_branch, 'g--^', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Branch Random');
plot(cluster_sizes, mean_adjacent, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'DisplayName', 'Spatially Adjacent');

title({'synapse weight similarity within cluster', '(log(Head_Volume)as proxy )'}, 'FontSize', 14);
xlabel('Cluster Size (N spines)', 'FontSize', 12);
ylabel('Average Pairwise Distance (Normalized Log-Vol)', 'FontSize', 12);
legend('Location', 'best');
grid on;

msg = sprintf('如果 蓝线 < 绿线，说明相邻的 Spine 在大小上确实更相似。');
text(min(cluster_sizes), max([mean_global; mean_branch])*1.05, msg, 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');

hold off;

end

function d = calc_dist(ids, all_ids, values)
    [~, idx] = ismember(ids, all_ids);
    vals = values(idx);
    dist_matrix = pdist(vals, 'euclidean');
    d = mean(dist_matrix);
end