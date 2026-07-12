function [total_path_length, total_spines_on_path, overall_spine_density, bins_info]=calculate_spine_density_on_path(tree, path_start_node, path_end_node, bin_length)
% Input:
% Output:

parent_paths_matrix = ipar_tree(tree); % [cite: 236]
path_to_root_for_end_node = parent_paths_matrix(path_end_node, :);
path_to_root_for_end_node = path_to_root_for_end_node(path_to_root_for_end_node > 0);

if path_start_node ~= 1
    start_node_index_in_path = find(path_to_root_for_end_node == path_start_node);
    if isempty(start_node_index_in_path)
        error('指定的 path_start_node 不在 path_end_node 到根节点的路径上。');
    end
    actual_path_nodes = fliplr(path_to_root_for_end_node(1:start_node_index_in_path));
else
    actual_path_nodes = fliplr(path_to_root_for_end_node);
end

if isempty(actual_path_nodes) || length(actual_path_nodes) < 2
    disp('定义的路径无效或过短。');
    total_path_length = 0;
    total_spines_on_path = 0;
    overall_spine_density = 0;
    bins_info = struct('start_distance', {}, 'end_distance', {}, 'spine_count', {}, 'density', {});
    return;
end

disp(['路径上的主干节点数量: ', num2str(length(actual_path_nodes))]);
% disp(actual_path_nodes);

% --- 2. Computepath length ---
segment_lengths = len_tree(tree); % [cite: 341]
path_segment_indices = actual_path_nodes(2:end);
total_path_length = sum(segment_lengths(path_segment_indices)); % [cite: 341]

if total_path_length == 0
    disp('路径总长度为0，无法计算密度。');
    total_spines_on_path = 0;
    overall_spine_density = 0;
    bins_info = struct('start_distance', {}, 'end_distance', {}, 'spine_count', {}, 'density', {});
    return;
end
disp(['路径总长度: ', num2str(total_path_length), ' µm']);

direct_parent_indices = idpar_tree(tree); % [cite: 233]
spine_type = 4;
all_spine_indices = find(tree.R == spine_type);

spines_on_path_indices = [];
for i = 1:length(actual_path_nodes)
    current_path_node = actual_path_nodes(i);
    attached_spines = find(direct_parent_indices == current_path_node & tree.R == spine_type);
    spines_on_path_indices = [spines_on_path_indices; attached_spines];
end
spines_on_path_indices = unique(spines_on_path_indices);
total_spines_on_path = length(spines_on_path_indices);
disp(['路径上的总树突棘数量: ', num2str(total_spines_on_path)]);

if total_path_length > 0
    overall_spine_density = total_spines_on_path / total_path_length; % [cite: 150]
else
    overall_spine_density = 0;
end
disp(['路径上的总体树突棘密度: ', num2str(overall_spine_density), ' spines/µm']);


all_nodes_metric_path_length_to_root = Pvec_tree(tree, len_tree(tree)); % [cite: 245]

path_start_node_metric_length_to_root = 0;
if path_start_node ~= 1
   path_start_node_metric_length_to_root = all_nodes_metric_path_length_to_root(path_start_node);
end

spine_attachment_points_on_path = direct_parent_indices(spines_on_path_indices);

spine_distances_along_path = all_nodes_metric_path_length_to_root(spine_attachment_points_on_path) - path_start_node_metric_length_to_root;
spine_distances_along_path = max(0, spine_distances_along_path);
spine_distances_along_path = min(total_path_length, spine_distances_along_path);

num_bins = floor(total_path_length / bin_length);
if num_bins == 0 && total_path_length > 0
    num_bins = 1;
    actual_bin_edges = [0, total_path_length];
else
    actual_bin_edges = 0:bin_length:(num_bins*bin_length);
    if actual_bin_edges(end) < total_path_length
        actual_bin_edges = [actual_bin_edges, total_path_length];
        num_bins = num_bins + 1;
    end
end


if total_spines_on_path > 0 && ~isempty(actual_bin_edges) && length(actual_bin_edges) > 1
    spine_counts_in_bins = histcounts(spine_distances_along_path, actual_bin_edges);
else
    spine_counts_in_bins = zeros(1, num_bins);
end

bins_info = struct('start_distance', {}, 'end_distance', {}, 'spine_count', {}, 'density', {});
disp('--- 分箱统计 ---');
for i = 1:num_bins
    bins_info(i).start_distance = actual_bin_edges(i);
    bins_info(i).end_distance = actual_bin_edges(i+1);
    current_bin_actual_length = bins_info(i).end_distance - bins_info(i).start_distance;
    bins_info(i).spine_count = spine_counts_in_bins(i);
    if current_bin_actual_length > 0
        bins_info(i).density = bins_info(i).spine_count / current_bin_actual_length;
    else
        bins_info(i).density = 0;
    end
    disp(['Bin ', num2str(i), ': ', ...
          num2str(bins_info(i).start_distance), ' - ', num2str(bins_info(i).end_distance), ' µm, ', ...
          'Spine Count: ', num2str(bins_info(i).spine_count), ', ', ...
          'Density: ', num2str(bins_info(i).density), ' spines/µm']);
end

figure;
hold on; 
plot_tree(tree, [0.8 0.8 0.8]);
pointer_tree(tree, actual_path_nodes, [], [1 0 0]);

hold off;

end
