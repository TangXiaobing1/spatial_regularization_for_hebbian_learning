function hotspot_info = find_spine_hotspots(tree, spine_type, window_size, spine_threshold, step_size)
% Dependencies:
% TREES toolbox
% function find_biological_branches
% Input:
% Output:

% --- 1. Default parameters ---
if nargin < 2, spine_type = 4; end
if nargin < 3, window_size = 5; end
if nargin < 4, spine_threshold = 20; end
if nargin < 5, step_size = 1; end

biological_sect = find_biological_branches(tree, spine_type);
num_branches = size(biological_sect, 1);
direct_parent_indices = idpar_tree(tree);
all_path_lengths_from_root = Pvec_tree(tree, len_tree(tree));
all_segment_lengths = len_tree(tree);
parent_paths_matrix = ipar_tree(tree);
spine_indices = tree.R == spine_type;
spine_attachment_points = direct_parent_indices(spine_indices);
backbone_nodes = [];

hotspot_info = struct('branch_id', {}, 'hotspot_start', {}, 'hotspot_end', {}, ...
                      'spine_count', {}, 'density', {}, 'backbone_nodes_in_hotspot', {});

for i = 1:num_branches
    start_node = biological_sect(i, 1);
    end_node = biological_sect(i, 2);
    
    path_to_root = parent_paths_matrix(end_node, :);
    path_to_root = path_to_root(path_to_root > 0);
    start_idx_in_path = find(path_to_root == start_node);
    if isempty(start_idx_in_path), continue; end
    branch_nodes = fliplr(path_to_root(1:start_idx_in_path));
    
    if length(branch_nodes) < 2, branch_length = 0; else
        branch_length = sum(all_segment_lengths(branch_nodes(2:end)));
    end
    
    if branch_length < window_size, continue; end
    
    spines_on_this_branch_mask = ismember(spine_attachment_points, branch_nodes);
    attachment_points_on_branch = spine_attachment_points(spines_on_this_branch_mask);
    if isempty(attachment_points_on_branch), continue; end
    
    path_start_node_metric_length = all_path_lengths_from_root(start_node);
    spine_positions_on_branch = all_path_lengths_from_root(attachment_points_on_branch) - path_start_node_metric_length;
    
    backbone_nodes_path_lengths_on_branch = all_path_lengths_from_root(branch_nodes) - path_start_node_metric_length;
    
    window_starts = 0:step_size:(branch_length - window_size);
    
    for w_start = window_starts
        w_end = w_start + window_size;
        spine_count_in_window = sum(spine_positions_on_branch >= w_start & spine_positions_on_branch < w_end);
        
        if spine_count_in_window >= spine_threshold
            nodes_in_hotspot_mask = (backbone_nodes_path_lengths_on_branch >= w_start & backbone_nodes_path_lengths_on_branch <= w_end);
            hotspot_backbone_nodes = branch_nodes(nodes_in_hotspot_mask);

            hotspot_data.branch_id = i;
            hotspot_data.hotspot_start = w_start;
            hotspot_data.hotspot_end = w_end;
            hotspot_data.spine_count = spine_count_in_window;
            hotspot_data.density = spine_count_in_window / window_size;
            hotspot_data.backbone_nodes_in_hotspot = hotspot_backbone_nodes'; % 转置为行向量方便查看
            
            % new
            backbone_nodes = [backbone_nodes,hotspot_backbone_nodes];
            hotspot_info(end+1) = hotspot_data;
        end
    end
end
disp(['分析完成。共找到 ', num2str(length(hotspot_info)), ' 个树突棘热点区域。']);

figure;
hold on; 
plot_tree(tree, [0.8 0.8 0.8]);
pointer_tree(tree, backbone_nodes, [], [0.6627, 0.2000, 0.1490]);
title('Hotspot Region on Tree');
hold off;


end