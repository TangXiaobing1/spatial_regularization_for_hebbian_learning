function [normalized_features, raw_features_matrix, spine_base_ids_in_order, feature_names, feature_mean, feature_std] = extract_and_normalize_spine_features(tree, spine_type)

% head/neck ratio
% Input:
% Output:


if nargin < 2
    spine_type = 4;
end

is_terminal_vec = T_tree(tree);
is_spine_vec = (tree.R == spine_type);
spine_head_indices = find(is_terminal_vec & is_spine_vec);

num_spines = length(spine_head_indices);
if num_spines == 0
    [normalized_features, raw_features_matrix, spine_base_ids_in_order, feature_names, feature_mean, feature_std] = deal([]);
    feature_names = {};
    disp('在树中未找到任何作为末梢点的棘突。');
    return;
end

spine_base_ids_in_order = zeros(num_spines, 1);

direct_parent_indices = idpar_tree(tree);
all_lengths_vec = len_tree(tree);
all_surfaces_vec = surf_tree(tree);
all_volumes_vec = vol_tree(tree);
raw_features_matrix = zeros(num_spines, 6);
feature_names = {'Length', 'HeadDiameter', 'SurfaceArea', 'Volume', 'NeckDiameter', 'HeadNeckRatio'};

fprintf('正在为 %d 个独立棘突提取形态学特征...\n', num_spines);

for i = 1:num_spines
    current_spine_head = spine_head_indices(i);
    current_node = current_spine_head;
    
    current_spine_path_nodes = [];
    attachment_point_on_backbone = -1;
    while true
        if tree.R(current_node) ~= spine_type
            attachment_point_on_backbone = current_node;
            break;
        end
        current_spine_path_nodes = [current_node, current_spine_path_nodes];
        parent_node = direct_parent_indices(current_node);
        if parent_node == current_node || parent_node == 0, break; end
        current_node = parent_node;
    end
    
    if ~isempty(current_spine_path_nodes)
        spine_base_ids_in_order(i) = current_spine_path_nodes(1);
    else
        spine_base_ids_in_order(i) = current_spine_head;
    end
    
    path_with_attachment = [current_spine_path_nodes, attachment_point_on_backbone];
    
    feature_length = sum(all_lengths_vec(current_spine_path_nodes));
    feature_head_diameter = tree.D(current_spine_head);
    feature_surface_area = sum(all_surfaces_vec(current_spine_path_nodes));
    feature_volume = sum(all_volumes_vec(current_spine_path_nodes));
    feature_neck_diameter = min(tree.D(path_with_attachment));
    feature_head_neck_ratio = feature_head_diameter / (feature_neck_diameter + eps);
    
    raw_features_matrix(i, :) = [
        feature_length, feature_head_diameter, feature_surface_area, ...
        feature_volume, feature_neck_diameter, feature_head_neck_ratio
    ];
end

[normalized_features, feature_mean, feature_std] = zscore(raw_features_matrix);
fprintf('特征提取与归一化完成。\n');
end