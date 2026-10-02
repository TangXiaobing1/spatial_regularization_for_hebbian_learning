function result = analyze_pair_level_weight_class_permutations( ...
    tree, thresholds, n_shuffles, seed)
% Pair-level 2-14 um counts for five quintiles plus Top10.
% Returns branch-preserving and within-neuron global permutation samples.

if nargin < 3, n_shuffles = 1000; end
if nargin < 4, seed = 20260804; end
thresholds = thresholds(:)'; % P20, P40, P60, P80, P90
n_classes = 6;
result = struct('status', 'not_processed', 'agglomeration_id', '', ...
    'n_spines', 0, 'observed_pair_count', 0, ...
    'observed_class_count', zeros(1, n_classes), ...
    'branch_shuffle_class_count', zeros(n_shuffles, n_classes), ...
    'global_shuffle_class_count', zeros(n_shuffles, n_classes));
result.agglomeration_id = extract_agglomeration_id(tree);

[spine_group_ids, head_node_indices] = identify_spine_heads_by_diameter( ...
    tree, 'DiameterThreshold', 3);
unique_spine_ids = unique(spine_group_ids);
unique_spine_ids(unique_spine_ids == 0) = [];
n_spines = numel(unique_spine_ids);
result.n_spines = n_spines;
if n_spines < 2
    result.status = 'fewer_than_two_spines';
    return;
end

node_volumes = vol_tree(tree);
parents = idpar_tree(tree);
path_lengths_to_root = Pvec_tree(tree);
spine_volumes = zeros(n_spines, 1);
dendrite_anchors = zeros(n_spines, 1);
spine_locations = zeros(n_spines, 1);
for spine_index = 1:n_spines
    spine_id = unique_spine_ids(spine_index);
    nodes = find(spine_group_ids == spine_id);
    head_indices = head_node_indices(ismember(head_node_indices, nodes));
    spine_volumes(spine_index) = sum(node_volumes(head_indices));
    anchor_node = [];
    for node_index = 1:numel(nodes)
        node = nodes(node_index);
        parent = parents(node);
        if parent > 0 && spine_group_ids(parent) ~= spine_id
            anchor_node = node;
            break;
        end
    end
    if isempty(anchor_node), anchor_node = nodes(1); end
    dendrite_parent = parents(anchor_node);
    dendrite_anchors(spine_index) = dendrite_parent;
    if dendrite_parent > 0 && dendrite_parent <= numel(path_lengths_to_root)
        spine_locations(spine_index) = path_lengths_to_root(dendrite_parent);
    else
        spine_locations(spine_index) = path_lengths_to_root(anchor_node);
    end
end

log_volumes = log10(spine_volumes + eps);
log_sd = std(log_volumes);
if ~isfinite(log_sd) || log_sd == 0
    result.status = 'invalid_volume_variance';
    return;
end
z_volumes = (log_volumes - mean(log_volumes)) ./ log_sd;
sorted_volumes = sort(z_volumes, 'descend');
center_threshold = sorted_volumes(max(1, ceil(n_spines * 0.20)));
is_center = z_volumes >= center_threshold;

evalc('biological_sections = find_biological_branches(tree);');
parent_paths = ipar_tree(tree);
global_target_weights = zeros(n_spines, 1);
batch_size = 100;
rng(seed, 'twister');

for branch = 1:size(biological_sections, 1)
    path_nodes = parent_paths(biological_sections(branch, 2), :);
    branch_nodes = path_nodes(path_nodes > 0);
    branch_mask = ismember(dendrite_anchors, branch_nodes);
    if sum(branch_mask) < 2, continue; end
    branch_spine_indices = find(branch_mask);
    local_locations = spine_locations(branch_mask);
    local_volumes = z_volumes(branch_mask);
    center_indices = find(is_center(branch_mask));
    if isempty(center_indices), continue; end
    [source_matrix, target_matrix] = ndgrid(center_indices, 1:numel(local_volumes));
    distances = abs(local_locations(source_matrix) - local_locations(target_matrix));
    target_indices = target_matrix(distances >= 2 & distances < 14);
    if isempty(target_indices), continue; end

    local_weights = accumarray(target_indices(:), 1, ...
        [numel(local_volumes), 1]);
    result.observed_pair_count = result.observed_pair_count + sum(local_weights);
    global_target_weights(branch_spine_indices) = ...
        global_target_weights(branch_spine_indices) + local_weights;
    for class_index = 1:n_classes
        result.observed_class_count(class_index) = ...
            result.observed_class_count(class_index) + ...
            sum(local_weights(class_mask(local_volumes, thresholds, class_index)));
    end

    n_local = numel(local_volumes);
    for first_shuffle = 1:batch_size:n_shuffles
        last_shuffle = min(first_shuffle + batch_size - 1, n_shuffles);
        rows = first_shuffle:last_shuffle;
        [~, permutations] = sort(rand(numel(rows), n_local), 2);
        shuffled_values = local_volumes(permutations);
        for class_index = 1:n_classes
            flags = class_mask(shuffled_values, thresholds, class_index);
            result.branch_shuffle_class_count(rows, class_index) = ...
                result.branch_shuffle_class_count(rows, class_index) + ...
                flags * local_weights;
        end
    end
end

% A separate deterministic stream makes global shuffles reproducible and
% independent of how many biological paths were traversed above.
rng(seed + 500000000, 'twister');
for first_shuffle = 1:batch_size:n_shuffles
    last_shuffle = min(first_shuffle + batch_size - 1, n_shuffles);
    rows = first_shuffle:last_shuffle;
    [~, permutations] = sort(rand(numel(rows), n_spines), 2);
    shuffled_values = z_volumes(permutations);
    for class_index = 1:n_classes
        flags = class_mask(shuffled_values, thresholds, class_index);
        result.global_shuffle_class_count(rows, class_index) = ...
            flags * global_target_weights;
    end
end
result.status = 'ok';
end


function flags = class_mask(values, thresholds, class_index)
switch class_index
    case 1
        flags = values <= thresholds(1);
    case 2
        flags = values > thresholds(1) & values <= thresholds(2);
    case 3
        flags = values > thresholds(2) & values <= thresholds(3);
    case 4
        flags = values > thresholds(3) & values <= thresholds(4);
    case 5
        flags = values > thresholds(4);
    case 6
        flags = values > thresholds(5);
end
end


function id = extract_agglomeration_id(tree)
id = '';
if isfield(tree, 'name') && ~isempty(tree.name)
    token = regexp(char(tree.name), '^\d+', 'match', 'once');
    if ~isempty(token), id = token; end
end
end
