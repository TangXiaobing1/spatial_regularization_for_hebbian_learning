function relabeled_tree = relabel_neuron_nodes(original_tree)
% Input:
% Output:

%% 1. Initialize
fprintf('启动节点重标记过程 (v5 - 最终修复版)...\n');
relabeled_tree = original_tree;
num_nodes = length(relabeled_tree.X);

category_map = struct('unclassified', 0, 'soma', 1, 'neurite', 2, 'spine', 3);
category_names = {'soma', 'neurite', 'spine'};

node_labels = zeros(num_nodes, 1);
relabeled_tree.R = node_labels;

fprintf('标记Soma节点...\n');
soma_radius_threshold = 15; % um
soma_indices = find(eucl_tree(relabeled_tree, 1) < soma_radius_threshold);
node_labels(soma_indices) = category_map.soma;
fprintf('找到 %d 个节点被标记为 "soma"。\n', length(soma_indices));

fprintf('阶段一：执行启发式初始棘突标记...\n');

all_terminal_indices_h = find(T_tree(relabeled_tree));
segment_lengths_h = len_tree(relabeled_tree);
parent_indices_h = idpar_tree(relabeled_tree);
spine_length_threshold_h = 20; % um

for i = 1:length(all_terminal_indices_h)
    node_idx = all_terminal_indices_h(i);
    if node_labels(node_idx) == category_map.soma
        continue;
    end
    if segment_lengths_h(node_idx) < spine_length_threshold_h
        node_labels(node_idx) = category_map.spine;
        parent_idx = parent_indices_h(node_idx);
        if node_labels(parent_idx) ~= category_map.soma
            node_labels(parent_idx) = category_map.spine;
        end
    end
end
fprintf('初步标记了 %d 个节点为 "spine"。\n', sum(node_labels == category_map.spine));

fprintf('阶段二：开始迭代优化与回溯标记...\n');
relabeled_tree.R = node_labels;

segment_lengths = len_tree(relabeled_tree);
parent_indices = idpar_tree(relabeled_tree);
is_branch_point_vec = B_tree(relabeled_tree); 
max_total_spine_length = 10;
max_spine_neck_diameter = 2;
num_iterations = 2;
for iter = 1:num_iterations
    fprintf('  - 第%d轮优化...\n', iter);
    
    true_terminal_vec = find_non_spine_termination_points(relabeled_tree, category_map.spine);
    all_terminals_vec = T_tree(relabeled_tree);
    spine_head_candidates = find(all_terminals_vec & ~true_terminal_vec);
    
    node_labels(node_labels == category_map.spine) = category_map.unclassified;
    
    for i = 1:length(spine_head_candidates)
        current_node_idx = spine_head_candidates(i);
        path_walked_length = 0;
        
        while true
            if node_labels(current_node_idx) == category_map.soma, break; end
            
            parent_node_idx = parent_indices(current_node_idx);
            if parent_node_idx == current_node_idx || ...
               is_branch_point_vec(parent_node_idx) || ...
               relabeled_tree.D(current_node_idx) > max_spine_neck_diameter || ...
               path_walked_length > max_total_spine_length
                node_labels(current_node_idx) = category_map.spine;
                break;
            end
            
            node_labels(current_node_idx) = category_map.spine;
            path_walked_length = path_walked_length + segment_lengths(current_node_idx);
            current_node_idx = parent_node_idx;
        end
    end
    
    relabeled_tree.R = node_labels;
end
fprintf('迭代完成，最终有 %d 个节点被标记为 "spine"。\n', sum(node_labels == category_map.spine));

fprintf('标记Neurite节点...\n');
neurite_indices = find(node_labels == category_map.unclassified);
node_labels(neurite_indices) = category_map.neurite;
fprintf('找到 %d 个节点被标记为 "neurite"。\n', length(neurite_indices));

relabeled_tree.R = node_labels;
relabeled_tree.rnames = category_names;
fprintf('节点重标记完成。\n');
end