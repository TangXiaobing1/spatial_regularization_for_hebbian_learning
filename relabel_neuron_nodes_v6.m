function relabeled_tree = relabel_neuron_nodes_v6(original_tree)
% 2 = Dendrite (dendrite)
% 4. Axon Identify：

fprintf('启动节点重标记 (v10 - 轴突近端检测版)...\n');

[relabeled_tree, order_indices] = sort_tree(original_tree, '-LO');

num_nodes = length(relabeled_tree.X);

category_map = struct('unclassified', 0, 'soma', 1, 'dendrite', 2, 'axon', 3, 'spine', 4);
category_names = {'soma', 'dendrite', 'axon', 'spine'};

node_labels = zeros(num_nodes, 1);

fprintf('阶段一：标记 Soma...\n');
soma_radius_threshold = 15; % um
dists_to_root = eucl_tree(relabeled_tree); % [cite: 742]
soma_indices = dists_to_root < soma_radius_threshold;
node_labels(soma_indices) = category_map.soma;

fprintf('阶段二：计算拓扑质量 (Mass)...\n');

seg_len = len_tree(relabeled_tree); % [cite: 754]
idpar = idpar_tree(relabeled_tree); % [cite: 761]

node_mass = seg_len; 
max_child_mass_of_parent = zeros(num_nodes, 1);

% (N -> 2)
for i = num_nodes:-1:2
    pid = idpar(i);
    if pid > 0
        node_mass(pid) = node_mass(pid) + node_mass(i);
        if node_mass(i) > max_child_mass_of_parent(pid)
            max_child_mass_of_parent(pid) = node_mass(i);
        end
    end
end

fprintf('阶段三：初步识别 Spine...\n');

% Parameter setup
spine_abs_max_mass = 10;
neurite_heavy_mass = 15;
sibling_ratio_thr = 0.3;
neck_ratio_thr = 0.7;

for i = 2:num_nodes
    if node_labels(i) == category_map.soma; continue; end
    
    pid = idpar(i);
    parent_label = node_labels(pid);
    
    if parent_label == category_map.spine
        node_labels(i) = category_map.spine;
        continue;
    end
    
    if node_mass(i) > spine_abs_max_mass
        continue;
    end
    
    is_spine = false;
    
    if parent_label == category_map.soma
        is_spine = true;
    elseif node_mass(pid) > neurite_heavy_mass
        is_spine = true;
    else
        max_sibling_mass = max_child_mass_of_parent(pid);
        if node_mass(i) < (max_sibling_mass * sibling_ratio_thr)
            is_spine = true;
        end
        if (relabeled_tree.D(i) / relabeled_tree.D(pid)) < neck_ratio_thr
             is_spine = true;
        end
    end
    
    if is_spine
        node_labels(i) = category_map.spine;
    end
end

fprintf('阶段四：识别 Axon 并修正末端误标...\n');

path_distances = Pvec_tree(relabeled_tree); 

primary_neurite_indices = [];
for i = 2:num_nodes
    if node_labels(idpar(i)) == category_map.soma && ...
       node_labels(i) ~= category_map.spine
        primary_neurite_indices = [primary_neurite_indices, i];
    end
end

axon_start_node = -1;
max_axon_mass = -1;
proximal_check_dist = 100;

fprintf('  - 分析 %d 个一级分支的近端 %.0f um 区域...\n', length(primary_neurite_indices), proximal_check_dist);

for idx = primary_neurite_indices
    subtree_mask = sub_tree(relabeled_tree, idx); % [cite: 1572]
    
    start_dist = path_distances(idx);
    proximal_mask = subtree_mask & (path_distances < (start_dist + proximal_check_dist));
    
    proximal_labels = node_labels(proximal_mask);
    spine_count = sum(proximal_labels == category_map.spine);
    
    if spine_count == 0
        if node_mass(idx) > max_axon_mass
            max_axon_mass = node_mass(idx);
            axon_start_node = idx;
        end
    end
end

if axon_start_node ~= -1
    fprintf('  - 锁定 Axon: 起始点 %d, 长度 %.2f um\n', axon_start_node, max_axon_mass);
    
    is_axon_tree = sub_tree(relabeled_tree, axon_start_node);
    
    update_mask = is_axon_tree & (node_labels ~= category_map.soma);
    node_labels(update_mask) = category_map.axon;
    
    fprintf('  - 已清除 Axon 上误标的 Spine，并统一标记。\n');
else
    fprintf('  - 警告: 未找到符合条件的 Axon (可能所有分支近端都有 Spine)。\n');
end

fprintf('阶段五：剩余节点归类为 Dendrite...\n');
dendrite_mask = (node_labels == 0);
node_labels(dendrite_mask) = category_map.dendrite;

count_soma = sum(node_labels == category_map.soma);
count_dend = sum(node_labels == category_map.dendrite);
count_axon = sum(node_labels == category_map.axon);
count_spine = sum(node_labels == category_map.spine);

fprintf('分类统计:\n Soma: %d\n Dendrite: %d\n Axon: %d\n Spine: %d\n', ...
    count_soma, count_dend, count_axon, count_spine);

relabeled_tree.R = node_labels;
relabeled_tree.rnames = category_names;

end