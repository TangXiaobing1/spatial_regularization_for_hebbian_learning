function [total_spine_count, nodes_per_spine, spine_base_nodes] = analyze_spine_node_num(relabeled_tree)
% Input:
% Output:


if ~isfield(relabeled_tree, 'R') || ~isfield(relabeled_tree, 'rnames')
    error('输入的树结构体不完整，缺少 .R 或 .rnames 字段。');
end

rnames_list = relabeled_tree.rnames;

if ismember('soma', rnames_list)
    soma_id = find(strcmp(rnames_list, 'soma'));
else
    fprintf('警告: 在 .rnames 中未找到 "soma"。将假定 soma_id = 1。\n');
    soma_id = 1;
end

if ismember('neurite', rnames_list)
    neurite_id = find(strcmp(rnames_list, 'neurite'));
else
    fprintf('警告: 在 .rnames 中未找到 "neurite"。将假定 neurite_id = 3。\n');
    neurite_id = 3;
end

if ismember('spine', rnames_list)
    spine_id = find(strcmp(rnames_list, 'spine'));
else
    fprintf('警告: 在 .rnames 中未找到 "spine"。将假定 spine_id = 4。\n');
    spine_id = 4;
end

if isempty(spine_id)
    fprintf('树中没有找到 "spine" 类别。\n');
    total_spine_count = 0;
    nodes_per_spine = [];
    spine_base_nodes = [];
    return;
end

all_spine_indices = find(relabeled_tree.R == spine_id);

if isempty(all_spine_indices)
    fprintf('没有节点被标记为棘突。\n');
    total_spine_count = 0;
    nodes_per_spine = [];
    spine_base_nodes = [];
    return;
end

parent_indices = idpar_tree(relabeled_tree);
child_indices_map = create_child_map(relabeled_tree);


parents_of_spines = parent_indices(all_spine_indices);
parent_labels = relabeled_tree.R(parents_of_spines);

is_base_node_mask = (parent_labels ~= spine_id);
spine_base_nodes = all_spine_indices(is_base_node_mask);

total_spine_count = length(spine_base_nodes);
fprintf('共找到 %d 个独立的spine。\n', total_spine_count);

nodes_per_spine = zeros(total_spine_count, 1);

for i = 1:total_spine_count
    base_node = spine_base_nodes(i);
    
    count = 0;
    queue = base_node;
    visited = [];
    
    while ~isempty(queue)
        current_node = queue(1);
        queue(1) = [];
        
        if ismember(current_node, visited)
            continue;
        end
        
        visited = [visited; current_node];
        count = count + 1;
        
        children = child_indices_map{current_node};
        
        for j = 1:length(children)
            child_node = children(j);
            if relabeled_tree.R(child_node) == spine_id
                queue = [queue; child_node];
            end
        end
    end
    
    nodes_per_spine(i) = count;
end

fprintf('spine大小分析完成。\n');

end


function child_map = create_child_map(tree)
    num_nodes = length(tree.X);
    child_map = cell(num_nodes, 1);
    parent_indices = idpar_tree(tree);
    for i = 1:num_nodes
        children = find(parent_indices == i);
        child_map{i} = children;
    end
end