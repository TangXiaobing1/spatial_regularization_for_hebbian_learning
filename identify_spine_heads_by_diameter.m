function [spine_ids, head_indices] = identify_spine_heads_by_diameter(tree, varargin)
% Input:
% Output:

    % --- 0. Parse parameters ---
    p = inputParser;
    addRequired(p, 'tree');
    addParameter(p, 'DiameterThreshold', 3);
    addParameter(p, 'ContinuousSteps', 3);
    parse(p, tree, varargin{:});
    
    thresh_multiplier = p.Results.DiameterThreshold;
    continuous_steps_req = p.Results.ContinuousSteps;
    
    N = length(tree.X);
    
    % InitializeOutput
    spine_ids = zeros(N, 1);
    head_indices = [];
    
    idpar = idpar_tree(tree);
    is_terminal = T_tree(tree);
    
    child_indices = cell(N, 1);
    for i = 1:N
        p = idpar(i);
        if p > 0 && p <= N
            child_indices{p} = [child_indices{p}, i];
        end
    end
    
    anchor_nodes = [];
    for i = 1:N
        if tree.R(i) == 4
            p = idpar(i);
            if p == 1 || (p > 0 && tree.R(p) ~= 4)
                anchor_nodes = [anchor_nodes, i];
            end
        end
    end
    
    current_spine_index = 0;
    
    for s = 1:length(anchor_nodes)
        root_spine_node = anchor_nodes(s);
        
        spine_nodes = [];
        queue = root_spine_node;
        
        while ~isempty(queue)
            curr = queue(1);
            queue(1) = [];
            spine_nodes = [spine_nodes, curr];
            
            childs = child_indices{curr};
            for k = 1:length(childs)
                if tree.R(childs(k)) == 4
                    queue = [queue, childs(k)];
                end
            end
        end
        
        if isempty(spine_nodes), continue; end
        
        [~, min_idx_local] = min(tree.D(spine_nodes));
        neck_center = spine_nodes(min_idx_local);
        
        valid_nodes = [neck_center];
        
        p = idpar(neck_center);
        if p > 0 && tree.R(p) == 4
            valid_nodes = [valid_nodes, p];
        end
        
        childs = child_indices{neck_center};
        if ~isempty(childs)
            valid_nodes = [valid_nodes, childs(1)];
        end
        
        neck_diameter = mean(tree.D(valid_nodes));
        
        curr = neck_center;
        found_head_start = false;
        head_start_node = [];
        max_steps = 100; 
        step = 0;
        
        continuous_queue = []; 
        
        while step < max_steps
            step = step + 1;
            childs = child_indices{curr};
            
            if isempty(childs)
                break;
            end
            
            next_node = childs(1);
            
            if tree.D(next_node) > thresh_multiplier * neck_diameter
                continuous_queue = [continuous_queue, next_node];
                
                if length(continuous_queue) >= continuous_steps_req
                    head_start_node = continuous_queue(1);
                    found_head_start = true;
                    break;
                end
            else
                continuous_queue = [];
            end
            
            curr = next_node;
        end
        
        current_spine_index = current_spine_index + 1;
        spine_ids(spine_nodes) = current_spine_index;
        
        % Collect Head nodeindex
        if found_head_start
            queue = head_start_node;
            while ~isempty(queue)
                curr = queue(1);
                queue(1) = [];
                
                head_indices = [head_indices, curr];
                
                childs = child_indices{curr};
                for k = 1:length(childs)
                    if tree.R(childs(k)) == 4
                        queue = [queue, childs(k)];
                    end
                end
            end
        else
            is_spine_term = ismember(spine_nodes, find(is_terminal));
            term_indices = spine_nodes(is_spine_term);
            head_indices = [head_indices, term_indices];
        end
    end
    
    if ~isempty(head_indices)
        head_indices = head_indices(:);
    end
    
end