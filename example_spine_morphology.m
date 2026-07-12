Hp = load('H01_pyramidal.mat');

num_neurons = length(Hp.trees_rearranged_cell);
db_normalized_features = cell(num_neurons, 1);
db_spine_base_ids      = cell(num_neurons, 1);
db_valid_neuron_flags  = false(num_neurons, 1);

tic;
for i = 1:num_neurons
    tree = Hp.trees_rearranged_cell{i};
    if isempty(tree) || ~isfield(tree, 'R')
        continue;
    end
    
    [norm_features, ~, base_ids, ~, ~, ~] = extract_and_normalize_spine_features(tree, 4);
    
    if isempty(norm_features) || size(norm_features, 1) < 2
        continue;
    end
    
    db_normalized_features{i} = norm_features;
    db_spine_base_ids{i}      = base_ids;
    db_valid_neuron_flags(i)  = true;
        
end


save('H01_pyramidal_spine_morphology.mat', ...
     'db_normalized_features', 'db_spine_base_ids', 'db_valid_neuron_flags', '-v7.3');
