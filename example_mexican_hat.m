HP = load('H01_pyramidal.mat');


num_neurons = 2992;
j_values = 0.05:0.05:0.25;
num_params = length(j_values);
all_results_cell = cell(num_neurons, num_params);


p = parpool('local', 16);
HP_const = parallel.pool.Constant(HP);

for j_idx = 1:num_params
    current_j = j_values(j_idx);
    fprintf('处理参数 j = %.2f...\n', current_j);
    
    parfor i = 1:num_neurons
        tree = HP_const.Value.trees_rearranged_cell{i};
        data = analyze_mexican_hat_profile_v2(tree, 'DiameterThreshold', 3, 'LargeThreshold', current_j);
        data.original_filename = tree.original_filename;
        all_results_cell(i, j_idx) = {data};
    end
end

clear HP_const;
delete(p);


save('H01_MexicanHat_pyramidal_v2_new_head.mat', 'all_results_cell', 'j_values', '-v7.3');

