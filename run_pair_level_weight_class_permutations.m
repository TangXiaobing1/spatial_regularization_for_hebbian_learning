function run_pair_level_weight_class_permutations(varargin)
% Five pooled-pair quintiles plus Top10 under branch/global shuffles.

p = inputParser;
addParameter(p, 'CellTypes', {'pyramidal', 'interneuron', ...
    'spiny_stellate_neuron', 'excitatory_spiny_neuron'});
addParameter(p, 'NShuffles', 1000);
addParameter(p, 'Workers', 16);
addParameter(p, 'Seed', 20260804);
addParameter(p, 'Overwrite', false);
addParameter(p, 'OutputDir', ...
    ['/home/ubuntu/Documents/Neuron_morphology_function_and_connectivity_analysis/' ...
    'data/H01/pair_level_weight_classes_20260804']);
parse(p, varargin{:});
options = p.Results;

addpath(genpath('/home/ubuntu/Documents/Neuron_morphology_function_and_connectivity_analysis'));
if ~exist(options.OutputDir, 'dir'), mkdir(options.OutputDir); end
base_dir = ['/home/ubuntu/Documents/Neuron_morphology_function_and_connectivity_analysis/' ...
    'data/H01'];
specs = struct( ...
    'label', {'pyramidal', 'interneuron', 'spiny_stellate_neuron', ...
    'excitatory_spiny_neuron'}, ...
    'tree_file', {'H01_level1_pyramidal_rearranged.mat', ...
    'H01_level1_interneuron_rearranged.mat', ...
    'H01_level1_spiny_stellate_neuron_rearranged.mat', ...
    'H01_level1_excitatory_spiny_neuron_with_atypical_tree_rearranged.mat'});
specs = specs(ismember({specs.label}, options.CellTypes));
class_labels = {'Q00_20', 'Q20_40', 'Q40_60', 'Q60_80', 'Q80_100', 'Top10'};

pool = gcp('nocreate');
if isempty(pool)
    parpool('local', options.Workers);
elseif pool.NumWorkers ~= options.Workers
    delete(pool);
    parpool('local', options.Workers);
end

for spec_index = 1:numel(specs)
    spec = specs(spec_index);
    output_path = fullfile(options.OutputDir, ...
        sprintf('weight_class_permutations_%s.mat', spec.label));
    summary_path = fullfile(options.OutputDir, ...
        sprintf('weight_class_summary_%s.csv', spec.label));
    if isfile(output_path) && ~options.Overwrite
        fprintf('Skipping existing result: %s\n', output_path);
        continue;
    end

    intermediate = load(fullfile(base_dir, 'data', ...
        sprintf('mexicanhat_data_H01_%s.mat', spec.label)), ...
        'j_values', 'all_probability_data');
    param_index = find(abs(intermediate.j_values - 0.20) < 1e-12, 1);
    pair_data = intermediate.all_probability_data.(sprintf('param_%d', param_index));
    distances = pair_data.global_dists(:);
    volumes = pair_data.global_real_vols(:);
    thresholds = prctile(volumes, [20, 40, 60, 80, 90]);
    zone = distances >= 2 & distances < 14;
    reference_pair_count = sum(zone);
    reference_counts = zeros(1, 6);
    for class_index = 1:6
        reference_counts(class_index) = sum( ...
            class_mask(volumes(zone), thresholds, class_index));
    end
    reference_observed = reference_counts ./ reference_pair_count;
    clear intermediate pair_data distances volumes;

    source = load(fullfile(base_dir, spec.tree_file), 'trees_rearranged_cell');
    source_indices = find(~cellfun('isempty', source.trees_rearranged_cell));
    trees = source.trees_rearranged_cell(source_indices);
    clear source;
    n_neurons = numel(trees);
    fprintf('%s: %d neurons, %d branch + %d global shuffles.\n', ...
        spec.label, n_neurons, options.NShuffles, options.NShuffles);
    observed_pair_count = zeros(n_neurons, 1);
    observed_class_count = zeros(n_neurons, 6);
    branch_shuffle_class_count = zeros(n_neurons, options.NShuffles, 6);
    global_shuffle_class_count = zeros(n_neurons, options.NShuffles, 6);
    status = cell(n_neurons, 1);
    error_message = cell(n_neurons, 1);
    agglomeration_id = cell(n_neurons, 1);
    tree_constant = parallel.pool.Constant(trees);
    clear trees;

    parfor neuron = 1:n_neurons
        try
            local = analyze_pair_level_weight_class_permutations( ...
                tree_constant.Value{neuron}, thresholds, options.NShuffles, ...
                options.Seed + 1000000 + neuron);
            observed_pair_count(neuron) = local.observed_pair_count;
            observed_class_count(neuron, :) = local.observed_class_count;
            branch_shuffle_class_count(neuron, :, :) = reshape( ...
                local.branch_shuffle_class_count, ...
                [1, options.NShuffles, 6]);
            global_shuffle_class_count(neuron, :, :) = reshape( ...
                local.global_shuffle_class_count, ...
                [1, options.NShuffles, 6]);
            status{neuron} = local.status;
            agglomeration_id{neuron} = local.agglomeration_id;
        catch exception
            status{neuron} = 'error';
            error_message{neuron} = getReport(exception, 'extended', ...
                'hyperlinks', 'off');
        end
        if mod(neuron, 100) == 0
            fprintf('%s completed neuron %d/%d\n', spec.label, neuron, n_neurons);
        end
    end
    clear tree_constant;
    assert(all(strcmp(status, 'ok')), '%s contained analysis errors.', spec.label);
    total_pairs = sum(observed_pair_count);
    observed_counts = sum(observed_class_count, 1);
    assert(total_pairs == reference_pair_count);
    assert(all(observed_counts == reference_counts));
    assert(sum(observed_counts(1:5)) == total_pairs);
    observed_statistics = observed_counts ./ total_pairs;

    summary_rows = {};
    branch_statistics = squeeze(sum(branch_shuffle_class_count, 1)) ./ total_pairs;
    global_statistics = squeeze(sum(global_shuffle_class_count, 1)) ./ total_pairs;
    null_names = {'branch_preserving', 'within_neuron_global'};
    null_matrices = {branch_statistics, global_statistics};
    for null_index = 1:2
        null_values = null_matrices{null_index};
        for class_index = 1:6
            samples = null_values(:, class_index);
            observed = observed_statistics(class_index);
            n_less = sum(samples < observed);
            p_value = n_less / options.NShuffles;
            quantiles = prctile(samples, [2.5, 50, 97.5]);
            absolute_suppression = mean(samples) - observed;
            relative_suppression = 100 * absolute_suppression / mean(samples);
            summary_rows(end + 1, :) = {spec.label, null_names{null_index}, ...
                class_labels{class_index}, total_pairs, observed, ...
                mean(samples), std(samples), quantiles(1), quantiles(2), ...
                quantiles(3), n_less, options.NShuffles, p_value, ...
                absolute_suppression, relative_suppression}; %#ok<AGROW>
        end
    end
    summary_table = cell2table(summary_rows, 'VariableNames', { ...
        'cell_type', 'shuffle_type', 'weight_class', 'n_pairs_2_14um', ...
        'observed_statistic', 'shuffle_mean', 'shuffle_sd', 'shuffle_q025', ...
        'shuffle_median', 'shuffle_q975', 'n_shuffle_less_than_observed', ...
        'n_shuffles', 'p_value_strict_lower_tail', ...
        'absolute_suppression', 'relative_suppression_percent'});
    writetable(summary_table, summary_path);
    metadata = struct('zone_um', [2, 14], 'center_top_percent', 20, ...
        'thresholds', thresholds, 'class_labels', {class_labels}, ...
        'n_shuffles_per_null', options.NShuffles, 'seed', options.Seed, ...
        'p_value', 'sum(shuffle_statistic < observed_statistic) / n_shuffles');
    save(output_path, 'metadata', 'source_indices', 'agglomeration_id', ...
        'status', 'error_message', 'observed_pair_count', ...
        'observed_class_count', 'branch_shuffle_class_count', ...
        'global_shuffle_class_count', 'branch_statistics', ...
        'global_statistics', 'observed_statistics', 'summary_table', '-v7.3');
    fprintf('Saved %s and %s\n', output_path, summary_path);
end
end


function flags = class_mask(values, thresholds, class_index)
switch class_index
    case 1, flags = values <= thresholds(1);
    case 2, flags = values > thresholds(1) & values <= thresholds(2);
    case 3, flags = values > thresholds(2) & values <= thresholds(3);
    case 4, flags = values > thresholds(3) & values <= thresholds(4);
    case 5, flags = values > thresholds(4);
    case 6, flags = values > thresholds(5);
end
end
