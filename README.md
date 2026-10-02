# Topology and Morphology Analysis Toolbox for Neurons

## Overview

This toolbox provides MATLAB functions for analyzing reconstructed neurons. It focuses on dendritic topology, spine morphology, spine distribution, and statistical analysis of spine organization. It also includes utilities for feature extraction, normalization, visualization, and related simulation workflows.

## Core Modules

- **Dendritic topology analysis**: Identify biological branches, branch points, and non-spine termination points.
- **Spine morphology analysis**: Extract, normalize, classify, and summarize spine features.
- **Spine distribution analysis**: Measure spine density, inter-spine intervals, nearest-neighbor distances, and spatial profiles.
- **Statistical analysis**: Test spatial clustering, morphological similarity, and relationships between spine size and spatial organization.
- **Utilities and supporting analyses**: Relabel neuron nodes, identify spine heads, and analyze spine volume distributions.

## Dependencies

- [TREES toolbox](https://www.treestoolbox.org/), version 1.16.

Other versions may work but have not been tested.

## Dendritic Topology Analysis

This module analyzes the neurite backbone while excluding spine branches.

### Find biological branches

```matlab
biological_sect = find_biological_branches(tree, spine_type, min_branch_length)
```

Similar to `dissect_tree` in the TREES toolbox, but designed for neurons with many spines. The `min_branch_length` parameter helps prevent short, incorrectly labeled spine branches from being treated as biological branches.

### Find true branch points

```matlab
is_true_branch_point_vec = find_true_branch_points(tree, spine_type)
```

Similar to `B_tree` in the TREES toolbox, but excludes spine-related branch points.

### Find non-spine termination points

```matlab
is_true_termination_point_vec = find_non_spine_termination_points(tree, spine_type)
```

Similar to `T_tree` in the TREES toolbox, but excludes spine termination points.

## Spine Morphology Analysis

### Calculate spine morphological features

```matlab
neck_diameters = calculate_spine_neck_diameter(tree, spine_type)

[head_volumes_terminal, spine_volumes_total, dist_from_root] = ...
    calculate_spine_head_volumes(tree, spine_type)

[normalized_features, raw_features_matrix, spine_ids_in_order, ...
    feature_names, feature_mean, feature_std] = ...
    extract_and_normalize_spine_features(tree)
```

`extract_and_normalize_spine_features` extracts six morphological features—length, head diameter, neck diameter, surface area, volume, and head-to-neck ratio—and returns their z-score-normalized values.

### Extract spine-head features

```matlab
head_volumes = extract_spine_head_volume_feature(tree, spine_type)
```

Extracts spine-head volume features for downstream analysis.

### Relabel neuron nodes

```matlab
relabeled_tree = relabel_neuron_nodes(original_tree)

relabeled_tree = relabel_neuron_nodes_v6(original_tree)
```

These functions classify unlabeled nodes for further analysis. The `v6` version additionally identifies axon nodes and uses the following labels:

- `1`: soma
- `2`: dendrite
- `3`: axon
- `4`: spine

### Identify spine heads by diameter

```matlab
spine_head_ids = identify_spine_heads_by_diameter(tree, spine_type, threshold)
```

Identifies spine heads using a diameter-based criterion.

### Count nodes per spine

```matlab
[total_spine_count, nodes_per_spine, spine_base_nodes] = ...
    analyze_spine_node_num(relabeled_tree)
```

Some spines contain multiple nodes, particularly in the H01 dataset. Therefore, the number of nodes labeled as spines may exceed the number of individual spines. This function also identifies each spine's base node on the neuron backbone.

## Spine Distribution Analysis

### Calculate spine density

```matlab
[total_path_length, total_spines_on_path, overall_spine_density, bins_info] = ...
    calculate_spine_density_on_path(tree, path_start_node, path_end_node, bin_length)

branch_info = calculate_spine_density_per_branch(tree, spine_type)
```

The first function calculates spine density along a path. The second summarizes spine density by branch order.

### Calculate inter-spine intervals

```matlab
all_isi_values = calculate_isi(tree, spine_type)
```

Calculates inter-spine intervals (ISI).

### Calculate nearest-neighbor distances

```matlab
mean_nnd = calculate_nnd_for_subset(tree, subset_spine_indices, biological_sect)
```

Calculates the mean nearest-neighbor distance (NND) for a selected spine set.

### Analyze spatial organization

```matlab
spines_sorted_by_branch = get_spines_sorted_by_branch(...)

profile_data = analyze_mexican_hat_profile(...)

[results_struct, elbow_radius] = analyze_spatial_window_similarity(...)

[results_struct, elbow_radius] = analyze_spatial_window_similarity_1d(...)

results = analyze_synapse_weight_spatial_clustering(...)
```

These functions sort spines along biological branches and analyze spatial relationships using branch topology, physical-distance windows, spine-head volume, and density-matched random controls.

### Analyze spine volume distributions

```matlab
[U, S, V, volume_matrix_normalized, bin_edges] = ...
    analyze_spine_volume_svd(tree, spine_type, bin_length, min_branch_length)
```

Aligns branches by physical distance, bins spine-head volume along each branch, normalizes the resulting matrix, and performs singular value decomposition (SVD).

## Statistical Analysis

### Find spine hotspots

```matlab
hotspot_info = find_spine_hotspots(tree, spine_type)
```

Identifies spatial regions with concentrated spines. The repository also contains a `find_spine_hotspots_v3` implementation.

### Measure within-cluster similarity

```matlab
result = calculate_intra_cluster_similarity(feature_matrix_for_cluster, varargin)
```

Calculates average pairwise Euclidean distance or average variance within clusters.

### Test spatial clustering

```matlab
[nnd_distribution, observed_nnd, p_value] = ...
    permutation_test_nnd(tree, target_spine_indices, all_spine_indices, num_permutations)

results = within_branch_permutation_test(tree, ...)

results = within_branch_permutation_test_v2(tree, ...)

run_pair_level_weight_class_permutations(...)
```

These functions use permutation tests to compare observed spatial organization with random or within-branch null models.

`run_pair_level_weight_class_permutations` tests pooled pair-level spine-volume classes around fixed observed large-synapse centers using branch-restricted and neuron-wide volume permutations.

### Test cluster homogeneity

```matlab
results = test_cluster_homogeneity(...)

results = compare_spatial_homogeneity(...)

results = compare_spatial_homogeneity_v2(...)
```

These functions test whether spatially clustered spines are morphologically more similar than randomly selected spines. The `v2` comparison includes global-random, branch-random, and branch-adjacent models.

### Test spine-head similarity

```matlab
results = verify_spine_head_similarity(tree, spine_type, max_cluster_size)
```

Uses spine-head volume as a proxy for synaptic weight and tests whether adjacent synapses tend to have similar weights.

### Fit ISI distributions

```matlab
[pdTable, bestPd] = compute_ISI_models(ISI, plotFlag)
```

Fits candidate probability distributions to inter-spine interval data.

## Supporting Scripts

`Subtree_distribution_with_different_types_of_spine.m` compares subtree distributions across spine categories.
`example_spine_morphology.m` and `example_mexican_hat_profile.m` give a quick start.

## Data and File Conventions

Most functions operate on tree structures from the TREES toolbox. The expected node labels depend on the function; when relabeled trees are used, the `R` field commonly contains soma, dendrite, axon, and spine labels `1`, `2`, `3`, and `4`, respectively.
