function results = within_branch_permutation_test(tree, target_percentile, neighbor_percentile, ...
        n_permutations, distance_range, P_obs_external)
% WITHIN_BRANCH_PERMUTATION_TEST
% Permutation test for the spatial arrangement of large spines, using the
% one-dimensional "path distance from the soma" (Pvec) axis together with a
% same-branch constraint.
%
% Neighbour definition. Two spines are neighbours if and only if
%   (a) they belong to the same branch returned by find_biological_branches, and
%   (b) the absolute difference of their Pvec values lies in distance_range.
% Within a single branch Pvec increases monotonically along the path, so under
% (a) the quantity |dPvec| equals the arc length between the two attachment
% points along that branch. Without (a), two spines sitting on different branches
% at a similar depth from the soma are treated as neighbours even though they may
% be tens of microns apart.
%
% Null model. Spine volumes are shuffled within each branch, which preserves the
% volume multiset of every branch and therefore the number of "large" spines in
% it. The test therefore asks whether large spines are arranged differently
% within a branch than a random rearrangement of that branch's own volumes.
%
% Input
%   tree                : TREES toolbox tree structure (R==4 marks spine nodes)
%   target_percentile   : percentile range selecting the central spines,
%                         default [80 100], i.e. the top 20% by volume
%   neighbor_percentile : percentile range for qualifying neighbours,
%                         default [90 100], i.e. the top 10%. Passing the same
%                         range as target_percentile gives a single threshold.
%   n_permutations      : number of permutations, default 1000
%   distance_range      : [min max] along-path distance in um, default [4 12]
%   P_obs_external      : optional externally supplied observed statistic
%
% Output: results struct with fields
%   .P_obs / .P_obs_computed / .P_null_distribution / .p_value
%   .target_percentile / .neighbor_percentile / .distance_range
%   .target_threshold / .neighbor_threshold
%   .n_spines / .n_target / .n_permutations
%   .n_branches / .n_branch_rows_valid / .n_branch_rows_reversed
%   .n_branch_rows_malformed / .n_assigned / .n_unassigned
%   .n_target_with_neighbors / .n_neighbor_pairs / .mean_neighbors_per_target
%
% p-value (one-sided):
%   p = ( #{null <= observed} + 1e-10 ) / ( n_permutations + 1e-10 )
% A small p means the observed statistic lies below the null distribution
% (repulsion / dispersion); a large p means it lies above it (clustering, in
% which case report 1 - p). Note that p is not defined when the statistic itself
% is undefined (P_obs is NaN); such cells must be filtered out before use.
%
% Requires: identify_spine_heads_by_diameter, vol_tree, Pvec_tree, idpar_tree,
%           find_biological_branches, assign_spines_to_branches

    %% Default parameters
    if nargin < 2 || isempty(target_percentile)
        target_percentile = [80, 100];
    end
    if nargin < 3 || isempty(neighbor_percentile)
        neighbor_percentile = [90, 100];
    end
    if nargin < 4 || isempty(n_permutations)
        n_permutations = 1000;
    end
    if nargin < 5 || isempty(distance_range)
        distance_range = [4, 12];
    end
    if nargin < 6
        P_obs_external = [];
    end

    %% 1. Spine volumes and positions
    [spine_ids, head_indices] = identify_spine_heads_by_diameter(tree);
    if isempty(head_indices) || isempty(spine_ids)
        error('No spine head could be identified.');
    end

    node_volumes      = vol_tree(tree);
    head_node_volumes = node_volumes(head_indices);
    head_spine_ids    = spine_ids(head_indices);
    num_spines_all    = max(spine_ids);

    spine_volumes = accumarray(head_spine_ids, head_node_volumes, [num_spines_all, 1], @sum, 0);

    valid_spine_idx = spine_volumes > 0;
    spine_volumes   = spine_volumes(valid_spine_idx);
    n_spines        = numel(spine_volumes);
    if n_spines < 10
        error('Too few valid spines (fewer than 10).');
    end

    parents              = idpar_tree(tree);
    path_lengths_to_root = Pvec_tree(tree);

    spine_positions = zeros(n_spines, 1);
    attach_nodes    = zeros(n_spines, 1);
    valid_spine_ids = find(valid_spine_idx);

    for s = 1:n_spines
        spine_id           = valid_spine_ids(s);
        all_nodes_in_spine = find(spine_ids == spine_id);

        % anchor node: the base of the spine, i.e. the first node whose parent
        % belongs to the dendrite rather than to this spine
        anchor_node = [];
        for i = 1:numel(all_nodes_in_spine)
            curr = all_nodes_in_spine(i);
            p    = parents(curr);
            if p > 0 && spine_ids(p) ~= spine_id
                anchor_node = curr;
                break;
            end
        end
        if isempty(anchor_node)
            anchor_node = all_nodes_in_spine(1);
        end

        dendrite_parent = parents(anchor_node);
        if dendrite_parent > 0 && dendrite_parent <= numel(path_lengths_to_root)
            spine_positions(s) = path_lengths_to_root(dendrite_parent);
            attach_nodes(s)    = dendrite_parent;
        else
            spine_positions(s) = path_lengths_to_root(anchor_node);
            attach_nodes(s)    = anchor_node;
        end
    end

    %% 2. Branch membership
    biological_sect = find_biological_branches(tree);
    if isempty(biological_sect)
        error('find_biological_branches returned no branch; the same-branch constraint cannot be applied.');
    end
    [branch_label, binfo] = assign_spines_to_branches(tree, attach_nodes, biological_sect);

    %% 3. Thresholds and labels
    target_threshold   = prctile(spine_volumes, target_percentile(1));
    neighbor_threshold = prctile(spine_volumes, neighbor_percentile(1));

    is_target          = spine_volumes >= target_threshold;
    is_neighbor_target = spine_volumes >= neighbor_threshold;
    n_target           = sum(is_target);

    if n_target < 2
        error('Too few target synapses (fewer than 2).');
    end

    %% 4. Neighbour adjacency matrix
    % Positions and branch labels are invariant under permutation; only the set
    % of "large" spines changes. The adjacency therefore has to be built once,
    % which reduces every permutation to a single matrix-vector product.
    nbr = false(n_spines, n_spines);
    for i = 1:n_spines
        if branch_label(i) == 0
            continue;                       % unassigned spines have no same-branch neighbour
        end
        dd = abs(spine_positions - spine_positions(i));
        nbr(i, :) = (branch_label == branch_label(i)) & ...
                    (dd >= distance_range(1)) & (dd <= distance_range(2));
    end
    nbr(1:(n_spines + 1):end) = false;      % remove self-loops

    tgt   = find(is_target);
    Nmat  = double(nbr(tgt, :));
    denom = sum(Nmat, 2);                   % number of same-branch neighbours per target spine

    %% 5. Observed statistic
    P_obs_computed = stat_from_matrix(Nmat, denom, is_neighbor_target);
    if ~isempty(P_obs_external)
        P_obs = P_obs_external;
    else
        P_obs = P_obs_computed;
    end

    ok_center = denom >= 1;

    %% 6. Within-branch permutation null model
    % Shuffling within branches is a global permutation of the volume vector, so
    % the multiset of volumes is unchanged and the percentile thresholds computed
    % on the shuffled data would be identical to the observed ones. Recomputing
    % them would therefore be a no-op, and the observed thresholds are reused.
    P_null_distribution = zeros(n_permutations, 1);
    for perm = 1:n_permutations
        volumes_shuffled = shuffle_within_branch_groups(spine_volumes, branch_label);
        is_neighbor_shuffled = volumes_shuffled >= neighbor_threshold;
        P_null_distribution(perm) = stat_from_matrix(Nmat, denom, is_neighbor_shuffled);
    end

    %% 7. p-value and output
    eps_val = 1e-10;
    p_value = (sum(P_null_distribution <= P_obs) + eps_val) / (n_permutations + eps_val);

    results = struct();
    results.P_obs               = P_obs;
    results.P_obs_computed      = P_obs_computed;
    results.P_null_distribution = P_null_distribution;
    results.p_value             = p_value;
    results.target_percentile   = target_percentile;
    results.neighbor_percentile = neighbor_percentile;
    results.distance_range      = distance_range;
    results.target_threshold    = target_threshold;
    results.neighbor_threshold  = neighbor_threshold;
    results.n_spines            = n_spines;
    results.n_target            = n_target;
    results.n_permutations      = n_permutations;

    results.n_branches              = binfo.n_branches;
    results.n_branch_rows_valid     = binfo.n_rows_valid;
    results.n_branch_rows_reversed  = binfo.n_rows_reversed;
    results.n_branch_rows_malformed = binfo.n_rows_malformed;
    results.n_assigned              = binfo.n_assigned;
    results.n_unassigned            = binfo.n_unassigned;
    results.n_target_with_neighbors = sum(ok_center);
    results.n_neighbor_pairs        = sum(denom);
    results.mean_neighbors_per_target = sum(denom) / max(1, sum(ok_center));
end

% -------------------------------------------------------------------------
function P = stat_from_matrix(Nmat, denom, is_neighbor)
% Statistic of one permutation: mean over target spines of the proportion of
% their neighbours that qualify as large.
    num = Nmat * double(is_neighbor(:));
    prop = nan(numel(denom), 1);
    ok = denom >= 1;
    prop(ok) = num(ok) ./ denom(ok);
    P = local_nanmean(prop);
end

function m = local_nanmean(x)
    x = x(~isnan(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

function v = shuffle_within_branch_groups(v, branch_label)
% Shuffle spine volumes within every branch; spines with label 0 stay untouched.
    groups = unique(branch_label(branch_label > 0));
    for g = groups(:)'
        idx = find(branch_label == g);
        if numel(idx) >= 2
            v(idx) = v(idx(randperm(numel(idx))));
        end
    end
end
