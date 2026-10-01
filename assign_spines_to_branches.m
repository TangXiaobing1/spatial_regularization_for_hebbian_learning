function [branch_label, info] = assign_spines_to_branches(tree, attach_nodes, biological_sect)
% ASSIGN_SPINES_TO_BRANCHES
% Map every spine to the index of the biological branch it is attached to, using
% the branch list returned by find_biological_branches.
%
% Note on biological_sect. The matrix is documented as holding [start end] pairs
% in which start is an ancestor of end. In practice, on the cells tested, 11 of
% 103 rows are stored in the reverse order and 5 rows have two endpoints with no
% ancestor relation at all (artefacts of the merging loop inside
% find_biological_branches). This function therefore does not trust the column
% order: it tests the ancestor relation explicitly and derives the node set of
% each branch from the corresponding root path.
%
% Rows whose two endpoints do not lie on a common root-to-leaf path are
% discarded. Expanding such a row would yield a Y-shaped path spanning a fork,
% which would merge sister branches and reintroduce exactly the cross-branch
% pseudo-neighbours that the same-branch constraint is meant to remove.
%
% Only spines whose attachment node lies on a branch segment receive a branch
% index. The remaining spines get 0 and take part in no neighbour relation, i.e.
% they neither act as neighbours nor contribute to the statistic.
%
% Input
%   tree            : TREES toolbox tree structure
%   attach_nodes    : n x 1 dendrite attachment node of every spine
%   biological_sect : m x 2 output of find_biological_branches
%
% Output
%   branch_label : n x 1 branch index (1..m); 0 = not assigned
%   info         : struct with fields
%                  n_branches / n_rows_valid / n_rows_reversed / n_rows_malformed
%                  / n_assigned / n_unassigned / node2branch

    N  = numel(tree.X);
    pp = ipar_tree(tree);
    m  = size(biological_sect, 1);

    % node2branch: label the nodes of every branch segment. Nodes shared by two
    % segments (branch points) are claimed by the first branch that reaches them,
    % which keeps the assignment deterministic.
    node2branch      = zeros(N, 1);
    n_rows_valid     = 0;
    n_rows_reversed  = 0;
    n_rows_malformed = 0;

    for b = 1:m
        a = biological_sect(b, 1);
        c = biological_sect(b, 2);
        if a < 1 || a > N || c < 1 || c > N
            n_rows_malformed = n_rows_malformed + 1;
            continue;
        end

        pa = pp(a, :); pa = pa(pa > 0);     % [a ... root]
        pc = pp(c, :); pc = pc(pc > 0);     % [c ... root]

        if any(pc == a)
            % a is an ancestor of c: segment = c -> a
            ic  = find(pc == a, 1);
            seg = pc(1:ic);
        elseif any(pa == c)
            % c is an ancestor of a (row stored in reverse order): segment = a -> c
            ia  = find(pa == c, 1);
            seg = pa(1:ia);
            n_rows_reversed = n_rows_reversed + 1;
        else
            % the two endpoints share no root path: discard the row
            n_rows_malformed = n_rows_malformed + 1;
            continue;
        end

        n_rows_valid = n_rows_valid + 1;
        for k = 1:numel(seg)
            if node2branch(seg(k)) == 0
                node2branch(seg(k)) = b;
            end
        end
    end

    branch_label = zeros(numel(attach_nodes), 1);
    for s = 1:numel(attach_nodes)
        nd = attach_nodes(s);
        if nd >= 1 && nd <= N
            branch_label(s) = node2branch(nd);
        end
    end

    info = struct();
    info.n_branches       = m;
    info.n_rows_valid     = n_rows_valid;
    info.n_rows_reversed  = n_rows_reversed;
    info.n_rows_malformed = n_rows_malformed;
    info.n_assigned       = sum(branch_label > 0);
    info.n_unassigned     = sum(branch_label == 0);
    info.node2branch      = node2branch;
end
