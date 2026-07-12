function avg_pairwise_dist = calculate_intra_cluster_similarity(feature_matrix_for_cluster)

if size(feature_matrix_for_cluster, 1) < 2
    avg_pairwise_dist = 0;
    return;
end

pairwise_distances = pdist(feature_matrix_for_cluster, 'euclidean');

avg_pairwise_dist = mean(pairwise_distances);
end