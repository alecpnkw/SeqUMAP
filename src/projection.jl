function project_sequences(seqs::Array{String,1}, ndim::Int;
    k = 2,
    lookup_dic = AA_DICT,
    n_neighbors = 30, 
    min_dist = 1e-3
    )
    vecs = kmer_embed.(seqs, k, kmer_contribute_identity!; lookup_dic = lookup_dic);
    X = hcat(vecs...);
    X = convert(Array{Float64,2}, X)
    proj = umap(X, ndim; n_neighbors = n_neighbors, min_dist = min_dist)
    return proj
end