#define custom Distances Metric...
#see discussion of req methods here: https://github.com/JuliaStats/Distances.jl/issues/95

struct CorrectedKmer <: Metric 
    k::Int
end 

correctedkmer(a::AbstractArray, b::AbstractArray, k::Int) = CorrectedKmer(k)(a, b)
correctedkmer(a::Number, b::Number, k::Int) = CorrectedKmer(k)(a, b)

#looking back at RAD paper not sure about this... check!
function evaluate(dist::CorrectedKmer,a,b)
    D = sqeuclidean(a,b)
    #dist.N*log(1 -  minimum([D / (dist.N * 2),1])) / - dist.k
    D / (dist.k * (sum(a) + sum(b)))
end

@eval @inline (dist::CorrectedKmer)(a::AbstractArray, b::AbstractArray) = evaluate(dist, a, b)
@eval @inline (dist::CorrectedKmer)(a::Number, b::Number) = evaluate(dist, a, b)

function sequmap(seqs::Array{String,1}, ndim::Int;
    k = 2,
    lookup_dic = AA_DICT,
    n_neighbors = 30, 
    min_dist = 1e-3
    )
    vecs = kmer_embed.(seqs, k, kmer_contribute_identity!; lookup_dic = lookup_dic);
    N = length(vecs[1])
    X = hcat(vecs...);
    X = convert(Array{Float64,2}, X)
    proj = umap(X, ndim; n_neighbors = n_neighbors, min_dist = min_dist, metric = CorrectedKmer(k))
    return proj
end
