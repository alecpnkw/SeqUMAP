#break into kmers for strings
function collect_kmers(seq::String, k::Int, lookup_dic::Dict{Char,Int})
    return [[lookup_dic[c] for c in seq[1+i:k+i]] for i in 0:length(seq)-k]
end

#weights vector length is the number of unique kmers
#each kmer is counted once, contributes to one pos
function kmer_contribute_identity!(kmer_counts::Array{UInt64,1}, kmer::Vector{Int}, n_chars::Int)
    kmer_counts[get_kmer_index(kmer,n_chars)] += unsigned(1)
end

#this function gets the index in lexicographic order...
#remember 1-based indexing in Julia... 
function get_kmer_index(kmer::Vector{Int}, n_chars::Int)
    kmer_ix = 0
    for c in kmer
        kmer_ix = kmer_ix * n_chars + (c - 1)
    end
    return kmer_ix + 1
end 

#with a `kmer_contribution()` function... 
function kmer_embed(seq, k::Int, kmer_contribution::Function; lookup_dic = AA_DICT)
    n = length(union(values(lookup_dic)))
    kmer_counts = zeros(UInt64, n^k) 
    
    #break into kmers
    kmers = collect_kmers(seq, k, lookup_dic)

    #store kmer count in array
    for kmer in kmers
        kmer_contribution(kmer_counts, kmer, n)
    end
    return kmer_counts
end
