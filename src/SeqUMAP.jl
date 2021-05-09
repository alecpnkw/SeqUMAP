module SeqUMAP

using UMAP

include("encoding.jl")
include("embedding.jl")
include("projection.jl")

#encoding.jl...
export AA_DICT,
NT_DICT,
IUPACbool,
resolve_base,
resolve_seq,
string2encoding,

#embedding.jl...
collect_kmers,
kmer_contribute_identity!,
get_kmer_index,
kmer_embed,

#projection.jl...
project_sequences

end