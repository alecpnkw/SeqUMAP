# SeqUMAP

[Work in progress.] Encode and embed sequences as kmer vectors (both NT and AA) and generate projections using UMAP.jl.

### Setup

```
using Pkg
Pkg.add(PackageSpec(;name="SeqUMAP",url="https://github.com/alecpnkw/SeqUMAP.git"))
```

### Usage

Kmer embed and project sequences using

```
proj = project_sequences(seqs, ndim; k = 2, lookup_dic = NT_DICT)
```

where `seqs` is an array of strings. For ease of use, the following lookup dictionaries are provided: 

```
AA_DICT = Dict(
    'A' => 1, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 
    'G' => 6, 'H' => 7, 'I' => 8, 'K' => 9, 'L' => 10, 
    'M' => 11, 'N' => 12, 'P' => 13, 'Q' => 14, 'R' => 15, 
    'S' => 16, 'T' => 17, 'V' => 18, 'W' => 19, 'Y' => 20, 
    '*' => 21
    )

NT_DICT = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)
```