#Amino acid...
const AA_DICT = Dict('M'=>11,'I'=>8,'Y'=>20,'L'=>10,'*'=>21,'F'=>5,'Q'=>14,
'D'=>3,'V'=>18,'E'=>4,'T'=>17,'H'=>7,'P'=>13,'G'=>6,'N'=>12,'K'=>9,'C'=>2,
'R'=>15,'W'=>19,'A'=>1,'S'=>16)

#possibly add ambig handling for aas...

#Nucleotide...
const NT_DICT = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)

#disambiguate! 
const IUPACbool = Dict{Char,Array{Bool,1}}(Dict())
IUPACbool['A']=[true,false,false,false]
IUPACbool['C']=[false,true,false,false]
IUPACbool['G']=[false,false,true,false]
IUPACbool['T']=[false,false,false,true]
IUPACbool['U']=[false,false,false,true]
IUPACbool['R']=[true,false,true,false]
IUPACbool['Y']=[false,true,false,true]
IUPACbool['S']=[false,true,true,false]
IUPACbool['W']=[true,false,false,true]
IUPACbool['K']=[false,false,true,true]
IUPACbool['M']=[true,true,false,false]
IUPACbool['B']=[false,true,true,true]
IUPACbool['D']=[true,false,true,true]
IUPACbool['H']=[true,true,false,true]
IUPACbool['V']=[true,true,true,false]
IUPACbool['N']=[true,true,true,true];

function resolve_base(c)
    uc = uppercase(c)
    if uc in keys(IUPACbool)
        return ['A','C','G','T'][rand((1:4)[IUPACbool[uc]])]
    else
        return uc
    end
end

function resolve_seq(s)
    return join(resolve_base.(collect(s)))
end

#encode string
function string2encoding(seq::String, lookup_dic::Dict{Char,Int})
    return [lookup_dic[c] for c in seq]
end

