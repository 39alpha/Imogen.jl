module Imogen

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 3
end

using DataFrames
using Distances
using LinearAlgebra
using NearestNeighbors
using Primes
using Reexport
using SpecialFunctions
using Statistics

@reexport using Random

export Method, Approximate, Exact, Kozachenko, Kraskov, Kraskov1, Kraskov2
export ExactResult, approximate
export InfoDist, observe!, clear!, estimate
export Entropy, entropy!, entropy
export MutualInfo, mutualinfo!, mutualinfo
export TotalCorrelation, totalcor!, totalcor
export InteractionInfo, interaction!, interaction
export ActiveInfo, activeinfo!, activeinfo
export SpecificInfo, specificinfo!, specificinfo
export TransferEntropy, transferentropy!, transferentropy
export AbstractVertex, AbstractUnnamedVertex, AbstractNamedVertex, id, name, payload, above, below
export UnnamedVertex, Vertex, clone
export Hasse, top, bottom, vertices, zero!, prune, graphviz
export WilliamsBeer, ExactWilliamsBeer, pid, pid!
export Significance, EmpiricalSig, AnalyticSig, sig, @sig
export histories, box, encodehistories

include("core.jl")
include("methods.jl")

include("util.jl")

include("entropy.jl")

include("mi.jl")
include("totalcor.jl")
include("interaction.jl")

include("ai.jl")

include("si.jl")

include("te.jl")

include("lattice.jl")
include("pid.jl")
include("williamsbeer.jl")

include("sig.jl")

end
