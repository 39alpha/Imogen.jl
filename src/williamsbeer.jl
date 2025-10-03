using Printf

mutable struct WilliamsBeer
    Imin::Float64
    Π::Float64
end
Base.zero(::Type{WilliamsBeer}) = WilliamsBeer(0.0, 0.0)

function Base.isequal(p::WilliamsBeer, q::WilliamsBeer)
    isequal(p.Imin, q.Imin) && isequal(p.Π, q.Π)
end

Base.:(==)(p::WilliamsBeer, q::WilliamsBeer) = p.Imin == q.Imin && p.Π == q.Π

Base.:≈(p::WilliamsBeer, q::WilliamsBeer) = p.Imin ≈ q.Imin && p.Π ≈ q.Π

Base.iszero(p::WilliamsBeer) = p.Π == zero(p.Π)

function Base.show(io::IO, ::MIME"text/dot", p::WilliamsBeer)
    @printf io "Π = %0.3f\\nImin = %0.3f" p.Π p.Imin
end

function pid!(h::Hasse{<:AbstractVertex{WilliamsBeer}},
             stimulus::AbstractVector{Int},
             responses::AbstractMatrix{Int};
             zero=true)
    if zero
        zero!(h)
    end

    bs = maximum(stimulus)

    L = size(responses, 1)
    ss = subsets(L)
    si = Vector{Vector{Float64}}(undef, length(ss))
    for i in eachindex(ss)
        si[i] = specificinfo(Approximate, stimulus, responses, ss[i])
    end

    sdist = observe!(Entropy(maximum(stimulus)), stimulus)

    for i in eachindex(h)
        α = h[i]
        for s in 1:bs
            x = si[α[1]][s]
            for k in 2:length(α)
                x = min(x, si[α[k]][s])
            end
            payload(α).Imin += sdist[s] * x;
        end
        payload(α).Imin /= sdist.N
    end

    for i in eachindex(h)
        α = h[i]
        for s in 1:bs
            u = -Inf
            for β in below(α)
                x = si[β[1]][s]
                for k in 2:length(β)
                    x = min(x, si[β[k]][s])
                end
                u = max(u, x)
            end
            if isinf(u)
                u = 0.0
            end
            payload(α).Π += sdist[s] * u
        end
        payload(α).Π = payload(α).Imin - payload(α).Π / sdist.N
    end

    h
end

mutable struct ExactWilliamsBeer
    Imin::ExactResult
    Π::ExactResult
end
Base.zero(::Type{ExactWilliamsBeer}) = ExactWilliamsBeer(zero(ExactResult), zero(ExactResult))

function Base.isequal(p::ExactWilliamsBeer, q::ExactWilliamsBeer)
    isequal(p.Imin, q.Imin) && isequal(p.Π, q.Π)
end

Base.:(==)(p::ExactWilliamsBeer, q::ExactWilliamsBeer) = p.Imin == q.Imin && p.Π == q.Π

Base.:≈(p::ExactWilliamsBeer, q::ExactWilliamsBeer) = p.Imin ≈ q.Imin && p.Π ≈ q.Π

Base.iszero(p::ExactWilliamsBeer) = p.Π == zero(p.Π)

function Base.show(io::IO, ::MIME"text/dot", p::ExactWilliamsBeer)
    @printf io "Π = %s\\nImin = %s" p.Π p.Imin
end

function pid!(h::Hasse{<:AbstractVertex{ExactWilliamsBeer}},
             stimulus::AbstractVector{Int},
             responses::AbstractMatrix{Int};
             zero=true)
    if zero
        zero!(h)
    end

    bs = maximum(stimulus)

    L = size(responses, 1)
    ss = subsets(L)
    si = Vector{Vector{ExactResult}}(undef, length(ss))
    for i in eachindex(ss)
        si[i] = specificinfo(Exact, stimulus, responses, ss[i])
    end

    sdist = observe!(Entropy(maximum(stimulus)), stimulus)

    for i in eachindex(h)
        α = h[i]
        for s in 1:bs
            x = si[α[1]][s]
            for k in 2:length(α)
                x = min(x, si[α[k]][s])
            end
            payload(α).Imin += sdist[s] * x;
        end
        payload(α).Imin /= sdist.N;
    end

    for i in eachindex(h)
        α = h[i]
        for s in 1:bs
            us = ExactResult[]
            for β in below(α)
                x = si[β[1]][s]
                for k in 2:length(β)
                    x = min(x, si[β[k]][s])
                end
                push!(us, x)
            end
            if !isempty(us)
                payload(α).Π += sdist[s] * maximum(us)
            end
        end
        payload(α).Π = payload(α).Imin - payload(α).Π / sdist.N
    end

    h
end
