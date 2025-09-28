mutable struct SpecificInfo <: InfoDist
    joint::Matrix{Int}
    m1::Vector{Int}
    m2::Vector{Int}
    b1::Int
    b2::Int
    N::Int
    function SpecificInfo(b1::Int, b2::Int)
        if b1 < 1 || b2 < 1
            throw(ArgumentError("the support of each random variable must be at least 2"))
        end
        new(zeros(Int, b1, b2), zeros(Int, b1), zeros(Int, b2), b1, b2, 0)
    end
end

function SpecificInfo(xs::AbstractVector{Int}, ys::AbstractVector{Int})
    if isempty(xs) || isempty(ys)
        throw(ArgumentError("arguments must not be empty"))
    end
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    if xmin < 1 || ymin < 1
        throw(ArgumentError("observations must be positive, nonzero"))
    end
    observe!(SpecificInfo(max(2, xmax), max(2, ymax)), xs, ys)
end

function observe!(dist::SpecificInfo, xs::AbstractVector{Int}, ys::AbstractVector{Int})
    if length(xs) != length(ys)
        throw(ArgumentError("arguments must have the same length"))
    end
    dist.N += length(xs)
    for i in eachindex(xs)
        dist.joint[xs[i], ys[i]] += 1
        dist.m1[xs[i]] += 1
        dist.m2[ys[i]] += 1
    end
    dist
end

function estimate(::Type{Approximate}, dist::SpecificInfo)
    si = zeros(dist.b1)
    N = dist.N
    for x1 in eachindex(dist.m1)
        n1 = dist.m1[x1]
        for x2 in eachindex(dist.m2)
            if !iszero(dist.joint[x1,x2])
                j, n2 = dist.joint[x1,x2], dist.m2[x2]
                si[x1] += j * log2((N * j) / (n1 * n2))
            end
        end
        si[x1] /= n1
    end
    si
end

function estimate(::Type{Exact}, dist::SpecificInfo)
    si = ExactResult[]
    for x1 in eachindex(dist.m1)
        d = Dict{Int,Int}()

        n1 = dist.m1[x1]
        for x2 in eachindex(dist.m2)
            j, n2 = dist.joint[x1,x2], dist.m2[x2]
            num, den = dist.N * j, n1 * n2

            for (f, n) in eachfactor(num)
                d[f] = get(d, f, 0) + j * n
            end
            for (f, n) in eachfactor(den)
                d[f] = get(d, f, 0) - j * n
            end
        end

        push!(si, ExactResult(d, n1))
    end
    si
end

@inline function clear!(dist::SpecificInfo)
    dist.joint[:] .= 0
    dist.m1[:] .= 0
    dist.m2[:] .= 0
    dist.N = 0
    dist
end

function specificinfo!(::Type{T}, si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractVector{Int}) where {T<:Method}
    estimate(T, observe!(dist, stimulus, responses))
end
specificinfo!(si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractVector{Int}) = specificinfo!(Approximate, si, stimulus, responses)

function specificinfo!(::Type{T}, si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}) where {T<:Method}
    specificinfo!(T, si, stimulus, box(responses))
end
specificinfo!(si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}) = specificinfo!(Approximate, si, stimulus, responses)

function specificinfo!(::Type{T}, si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}, subset::AbstractVector{Int}) where {T<:Method}
    specificinfo!(T, si, stimulus, @view responses[subset, :])
end
specificinfo!(si::SpecificInfo, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}, subset::AbstractVector{Int}) = specificinfo!(Approximate, si, stimulus, responses, subset)

function specificinfo(::Type{T}, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}) where {T<:Method}
    specificinfo(T, stimulus, box(responses))
end
specificinfo(stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}) = specificinfo(Approximate, stimulus, responses)

function specificinfo(::Type{T}, stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}, subset::AbstractVector{Int}) where {T<:Method}
    specificinfo(T, stimulus, @view responses[subset, :])
end
specificinfo(stimulus::AbstractVector{Int}, responses::AbstractMatrix{Int}, subset::AbstractVector{Int}) = specificinfo(Approximate, stimulus, responses, subset)

function specificinfo(::Type{T}, stimulus::AbstractVector{Int}, responses::AbstractVector{Int}) where {T<:Method}
    estimate(T, SpecificInfo(stimulus, responses))
end
specificinfo(stimulus::AbstractVector{Int}, responses::AbstractVector{Int}) = specificinfo(Approximate, stimulus, responses)
