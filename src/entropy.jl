mutable struct Entropy{D} <: InfoDist
    data::Vector{Int}
    bs::NTuple{D, Int}
    N::Int
    function Entropy(b::Int, bs::Int...)
        if b ≤ zero(b) || any(b -> b ≤ zero(b), bs)
            throw(ArgumentError("all bases must be greater than zero"))
        end
        new{length(bs) + 1}(zeros(Int, b * prod(bs)), tuple(b, bs...), 0)
    end
end

function Entropy(xs::AbstractArray{Int,3})
    iszero(length(xs)) && throw(ArgumentError("no observations provided"))
    any(b -> b < 1, xs) && throw(ArgumentError("observations must be 1 or greater"))

    dist = Entropy(maximum(xs; dims=(2,3))...)
    observe!(dist, xs)
end
function Entropy(xs::AbstractArray{Int,2})
    iszero(length(xs)) && throw(ArgumentError("no observations provided"))
    any(b -> b < 1, xs) && throw(ArgumentError("observations must be 1 or greater"))

    dist = Entropy(maximum(xs; dims=2)...)
    observe!(dist, xs)
end
function Entropy(xs::AbstractArray{Int,1})
    iszero(length(xs)) && throw(ArgumentError("no observations provided"))
    any(b -> b < 1, xs) && throw(ArgumentError("observations must be 1 or greater"))

    dist = Entropy(maximum(xs))
    observe!(dist, xs)
end

function observe!(dist::Entropy, xs::AbstractArray{Int,3})
    @views for i in 1:size(xs, 3)
        @inbounds observe!(dist, xs[:,:,i])
    end
    dist
end
function observe!(dist::Entropy, xs::AbstractArray{Int,2})
    dist.N += size(xs, 2)
    @views for t in 1:size(xs, 2)
        @inbounds x = index(xs[:,t], dist.bs)
        dist.data[x] += 1
    end
    dist
end
function observe!(dist::Entropy, xs::AbstractArray{Int,1})
    dist.N += length(xs)
    @views for t in 1:length(xs)
        @inbounds dist.data[xs[t]] += 1
    end
    dist
end

function entropy(::Type{Approximate}, xs::AbstractArray{Int}, N::Int)
    h = N * log2(N)
    for i in eachindex(xs)
        @inbounds n = xs[i]
        if !iszero(n)
            h -= n * log2(n)
        end
    end
    h / N
end

function entropy(::Type{Exact}, xs::AbstractArray{Int}, N::Int)
    d = Dict{Int, Int}()
    for (f, n) in eachfactor(N)
        d[f] = get(d, f, 0) + N * n
    end
    for m in xs, (f, n) in eachfactor(m)
        d[f] = get(d, f, 0) - m * n
    end
    ExactResult(d, N)
end

estimate(::Type{T}, dist::Entropy) where {T<:Method} = entropy(T, dist.data, dist.N)

entropy!(::Type{T}, dist::Entropy, xs::AbstractArray{Int}) where {T<:Method} = estimate(T, observe!(dist, xs))
entropy!(dist::Entropy, xs::AbstractArray{Int}) = entropy!(Approximate, dist, xs)

entropy(::Type{T}, xs::AbstractArray{Int}) where {T<:Method} = estimate(T, Entropy(xs))
entropy(xs::AbstractArray{Int}) = estimate(Approximate, Entropy(xs))

function clear!(dist::Entropy)
    fill!(dist.data, 0)
    dist.N = 0
end

Base.length(dist::Entropy) = prod(size(dist))

Base.size(dist::Entropy) = dist.bs

Base.getindex(dist::Entropy, idx...) = getindex(dist.data, idx...)

function entropy(::Type{Kozachenko}, xs::AbstractMatrix{Float64}; metric::Metric=Euclidean())
    D, T = size(xs)
    δ = minimumdistances(xs; metric=metric)
    D*(mean(log.(δ)) + log(2.0)) + sdterm(D) + eulermascheroni(T)
    D * (mean(log.(δ)) + log(2.0)) + sdterm(D) + eulermascheroni(T)
end

eulermascheroni(x) = (x < zero(x)) ? zero(Float64) : digamma(x) - digamma(1)

function sdterm(dim::Int)
    c = (π / 4.0)^(dim / 2)
    if iseven(dim)
        log(c) - sum(log.(2:(dim ÷ 2)))
    else
        c = (c * 2.0^((dim + 1) / 2.)) / sqrt(π)
        log(c) - sum(log.(3:2:dim))
    end
end

function minimumdistances(xs; metric::Metric=Euclidean())
    b = BallTree(xs)
    _, δ = knn(b, xs, 2)
    first.(δ)
end

entropy(xs::AbstractMatrix{Float64}; kwargs...) = entropy(Kozachenko, xs; kwargs...)

entropy(xs::AbstractVector{Float64}; kwargs...) = entropy(reshape(xs, 1, length(xs)))
