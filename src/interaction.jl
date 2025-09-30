mutable struct InteractionInfo <: InfoDist
    joint::Array{Int}
    bs::Vector{Vector{Int}}
    N::Int
    function InteractionInfo(b1::Tuple, b2::Tuple, args::Tuple...)
        bs = [b1, b2, args...]
        if any(isempty, bs)
            throw(MethodError(InteractionInfo, tuple(bs...)))
        end
        if any(x -> any(b -> b < 2, x), bs)
            throw(ArgumentError("the support of each random variable must be at least 2"))
        end
        joint = zeros(Int, prod.(bs)...)

        new(joint, collect.(bs), 0)
    end
end

function InteractionInfo(xs::AbstractArray{Int,3}, ys::AbstractArray{Int,3}, args::AbstractArray{Int,3}...; kwargs...)
    if isempty(xs) || isempty(ys) || any(isempty, args)
        throw(ArgumentError("arguments must not be empty"))
    end
    xmax = max.(2, maximum(xs; dims=(2,3)))
    ymax = max.(2, maximum(ys; dims=(2,3)))
    zmax = map(zs -> tuple(max.(2, maximum(zs; dims=(2,3)))...), args)
    observe!(InteractionInfo(tuple(xmax...), tuple(ymax...), zmax...; kwargs...), xs, ys, args...)
end

function InteractionInfo(xs::AbstractArray{Int,2}, ys::AbstractArray{Int,2}, args::AbstractArray{Int,2}...; kwargs...)
    if isempty(xs) || isempty(ys) || any(isempty, args)
        throw(ArgumentError("arguments must not be empty"))
    end
    xmax = max.(2, maximum(xs; dims=2))
    ymax = max.(2, maximum(ys; dims=2))
    zmax = map(zs -> tuple(max.(2, maximum(zs; dims=2))...), args)
    observe!(InteractionInfo(tuple(xmax...), tuple(ymax...), zmax...; kwargs...), xs, ys, args...)
end

function InteractionInfo(xs::AbstractVector{Int}, ys::AbstractVector{Int}, args::AbstractVector{Int}...; kwargs...)
    if isempty(xs) || isempty(ys) || any(isempty, args)
        throw(ArgumentError("arguments must not be empty"))
    end
    xmax = max(2, maximum(xs))
    ymax = max(2, maximum(ys))
    zmax = map(zs -> tuple(max(2, maximum(zs))...), args)
    observe!(InteractionInfo(tuple(xmax...), tuple(ymax...), zmax...; kwargs...), xs, ys, args...)
end

function estimate(::Type{T}, dist::InteractionInfo) where {T<:Method}
    sum((-1)^(length(dims)-1) * entropy(T, sum(dist.joint; dims), dist.N) for dims in subsets0(ndims(dist.joint)))
end

function observe!(dist::InteractionInfo, xs::AbstractArray{Int,3}, ys::AbstractArray{Int,3}, args::AbstractArray{Int,3}...)
    if size(xs, 3) != size(ys, 3) || any(zs -> size(zs, 3) != size(xs, 3), args)
        throw(ArgumentError("time series should have the same number of replicates"))
    end

    @views for i in 1:size(xs, 3)
        observe!(dist, xs[:,:,i], ys[:,:,i], (zs[:,:,i] for zs in args)...)
    end

    dist
end

function observe!(dist::InteractionInfo, xs::AbstractArray{Int,2}, ys::AbstractArray{Int,2}, args::AbstractArray{Int,3}...)
    if size(xs, 2) != size(ys, 2) || any(zs -> size(zs, 2) != size(xs, 2), args)
        throw(ArgumentError("time series should have the same number of timesteps"))
    elseif any(b -> b < 1, xs) || any(b -> b < 1, ys) || any(zs -> any(b -> b < 1, zs), args)
        throw(ArgumentError("observations must be positive, nonzero"))
    end

    N = size(xs, 2)

    dist.N += N
    @views for t in 1:N
        idx = [index(xs[:,t], dist.bs[1]), index(ys[:,t], dist.bs[2]), (index(zs[:,t], dist.bs[2 + i]) for (i, zs) in args)...]
        dist.joint[idx...] += 1
    end
    dist
end

function observe!(dist::InteractionInfo, xs::AbstractArray{Int,1}, ys::AbstractArray{Int,1}, args::AbstractArray{Int,1}...)
    if length(xs) != length(ys) || any(length.(args) .!= length(xs))
        throw(ArgumentError("time series should have the same number of timesteps"))
    elseif any(b -> b < 1, xs) || any(b -> b < 1, ys) || any(zs -> any(b -> b < 1, zs), args)
        throw(ArgumentError("observations must be positive, nonzero"))
    end

    N = length(xs)

    dist.N += N
    @views for t in 1:N
        idx = [xs[t], ys[t], (zs[t] for zs in args)...]
        dist.joint[idx...] += 1
    end
    dist
end

@inline function clear!(dist::InteractionInfo)
    dist.joint[:] .= 0
    dist.N = 0
    dist
end

interaction!(::Type{T}, dist::InteractionInfo, args...; kwargs...) where {T<:Method} = estimate(T, observe!(dist, args...; kwargs...))
interaction!(dist::InteractionInfo, args...; kwargs...) = interaction!(Approximate, dist, args...; kwargs...)
interaction(::Type{T}, args...; kwargs...) where {T<:Method} = estimate(T, InteractionInfo(args...; kwargs...))
interaction(args...; kwargs...) = interaction(Approximate, args...; kwargs...)
