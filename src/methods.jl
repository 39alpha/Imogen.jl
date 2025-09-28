abstract type Method end

abstract type Approximate <: Method end
abstract type Exact <: Method end

abstract type Kozachenko <: Method end

abstract type Kraskov <: Method  end
abstract type Kraskov1 <: Kraskov end
abstract type Kraskov2 <: Kraskov end

struct ExactResult
    coeff :: Vector{Int}
    logarg :: Vector{Int}
    divisor :: Int

    function ExactResult(coeff, logarg, divisor)
        g = gcd(divisor, coeff...)

        keep = (coeff .!= 0) .& (logarg .!= 0)

        logarg = logarg[keep]
        coeff = coeff[keep]

        perm = sortperm(logarg)

        logarg = logarg[perm]
        coeff = coeff[perm]

        coeff .÷= g
        divisor ÷= g

        new(coeff, logarg, divisor)
    end
end

ExactResult(d::Dict{Int,Int}, N::Int) = ExactResult(collect(values(d)), collect(keys(d)), N)

ExactResult(x::Int) = ExactResult([x], [2], 1)
ExactResult(x::Rational) = ExactResult([numerator(x)], [2], denominator(x))

approximate(r::ExactResult) = dot(r.coeff, log2.(r.logarg)) / r.divisor

Base.zero(::Type{ExactResult}) = ExactResult([], [], 1)
Base.one(::Type{ExactResult}) = ExactResult([1], [2], 1)

Base.iszero(a::ExactResult) = isempty(a.logargs)
Base.isone(a::ExactResult) = length(a.logargs) == 1 && first(a.logargs) == 1

Base.isless(a::ExactResult, b::ExactResult) = isless(approximate(a), approximate(b))
Base.isless(a::Number, b::ExactResult) = isless(a, approximate(b))
Base.isless(a::ExactResult, b::Number) = isless(approximate(a), b)

function Base.show(io::IO, r::ExactResult)
    positive = Tuple{Int,Int}[]
    negative = Tuple{Int,Int}[]
    for (c, a) in zip(r.coeff, r.logarg)
        if c > 0
            push!(positive, (c, a))
        else
            push!(negative, (-c, a))
        end
    end
    sort!(positive, by=last)
    sort!(negative, by=last, rev=true)

    fmt(c, a) = a == 2 ? "$(c)" : (isone(c) ? "log₂$(a)" : "$(c)log₂$(a)")
    fmt(p) = fmt(p...)

    p = join(map(fmt, positive), " + ")
    n = join(map(fmt, negative), " - ")

    if isempty(p) && isempty(n)
        print(io, 0)
    elseif isone(r.divisor)
        if isempty(p)
            print(io, n)
        elseif isempty(n)
            print(io, p)
        else
            print(io, p, " - ", n)
        end
    elseif isempty(p)
        if length(negative) == 1
            print(io, "-", n, "/", r.divisor)
        else
            print(io, "-(", replace(n, "-" => "+"), ") / ", r.divisor)
        end
    elseif isempty(n)
        if length(positive) == 1
            print(io, p, "/", r.divisor)
        else
            print(io, "(", p, ") / ", r.divisor)
        end
    else
        print(io, "(", p, " - ", n, ") / ", r.divisor)
    end
end

function Base.:+(x::ExactResult, y::ExactResult)
    u = Dict(zip(x.logarg, x.coeff .* y.divisor))
    v = Dict(zip(y.logarg, y.coeff .* x.divisor))
    return ExactResult(mergewith((+), u, v), x.divisor * y.divisor)
end

function Base.:-(x::ExactResult, y::ExactResult)
    u = Dict(zip(x.logarg, x.coeff .* y.divisor))
    v = Dict(zip(y.logarg, -y.coeff .* x.divisor))
    return ExactResult(mergewith((+), u, v), x.divisor * y.divisor)
end

Base.:+(n::Int, x::ExactResult) = x + ExactResult(n)
Base.:+(x::ExactResult, n::Int) = x + ExactResult(n)
Base.:+(r::Rational, x::ExactResult) = x + ExactResult(r)
Base.:+(x::ExactResult, r::Rational) = x + ExactResult(r)

Base.:*(n::Int, x::ExactResult) = ExactResult(n .* x.coeff, x.logarg, x.divisor)
Base.:*(x::ExactResult, n::Int) = n * x

Base.:/(x::ExactResult, n::Int) = ExactResult(x.coeff, x.logarg, n * x.divisor)
