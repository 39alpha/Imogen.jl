abstract type InfoDist end

"""
    observe!(d, args...)
"""
observe!(::InfoDist, args...)

"""
    clear!(d)
"""
clear!(::InfoDist)

"""
    estimate([::Type{<:Method}=Approximate], d)
"""
estimate(::Type{<:Method}, ::InfoDist)
estimate(d::InfoDist) = estimate(Approximate, d)
