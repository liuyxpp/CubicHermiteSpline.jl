struct UnivariateCHSInterpolation{T}
    x::Vector{T}
    y::Vector{T}
    dydx::Vector{T}

    function UnivariateCHSInterpolation(x::AbstractVector{T}, y::AbstractVector{T}, dydx::AbstractVector{T}) where T
        @argcheck length(x) == length(y)
        @argcheck length(x) == length(dydx)
        p = sortperm(x)
        return new{T}(collect(x[p]), collect(y[p]), collect(dydx[p]))
    end
end

const CubicHermiteSplineInterpolation = UnivariateCHSInterpolation

function Base.show(io::IO, spl::UnivariateCHSInterpolation)
    n = length(spl.x)
    xmin, xmax = extrema(spl.x)
    ymin, ymax = extrema(spl.y)
    gmin, gmax = extrema(spl.dydx)
    println(io, "Univariate Cubic Hermite Spline Interpolation")
    println(io, "Number of points: ", n)
    print(io, "x: [", xmin, ", ", xmax, "], ")
    println(io, "y: [", ymin, ", ", ymax, "]")
    print(io, "gradient: [", gmin, ", ", gmax, "]")
end

"""
    basis(t)

Compute the basis functions of cubic Hermite spline interpolation. They are:
    H_{00} = 2t^3 - 3t^2 + 1 = (1+2t)(t-1)^2
    H_{10} = t^3 - 2t^2 + t = t(t-1)^2
    H_{01} = -2t^3 + 3t^2 = t^2(3-2t)
    H_{11} = t^3 - t^2 = t^2(t - 1)
"""
function basis(t)
    t2 = t * t
    it = t - 1
    it2 = it * it
    tt = 2 * t
    h00 = (1 + tt) * it2
    h10 = t * it2
    h01 = t2 * (3 - tt)
    h11 = t2 * it
    return h00, h10, h01, h11
end

"""
    basis_derivative(t)

Compute the basis functions of cubic Hermite spline interpolation. They are:
    H_{00} = 6t^2 - 6t = 6t(t-1)
    H_{10} = 3t^2 - 4t + 1 = (3t-1)(t-1)
    H_{01} = -6t^2 + 6t = -6t(t-1)
    H_{11} = 3t^2 - 2t = t(3t - 2)
"""
function basis_derivative(t)
    t2 = t * t
    h00 = 6*t2 - 6*t
    h10 = 3*t2 - 4*t + 1
    h01 = -6*t2 + 6*t
    h11 = 3*t2 - 2*t
    return h00, h10, h01, h11
end

function findinterval(v, x)
    for i in eachindex(x)
        if x[i] == v
            return i, x[i], nothing
        elseif x[i] > v
            return i-1, x[i-1], x[i]
        end
    end
end

function _interp(spl::UnivariateCHSInterpolation, v; grad=false)
    x = spl.x
    y = spl.y
    dydx = spl.dydx
    @argcheck v >= x[1]
    @argcheck v <= x[end]

    idx, x1, x2 = findinterval(v, x)
    if x2 === nothing
       return grad ? dydx[idx] : y[idx]
    end

    # mapping (x1, x2) to (0, 1), h is the scaling constant
    t = (v - x1) / (x2 - x1)
    h = x2 - x1
    y1, y2 = y[idx], y[idx+1]
    k1, k2 = dydx[idx], dydx[idx+1]

    h00, h10, h01, h11 = grad ? basis_derivative(t) : basis(t)
    r = y1*h00 + h*k1*h10 + y2*h01 + h*k2*h11

    # Need to rescale the interpolated gradient before return
    return grad ? (r/h) : r
end

(spl::UnivariateCHSInterpolation)(v::Real; grad=false) = _interp(spl, v; grad=grad)
(spl::UnivariateCHSInterpolation)(x::AbstractVector; grad=false) = spl.(x; grad=grad)

# handy methods
interp(spl::UnivariateCHSInterpolation, p) = spl(p)
grad(spl::UnivariateCHSInterpolation, p) = spl(p; grad=true)