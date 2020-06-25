module CubicHermiteSpline

using ArgCheck

export CubicHermiteSplineInterpolation
export interp

struct CubicHermiteSplineInterpolation
    x
    y
    gradient
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

function findinterval(v, x)
    for i in eachindex(x)
        if x[i] == v
            return i, x[i], nothing
        elseif x[i] > v
            return i-1, x[i-1], x[i]
        end
    end
end

function interp(spl::CubicHermiteSplineInterpolation, v)
    x = spl.x
    y = spl.y
    gradient = spl.gradient
    @argcheck v >= x[1]
    @argcheck v <= x[end]

    idx, x1, x2 = findinterval(v, x)
    if x2 === nothing
       return y[idx]
    end
    t = (v - x1) / (x2 - x1)
    h = x2 - x1
    y1, y2 = y[idx], y[idx+1]
    k1, k2 = gradient[idx], gradient[idx+1]

    h00, h10, h01, h11 = basis(t)
    return y1*h00 + h*k1*h10 + y2*h01 + h*k2*h11
end

(spl::CubicHermiteSplineInterpolation)(x::Real) = interp(spl, x)
(spl::CubicHermiteSplineInterpolation)(x::AbstractVector) = [spl(v) for v in x]

end # module