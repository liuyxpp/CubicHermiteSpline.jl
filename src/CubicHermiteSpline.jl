module CubicHermiteSpline

using ArgCheck
using DelaunayTriangulation

include("univariate.jl")
include("bivariate.jl")

export
    CubicHermiteSplineInterpolation,
    UnivariateCHSInterpolation,
    BivariateCHSInterpolation

export
    interp,
    grad

end # module