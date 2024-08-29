struct BivariateCHSInterpolation{T, R}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    dzdx::Vector{T}
    dzdy::Vector{T}
    triangles::R
end

function BivariateCHSInterpolation(x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, dzdx::AbstractVector{T}, dzdy::AbstractVector{T}) where T
    @argcheck length(x) == length(y)
    @argcheck length(x) == length(z)
    @argcheck length(x) == length(dzdx)
    @argcheck length(x) == length(dzdy)
    triangles = triangulate([x'; y'])
    return BivariateCHSInterpolation(collect(x), collect(y), collect(z), collect(dzdx), collect(dzdy), triangles)
end

function Base.show(io::IO, spl::BivariateCHSInterpolation)
    n = length(spl.x)
    xmin, xmax = extrema(spl.x)
    ymin, ymax = extrema(spl.y)
    fmin, fmax = extrema(spl.z)
    fxmin, fxmax = extrema(spl.dzdx)
    fymin, fymax = extrema(spl.dzdy)
    nt = num_triangles(spl.triangles)
    println(io, "Bivariate Cubic Hermite Spline Interpolation")
    println(io, "Number of points: ", n)
    println(io, "Number of Delaunay triangles: ", nt)
    print(io, "x: [", xmin, ", ", xmax, "], ")
    println(io, "y: [", ymin, ", ", ymax, "]")
    println(io, "f: [", fmin, ", ", fmax, "]")
    print(io, "fx: [", fxmin, ", ", fxmax, "], ")
    println(io, "fy: [", fymin, ", ", fymax, "]")
end

function cartesian2baricentric(triangle, x, y)
    p, q, r = triangle
    x1, y1 = p
    x2, y2 = q
    x3, y3 = r
    S1 = x2*y3 - x3*y2
    S2 = x3*y1 - x1*y3
    S3 = x1*y2 - x2*y1
    S = S1 + S2 + S3
    l1 = (S1 + (y2 - y3)*x + (x3 - x2)*y)/S
    l2 = (S2 + (y3 - y1)*x + (x1 - x3)*y)/S
    l3 = (S3 + (y1 - y2)*x + (x2 - x1)*y)/S
    return l1, l2, l3
end

function baricentric2cartesian(triangle, l1, l2, l3)
    p, q, r = triangle
    x1, y1 = p
    x2, y2 = q
    x3, y3 = r
    x = l1*x1 + l2*x2 + l3*x3
    y = l1*y1 + l2*y2 + l3*y3
    return x, y
end

function bchs_coefficients(l1, l2, l3, triangle)
    p, q, r = triangle
    x1, y1 = p
    x2, y2 = q
    x3, y3 = r
    l111 = l1 * l1 * l1
    l112 = l1 * l1 * l2
    l113 = l1 * l1 * l3
    l222 = l2 * l2 * l2
    l221 = l2 * l2 * l1
    l223 = l2 * l2 * l3
    l333 = l3 * l3 * l3
    l331 = l3 * l3 * l1
    l332 = l3 * l3 * l2
    l123 = l1 * l2 * l3
    α, β, γ = zeros(3), zeros(3), zeros(3)
    α[1] = l111 + 3 * l112 + 3 * l113 + 2 * l123
    α[2] = l222 + 3 * l221 + 3 * l223 + 2 * l123
    α[3] = l333 + 3 * l331 + 3 * l332 + 2 * l123
    β[1] = (x2 - x1) * (l112 + 0.5*l123) + (x3 - x1) * (l113 + 0.5*l123)
    β[2] = (x1 - x2) * (l221 + 0.5*l123) + (x3 - x2) * (l223 + 0.5*l123)
    β[3] = (x1 - x3) * (l331 + 0.5*l123) + (x2 - x3) * (l332 + 0.5*l123)
    γ[1] = (y2 - y1) * (l112 + 0.5*l123) + (y3 - y1) * (l113 + 0.5*l123)
    γ[2] = (y1 - y2) * (l221 + 0.5*l123) + (y3 - y2) * (l223 + 0.5*l123)
    γ[3] = (y1 - y3) * (l331 + 0.5*l123) + (y2 - y3) * (l332 + 0.5*l123)

    return α, β, γ
end

function bchs_grad_coefficients(l1, l2, l3, triangle)
    p, q, r = triangle
    x1, y1 = p
    x2, y2 = q
    x3, y3 = r
    S1 = x2*y3 - x3*y2
    S2 = x3*y1 - x1*y3
    S3 = x1*y2 - x2*y1
    S = S1 + S2 + S3
    l1x = (y2 - y3) / S
    l2x = (y3 - y1) / S
    l3x = (y1 - y2) / S
    l1y = (x3 - x2) / S
    l2y = (x1 - x3) / S
    l3y = (x2 - x1) / S
    αx, βx, γx = zeros(3), zeros(3), zeros(3)
    αy, βy, γy = zeros(3), zeros(3), zeros(3)
    a1 = 3 * l1 * l1 + 6 * l1 * l2 + 6 * l1 * l3 + 2 * l2 * l3
    a2 = 3 * l1 * l1 + 2 * l1 * l3
    a3 = 3 * l1 * l1 + 2 * l1 * l2
    αx[1] = a1 * l1x + a2 * l2x + a3 * l3x
    αy[1] = a1 * l1y + a2 * l2y + a3 * l3y
    a1 = 3 * l2 * l2 + 2 * l2 * l3
    a2 = 3 * l2 * l2 + 6 * l1 * l2 + 6 * l2 * l3 + 2 * l1 * l3
    a3 = 3 * l2 * l2 + 2 * l1 * l2
    αx[2] = a1 * l1x + a2 * l2x + a3 * l3x
    αy[2] = a1 * l1y + a2 * l2y + a3 * l3y
    a1 = 3 * l3 * l3 + 2 * l2 * l3
    a2 = 3 * l3 * l3 + 2 * l1 * l3
    a3 = 3 * l3 * l3 + 6 * l1 * l3 + 6 * l2 * l3 + 2 * l1 * l2
    αx[3] = a1 * l1x + a2 * l2x + a3 * l3x
    αy[3] = a1 * l1y + a2 * l2y + a3 * l3y
    b1 = 2 * l1 * l2 + 0.5 * l2 * l3
    b2 = l1 * l1 + 0.5 * l1 * l3
    b3 = 0.5 * l1 * l2
    c1 = 2 * l1 * l3 + 0.5 * l2 * l3
    c2 = 0.5 * l1 * l3
    c3 = l1 * l1 + 0.5 * l1 * l2
    βx[1] = (x2 - x1) * (b1 * l1x + b2 * l2x + b3 * l3x) + (x3 - x1) * (c1 * l1x + c2 * l2x + c3 * l3x)
    βy[1] = (x2 - x1) * (b1 * l1y + b2 * l2y + b3 * l3y) + (x3 - x1) * (c1 * l1y + c2 * l2y + c3 * l3y)
    γx[1] = (y2 - y1) * (b1 * l1x + b2 * l2x + b3 * l3x) + (y3 - y1) * (c1 * l1x + c2 * l2x + c3 * l3x)
    γy[1] = (y2 - y1) * (b1 * l1y + b2 * l2y + b3 * l3y) + (y3 - y1) * (c1 * l1y + c2 * l2y + c3 * l3y)
    b1 = l2 * l2 + 0.5 * l2 * l3
    b2 = 2 * l1 * l2 + 0.5 * l1 * l3
    b3 = 0.5 * l1 * l2
    c1 = 0.5 * l2 * l3
    c2 = 2 * l2 * l3 + 0.5 * l1 * l3
    c3 = l2 * l2 + 0.5 * l1 * l2
    βx[2] = (x1 - x2) * (b1 * l1x + b2 * l2x + b3 * l3x) + (x3 - x2) * (c1 * l1x + c2 * l2x + c3 * l3x)
    βy[2] = (x1 - x2) * (b1 * l1y + b2 * l2y + b3 * l3y) + (x3 - x2) * (c1 * l1y + c2 * l2y + c3 * l3y)
    γx[2] = (y1 - y2) * (b1 * l1x + b2 * l2x + b3 * l3x) + (y3 - y2) * (c1 * l1x + c2 * l2x + c3 * l3x)
    γy[2] = (y1 - y2) * (b1 * l1y + b2 * l2y + b3 * l3y) + (y3 - y2) * (c1 * l1y + c2 * l2y + c3 * l3y)
    b1 = l3 * l3 + 0.5 * l2 * l3
    b2 = 0.5 * l1 * l3
    b3 = 2 * l1 * l3 + 0.5 * l1 * l2
    c1 = 0.5 * l2 * l3
    c2 = l3 * l3 + 0.5 * l1 * l3
    c3 = 2 * l2 * l3 + 0.5 * l1 * l2
    βx[3] = (x1 - x3) * (b1 * l1x + b2 * l2x + b3 * l3x) + (x2 - x3) * (c1 * l1x + c2 * l2x + c3 * l3x)
    βy[3] = (x1 - x3) * (b1 * l1y + b2 * l2y + b3 * l3y) + (x2 - x3) * (c1 * l1y + c2 * l2y + c3 * l3y)
    γx[3] = (y1 - y3) * (b1 * l1x + b2 * l2x + b3 * l3x) + (y2 - y3) * (c1 * l1x + c2 * l2x + c3 * l3x)
    γy[3] = (y1 - y3) * (b1 * l1y + b2 * l2y + b3 * l3y) + (y2 - y3) * (c1 * l1y + c2 * l2y + c3 * l3y)

    return αx, βx, γx, αy, βy, γy
end

function _interp(spl::BivariateCHSInterpolation, x, y)
    triangles = spl.triangles
    tri = find_triangle(triangles, [x, y])
    any(tri .< 0) && return NaN

    i, j, k = tri
    tri_points = get_point(triangles, i, j, k)
    l1, l2, l3 = cartesian2baricentric(tri_points, x, y)
    z = (spl.z[i], spl.z[j], spl.z[k])
    dzdx = (spl.dzdx[i], spl.dzdx[j], spl.dzdx[k])
    dzdy = (spl.dzdy[i], spl.dzdy[j], spl.dzdy[k])
    f = zero(z[1])
    α, β, γ = bchs_coefficients(l1, l2, l3, tri_points)
    for m in 1:3
        f += α[m] * z[m] + β[m] * dzdx[m] + γ[m] * dzdy[m]
    end

    return f
end

function _grad(spl::BivariateCHSInterpolation, x, y)
    triangles = spl.triangles
    tri = find_triangle(triangles, [x, y])
    any(tri .< 0) && return (NaN, NaN)

    i, j, k = tri
    tri_points = get_point(triangles, i, j, k)
    l1, l2, l3 = cartesian2baricentric(tri_points, x, y)
    z = (spl.z[i], spl.z[j], spl.z[k])
    dzdx = (spl.dzdx[i], spl.dzdx[j], spl.dzdx[k])
    dzdy = (spl.dzdy[i], spl.dzdy[j], spl.dzdy[k])
    fx, fy = zero(z[1]), zero(z[1])
    αx, βx, γx, αy, βy, γy = bchs_grad_coefficients(l1, l2, l3, tri_points)
    for m in 1:3
        fx += αx[m] * z[m] + βx[m] * dzdx[m] + γx[m] * dzdy[m]
        fy += αy[m] * z[m] + βy[m] * dzdx[m] + γy[m] * dzdy[m]
    end

    return fx, fy
end

(spl::BivariateCHSInterpolation)(x::Real, y::Real) = _interp(spl, x, y)
(spl::BivariateCHSInterpolation)(x::AbstractVector, y::AbstractVector) = spl.(x, y)

# handy methods
interp(spl::BivariateCHSInterpolation, x, y) = spl(x, y)
interp(spl::BivariateCHSInterpolation, p) = spl(p[1], p[2])
grad(spl::BivariateCHSInterpolation, x, y) = _grad(spl, x, y)
grad(spl::BivariateCHSInterpolation, p) = _grad(spl, p[1], p[2])
