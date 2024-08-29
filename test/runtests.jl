using CubicHermiteSpline
using Test

@testset "basis function" begin
    h00, h10, h01, h11 = CubicHermiteSpline.basis(0)
    @test h00 == 1
    @test h10 == 0
    @test h01 == 0
    @test h11 == 0

    h00, h10, h01, h11 = CubicHermiteSpline.basis(1)
    @test h00 == 0
    @test h10 == 0
    @test h01 == 1
    @test h11 == 0

    t = 0.5
    h00, h10, h01, h11 = CubicHermiteSpline.basis(t)
    @test h00 ≈ (1 + 2t) * (t - 1)^2
    @test h10 ≈ t * (t - 1)^2
    @test h01 ≈ t^2 * (3 - 2t)
    @test h11 ≈ t^2 * (t - 1)
end

@testset "basis_derivative function" begin
    h00, h10, h01, h11 = CubicHermiteSpline.basis_derivative(0)
    @test h00 == 0
    @test h10 == 1
    @test h01 == 0
    @test h11 == 0

    h00, h10, h01, h11 = CubicHermiteSpline.basis_derivative(1)
    @test h00 == 0
    @test h10 == 0
    @test h01 == 0
    @test h11 == 1

    t = 0.5
    h00, h10, h01, h11 = CubicHermiteSpline.basis_derivative(t)
    @test h00 ≈ 6t * (t - 1)
    @test h10 ≈ (3t - 1) * (t - 1)
    @test h01 ≈ -6t * (t - 1)
    @test h11 ≈ t * (3t - 2)
end

@testset "Interval enclosing the data point" begin
    x = [0, 3, 5, 8, 10, 100]
    idx, x1, x2 = CubicHermiteSpline.findinterval(5, x)
    @test idx == 3
    @test x1 == 5
    @test x2 === nothing
    idx, x1, x2 = CubicHermiteSpline.findinterval(11, x)
    @test idx == 5
    @test x1 == 10
    @test x2 == 100

    x = [0.0, 0.3, 0.5, 0.8, 10.5, 100.1]
    idx, x1, x2 = CubicHermiteSpline.findinterval(0.5, x)
    @test idx == 3
    @test x1 == 0.5
    @test x2 === nothing
    idx, x1, x2 = CubicHermiteSpline.findinterval(0.99, x)
    @test idx == 4
    @test x1 == 0.8
    @test x2 == 10.5
end

@testset "Interpolation" begin
    # Cubic Hermite spline interpolation should exactly interpolate the 3rd order polynomial.
    f(x) = x^3 - 3x^2 + 2x - 5
    df(x) = 3x^2 - 6x + 2
    x = range(0, 2.5, step=0.5)
    y = f.(x)
    dydx = df.(x)
    spl = CubicHermiteSplineInterpolation(x, y, dydx)

    # Interpolate at the input data point, simply return input data
    xi = 0.5
    @test spl(xi) ≈ f(xi)
    @test interp(spl, xi) ≈ f(xi)
    # Interpolate at arbitrary position
    xi = 1.2
    @test spl(xi) ≈ f(xi)
    @test interp(spl, xi) ≈ f(xi)

    # Test array interpolation
    xi = [0.5, 1.2]
    yi = spl(xi)
    @test yi[1] ≈ f(xi[1])
    @test yi[2] ≈ f(xi[2])
    yi = interp(spl, xi)
    @test yi[1] ≈ f(xi[1])
    @test yi[2] ≈ f(xi[2])
end

@testset "Interpolation grad" begin
    # Cubic Hermite spline interpolation should exactly interpolate the 3rd order polynomial.
    f(x) = x^3 - 3x^2 + 2x - 5
    df(x) = 3x^2 - 6x + 2
    # x = range(0, 2.5, step=0.5)
    x = [0, 0.3, 0.8, 1.5, 2.7, 3]
    y = f.(x)
    dydx = df.(x)
    spl = CubicHermiteSplineInterpolation(x, y, dydx)

    # Interpolate at the input data point, simply return input data
    xi = 0.5
    @test spl(xi; grad=true) ≈ df(xi)
    @test grad(spl, xi) ≈ df(xi)
    # Interpolate at arbitrary position
    xi = 2.1
    @test spl(xi; grad=true) ≈ df(xi)
    @test grad(spl, xi) ≈ df(xi)

    # Test array interpolation
    xi = [0.5, 2.1]
    yi = spl(xi; grad=true)
    @test yi[1] ≈ df(xi[1])
    @test yi[2] ≈ df(xi[2])
    yi = grad(spl, xi)
    @test yi[1] ≈ df(xi[1])
    @test yi[2] ≈ df(xi[2])
end