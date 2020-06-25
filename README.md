# CubicHermiteSpline.jl

**CubicHermiteSpline.jl** is a naive implementation of cubic Hermite spline interpolation for 1D data points in pure Julia. Currently, the 1st order gradient should be given by the user. It is most useful when the gradient happens to be available. When the function to be interpolated is smooth and the accuracy of the gradients is high, the cubic Hermite spline interpolation should perform extremely well. A demonstration of the power of this interpolation can be found [here](http://www.yxliu.group/2020/06/cubic-hermite-spline).

## Usage

Below shows a trivial example when the function is a cubic polynomial. In such case, the function can be exactly interpolated.

```console
julia> using CubicHermiteSpline

julia> f(x) = x^3 - 3x^2 + 2x - 5;

julia> df(x) = 3x^2 - 6x + 2;

julia> x = range(0, 2.5, step=0.5)
0.0:0.5:2.5

julia> y = f.(x)
6-element Array{Float64,1}:
 -5.0
 -4.625
 -5.0
 -5.375
 -5.0
 -3.125

julia> gradient = df.(x)
6-element Array{Float64,1}:
  2.0
 -0.25
 -1.0
 -0.25
  2.0
  5.75

julia> spl = CubicHermiteSplineInterpolation(x, y, gradient);

julia> xi = 1.2;

julia> yi = spl(xi)
-5.192

julia> xi = [0.5, 1.2];

julia> yi = spl(xi, yi)
2-element Array{Float64,1}:
 -4.625
 -5.192
```

Note that 1st order gradients at each data points should be provided by the user. Please see `doc/tutorial.ipynb` for a detailed demonstration of the usage of this package.

## TODO

* To support 2D and higher dimension data points.
* To allow computing gradients from data points when gradients are not provided by the user.

## Links

* [Source code](https://github.com/liuyxpp/CubicHermiteSpline.jl)
* [Documentation](http://www.yxliu.group/2020/06/cubic-hermite-spline)