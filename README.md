# CubicHermiteSpline.jl

**CubicHermiteSpline.jl** is a naive implementation of cubic Hermite spline interpolation for 1D data points in pure Julia. Currently, the 1st order gradient should be given by the user. It is most useful when the gradient happens to be available. When the function to be interpolated is smooth and the accuracy of the gradients is high, the cubic Hermite spline interpolation should perform extremely well. A demonstration of the power of this interpolation can be found [here](http://www.yxliu.group/2020/06/cubic-hermite-spline).

## Features

* Univariate cubic Hermite spline interpolation for 1D data points (regular and irregular grids are both supported).
* Gradient (1st order derivative) of the interpolation. (New in version 0.2.0)
* Bivariate cubic Hermite spline interpolation for 2D data points (regular and irregular grids are both supported). (New in version 0.3.0)

## Usage

### 1D Interpolation

Below shows a trivial example when the function is a cubic polynomial. In such case, the function can be exactly interpolated.

First, prepare a set of data points to be interpolated. Note that here we use a cubic polynomial function which can be exactly interpolated by the cubic Hermite spline method.

```julia
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
```

The gradients at each data points are also computed which is required by the cubic Hermite spline method.

```julia
julia> gradient = df.(x)
6-element Array{Float64,1}:
  2.0
 -0.25
 -1.0
 -0.25
  2.0
  5.75
```

Then, we construct a interpolation instance by using CubicHermiteSpline package.

```julia
julia> spl = UnivariateCHSInterpolation(x, y, gradient);
```

Perform interpolation for a single input x.

```julia
julia> xi = 1.2;

julia> yi = spl(xi)  # Or using interp(spl, xi)
-5.192
```

Perform interpolation for an array of input x.

```julia
julia> xi = [0.5, 1.2];

julia> yi = spl(xi)
2-element Array{Float64,1}:
 -4.625
 -5.192
```

The 1st order derivative of the interpolation can be obtained.

```julia
julia> xi = 1.2;

julia> ki = spl(xi; grad=true)  # Or using grad(spl, xi)
```

Note that 1st order derivatives at each data point should be provided by the user. Please see `doc/tutorial_univariate.ipynb` for a detailed example of the univariate cubic Hermite spline interpolation.

### 2D Interpolation

Construct a 2D interpolation instance:

```julia
spl2d = BivariateCHSInterpolation(rand(100), rand(100), rand(100), rand(100), rand(100))
```

Perform interpolation for a single input point (x, y):

```julia
spl2d(x, y)
# Or use
interp(spl2d, x, y)
```

Compute the interpolated 1st order derivatives at a single input point (x, y):

```julia
grad(spl2d, x, y)
```

Please see `doc/tutorial_bivariate.ipynb` for a concrete example of bivariate cubic Hermite spline interpolation.

## Contribute

* Star the package on [github.com](https://github.com/liuyxpp/CubicHermiteSpline.jl).
* File an issue or make a pull request on [github.com](https://github.com/liuyxpp/CubicHermiteSpline.jl).
* Contact the author via email <lyx@fudan.edu.cn>.

## Links

* [Source code](https://github.com/liuyxpp/CubicHermiteSpline.jl)
* [Documentation](http://www.yxliu.group/2020/06/cubic-hermite-spline)
* [Tutorial in IPython Notebook](https://github.com/liuyxpp/CubicHermiteSpline.jl/blob/master/doc/tutorial.ipynb)
