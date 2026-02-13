export QuadratureRule1D, integrate, gauss_legendre_quadrature
using GaussQuadrature

"""
    QuadratureRule1D

    Defines a quadrature rule exact for a certain polynomial degree.

    # Arguments
    - cell::ReferenceInterval - Reference cell used for the quadrature.
    - points::Vector{Point1D} - Points used for the quadrature rule.
    - weights::Vector{Float64} - weights in the function evalutation used for the integration
    - degree::Int - degree of the quadrature.
"""
struct QuadratureRule1D
    cell::ReferenceInterval
    points::Vector{Point1D}
    weights::Vector{Float64}
    degree::Int
end

"""
    integrate(quadrature::QuadratureRule1D, f::Function)

    Approximately integrate `f` on the [`ReferenceInterval`](@ref) using the [`QuadratureRule1D`](@ref) quadrature.
"""
function integrate(quadrature::QuadratureRule1D, f::Function)
    return sum(quadrature.weights .* map(f, quadrature.points))
end

"""
    gauss_legendre_quadrature(cell::ReferenceIntevral, degree::Int)

    Generates a [`QuadratureRule1D`](@ref) of a given `degree` on the [`ReferenceInterval`](@ref).
"""
function gauss_legendre_quadrature(cell::ReferenceInterval, degree::Int)
    # Number of points needed for the given degree
    n_points = ceil(Int, (degree + 1) / 2)

    # Get the Gauss-Legendre points and weights on the reference interval [-1, 1]
    points_ref, weights_ref = legendre(n_points)

    # Map points from [-1, 1] to [a, b]
    a, b = cell.a, cell.b
    points = [(b - a) / 2 * p + (a + b) / 2 for p in points_ref]
    weights = [(b - a) / 2 * w for w in weights_ref]

    return QuadratureRule1D(cell, points, weights, degree)
end