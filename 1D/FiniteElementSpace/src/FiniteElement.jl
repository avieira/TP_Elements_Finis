export FiniteElement1D, LagrangeElement1D, evaluate_basis, interpolate

"""
    VandermondeMatrix(nodes::Vector{Point2D}, degree::Int, grad::Bool=false)

    Computes the generalized Vandermonde matrix for the given `nodes` up to a certain `degree`.
    If `grad` is false, the output M will be a matrix of size (length(nodes) × (degree + 1)), 
                        where each column reprensents a monomial of a polynom.
    If `grad` is true, the output M will be a matrix of size (length(nodes) × (degree + 1)),
                        where M[:,:] will be the coefficients of the derivative of the polynom wrt to the variable `x`.
"""
function VandermondeMatrix(nodes::Vector{Point1D}, degree::Int, grad::Bool=false)
    V = Array{Float64}(undef, length(nodes), degree + 1)
    #TODO
    if !grad
    else 
    end
    return V
end

#############################################################
#              FINITE ELEMENT DEFINITION IN 1D              #
#############################################################
"""
    FiniteElement1D

    Structure representing a Finite Element on a reference interval.

    # Arguments
    - cell::ReferenceInterval
    - degree::Int
    - nodes::Vector{Point1D} - Nodes coordinates, whose number is defined by the degree.
    - basis_coeffs::Matrix{Float64} - Basis polyoms' coefficients ; namely, basis function ϕ_i(x) = basis_coeffs[i,:] * x
    - entity_nodes::Dict{Int,Dict{Int,Vector{Int}}} - Mapping from entity to local node indices ; 
                                                        entity_nodes[0] is a dictionnary of nodes associated to each vertex in the interval
                                                        entity_nodes[1] is a dictionnary of nodes associated to each edge in the interval (here, only 1 edge).
    - nb_nodes_per_entity::Vector{Int} - Number of nodes per entity (nb of nodes per vertex, nb of nodes per edge).
"""
mutable struct FiniteElement1D
    cell::ReferenceInterval
    degree::Int
    nodes::Vector{Point1D} # Nodes coordinates
    basis_coeffs::Matrix{Float64} # Basis function coefficients
    entity_nodes::Dict{Int, Dict{Int, Vector{Int}}} # Mapping from entity to local node indices
    nb_nodes_per_entity::Vector{Int} # Number of nodes per entity
end

"""
    evaluate_basis(fe::FiniteElement1D, x::Vector{Point1D}, grad::Bool=false)

    Evaluates the basis functions {hat{P}^i}_i in `fe` at the set of points `x`.
    If `grad = true`, evaluates the gradient of the basis function.
    It returns a matrix M of size (length(x) × fe.degree+1) where M_ij = hat{P}^j(x_i).
"""
function evaluate_basis(fe::FiniteElement1D, x::Vector{Point1D}, grad::Bool=false)
    #TODO
    return Matrix{Float64}(undef, 1, 1)
end

"""
    interpolate(fe::FiniteElement1D, f::Function)

    Returns a Vector{Float64} containing the coefficients of the interpolation of `f` based on the nodes defined in the reference finite element `fe`.
"""
function interpolate(fe::FiniteElement1D, f::Function)
    #TODO
    return Vector{Float64}(undef, 1)
end

#############################################################
#     SPECIALIZATION LAGRANGE ELEMENT DEFINITION IN 1D      #
#############################################################
"""
    LagrangeElement1D(cell::ReferenceInterval, degree::Int)

    Returns a [`FiniteElement1D`](@ref) on the [`ReferenceInterval`](@ref) of the given `degree` with Lagrange nodes and polynoms.
"""
function LagrangeElement1D(cell::ReferenceInterval, degree::Int)
    # Generate nodes 
    nodes = Vector{Point1D}(undef, 1) #TODO

    entity_nodes = Dict{Int, Dict{Int, Vector{Int}}}()
    entity_nodes[0] = Dict([1=>[]])  # vertices #TODO
    entity_nodes[1] = Dict([1 => []])  # edges #TODO

    nb_nodes_per_entity = length.([entity_nodes[d][1] for d in 0:1])

    # Generate basis coefficients (Vandermonde matrix)
    basis_coeffs = Matrix{Float64}(undef, 1, 1) #TODO

    return FiniteElement1D(cell, degree, nodes, basis_coeffs, entity_nodes, nb_nodes_per_entity)
end