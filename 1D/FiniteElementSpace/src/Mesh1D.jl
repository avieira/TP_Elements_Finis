using SparseArrays
import Gmsh

export Mesh1D, generate_mesh_from_successive_points, generate_mesh_structured
"""
    Mesh2D

    Structure representing a mesh in a 1 dimensional domain.
    It consists of 3 arrays representing the vertices and the edges of the mesh. It supposes that all segments are straight portions of lines, no curve is represented.
    Other arrays reprensent some connectivity in the mesh, staticstics on it and a precomputed jacobian.

    # Arguments
    - vertices: Vector{Point1D} - List of all vertices composing the mesh.
    - edges: SparseMatrixCSC - Matrix representing the edges, of size (nb of vertices x nb of edges).
                                A column represents an edge between two points.
                                For instance, if edge i links points a and b, the column i of this array will be all 0 except in line a where it is -1 and in line b where it is 1.
    - entity_counts: Vector{Int} - Vector of size 2 giving the number of vertices and edges.
    - gradient: Matrix{Float64} - Matrix used for jacobian computation.
"""
struct Mesh1D
    vertices::Vector{Point1D}
    edges::SparseMatrixCSC{Int, Int}
    entity_counts::Vector{Int}
    scaling::Float64
end

#Constructor
function Mesh1D(vertices::Vector{Point1D}, edges::SparseMatrixCSC{Int, Int}, entity_counts::Vector{Int}, ref_cell::ReferenceInterval)
    #TODO : change scaling
    scaling = 1.0
    return Mesh1D(vertices, edges, entity_counts, scaling)
end

"""
    generate_mesh_from_successive_points(vertices::Vector{Point1D})

    Generate a [`Mesh1D`](@ref) given the points in vertices.
    This function supposes that all points in vertices are ordered in ascending order, such that the subintervals are simply [vertices[i], vertices[i+1]].
"""
function generate_mesh_from_successive_points(vertices::Vector{Point1D}, ref_cell::ReferenceInterval = ReferenceInterval())
    edges = spzeros(Int, length(vertices), length(vertices)-1)
    for i in 1:length(vertices)-1
        edges[i, i] = -1
        edges[i+1, i] = 1
    end
    entity_counts = [length(vertices), size(edges, 2)]
    return Mesh1D(vertices, edges, entity_counts, ref_cell)
end

"""
    generate_mesh_structured(a::<:Real, b::<:Real, num_elements::Int, ref_cell::ReferenceInterval = ReferenceInterval())

    Generate a [`Mesh1D`](@ref) for an interval [a,b]. All subintervals will have the same length (b-a)/n_elements.
    This interval will be discretized using num_elements+1 points and num_elements subintervals.
"""
function generate_mesh_structured(a::T, b::T, num_elements::Int, ref_cell::ReferenceInterval = ReferenceInterval()) where {T<:Real}
    #TODO
    return Mesh1D(vertices, edges, entity_counts, ref_cell)
end

export adjacency
"""
    adjacency(mesh::Mesh1D, d1::Int, d2::Int, i::Int)

    Returns the list of entity indices of dimension d2 connected to the i-th entity of dimension d1.
    For exemple, adjacency(mesh, 1, 0, i) returns the vertices connected to the i-th edge.
                 adjacency(mesh, 0, 1, i) returns the edges connected to the i-th node.
"""
function adjacency(mesh::Mesh1D, d1::Int, d2::Int, i::Int)
    if d1 == 0 && d2 == 1
        # Node to Edge
        connected_edges = findnz(mesh.edges[i, :])[1]
        return connected_edges
    elseif d1 == 1 && d2 == 0
        # Edge to Node
        connected_vertices = rowvals(mesh.edges)[nzrange(mesh.edges, i)]
        return connected_vertices
    elseif d1 == 1 && d2 == 1
        # Edge to Edge 
        return [i]
    else
        error("Adjacency not defined for the given dimensions.")
    end
end

export jacobian
"""
    jacobian(cell::Mesh2D, i_cell::Int)
    Returns the derivative of the linear transform needed for the integration by substitution on the reference cell.
    
    # Arguments
    - cell::Mesh1D - a mesh
    - i_cell::Int - index of the cell where the integration occurs.
"""
function jacobian(mesh::Mesh1D, i_cell::Int)
    #TODO
    return 0
end