module FiniteElementSpace

include("ReferenceInterval.jl")

include("FiniteElement.jl")

include("Quadrature.jl")

include("Mesh1D.jl")

using Plots

export FiniteElementSpace1D
"""
    FiniteElementSpace1D

    Structure defining a finite element space on a 1D mesh.

    # Arguments
    - mesh::Mesh1D
    - fe::FiniteElement1D - Finite element defined on a [`ReferenceInterval`](@ref)
    - cell_node_mappings::Array{Int, 2} - Global numbering linking the local numbering of nodes on the [`ReferenceIntevral`](@ref) to the global numbering of nodes on the mesh.
    - node_count::Int - Total number of nodes in the mesh.
"""
struct FiniteElementSpace1D
    mesh::Mesh1D
    fe::FiniteElement1D
    cell_node_mappings::Array{Int, 2}
    node_count::Int
end

#Constructor
function FiniteElementSpace1D(mesh::Mesh1D, fe::FiniteElement1D)
    dim = 1
    cell_node_mappings = Array{Int, 2}(undef, mesh.entity_counts[dim+1], length(fe.nodes))
    #TODO
    
    node_count = 0
    return FiniteElementSpace1D(mesh, fe, cell_node_mappings, node_count)
end

export plot_mesh_and_global_num
"""
    plot_mesh_and_global_num(mesh::Mesh1D, fes::FiniteElementSpace1D)

    Plot the edges of the mesh, and, if there are not too many, all the nodes in `fes` with their global numbering.
"""
function plot_mesh_and_global_num(mesh::Mesh1D, fes::FiniteElementSpace1D)
    p = Plots.plot()
    pty = [0., 0.]
    for e in axes(mesh.edges, 2)
        connected_vertices = rowvals(mesh.edges)[nzrange(mesh.edges, e)]
        ptx = [mesh.vertices[connected_vertices[1]].x, mesh.vertices[connected_vertices[2]].x]
        Plots.plot!(p, ptx, pty, linewidth=2, color = :black, legend=false)
    end

    cg1 = LagrangeElement2D(fes.fe.cell, 1)
    coord_map = evaluate_basis(cg1, fes.fe.nodes)
    cg1fs = FiniteElementSpace2D(fes.mesh, cg1)

    if fes.node_count>100
        println("Too many nodes to print ; I ain't print all that...")
    else
        for c in 1:fes.mesh.entity_counts[end]
            v_coords = fes.mesh.vertices[cg1fs.cell_node_mappings[c, :]]
            ptx = coord_map * v_coords
            pty = repeat([0.], length(ptx))
            Plots.scatter!(p, ptx, pty, markersize=4, color=:red)
            [annotate!(p, ptx, pty-0.5, Plots.text(string(fes.cell_node_mappings[c, i]), 12)) for (i, ptx, pty) in zip(eachindex(ptx),ptx,pty)];
        end
    end
    display(p)
end

export FunctionSpace1D, interpolate, evaluate, integrate
"""
    FunctionSpace1D
    Structure defining functions defined using finite elements.

    # Arguments
    - fes::FiniteElementSpace1D
    - coefficients::Vector{Float64} - coefficients for the function basis of fes, of size fes.node_count.

"""
struct FunctionSpace1D
    fes::FiniteElementSpace1D
    coefficients::Vector{Float64}
    FunctionSpace1D(fes::FiniteElementSpace1D, coefficients::Vector{Float64}) = begin
        if length(coefficients) != fes.node_count
            error("Coefficient vector length does not match the number of nodes in the finite element space.")
        end
        return new(fes, coefficients)
    end
    FunctionSpace1D(fes::FiniteElementSpace1D) = new(fes, zeros(fes.node_count))
end

"""
    interpolate(fes::FiniteElementSpace1D, f::Function)

    Interpolation of a function `f` in the finite element space `fes`.
    This fonction returns a [`FunctionSpace1D`](@ref) based on `fes` and computes all the coefficients of the interpolation of `f` on this finite element space.
"""
function interpolate(fes::FiniteElementSpace1D, f::Function)
    fct = FunctionSpace1D(fes)

    #TODO
    #cg1 = LagrangeElement1D(fes.fe.cell, ???)
    #λ_array = evaluate_basis(cg1, ???)
    #cg1fs = FiniteElementSpace1D(fes.mesh, cg1)

    #for c in 1:fes.mesh.entity_counts[end] #for each cell
    #    v_coords = fes.mesh.vertices[???] #Hint: you should use cg1fs
    #    node_coords = ??? #Hint: it is a simple matrix multiplication
    #    fct.coefficients[???] = map(f, node_coords)
    #end

    return fct
end

function FunctionSpace1D(fes::FiniteElementSpace1D, f::Function)
    return interpolate(fes, f)
end

"""
    integrate(fs::FunctionSpace1D)

    Performs the integration of a function defined on a finite element space `fs` on the domain defined by `fs.fes.mesh`.
"""
function integrate(fs::FunctionSpace1D)
    total_integral = 0.0
    #TODO

    return total_integral
end

"""
    plot(fss::Vector{FunctionSpace2D})

    Plots all the functions in `fss`.
"""
function plot(fss::Vector{FunctionSpace1D}; num_points_per_element::Int=10)
    #plots = Vector{typeof(Plots.plot())}(undef, length(fss))
    xs = Vector{Vector{Float64}}(undef, length(fss))
    ys = Vector{Vector{Float64}}(undef, length(fss))
    labels = Vector{String}(undef, length(fss))
    for i in eachindex(fss)
        fs = fss[i]
        fes = fs.fes
        fe = fes.fe
        mesh = fes.mesh

        x_vals = Float64[]
        y_vals = Float64[]

        for c in 1:mesh.entity_counts[end]
            vertices_indices = adjacency(mesh, 1, 0, c)
            x_start = mesh.vertices[vertices_indices[1]]
            x_end = mesh.vertices[vertices_indices[2]]
            local_xs = range(x_start, x_end; length=num_points_per_element) |> collect
            phi = evaluate_basis(fe, (2*(local_xs .- x_start)/(x_end - x_start)) .- 1)
            local_coeffs = fs.coefficients[fes.cell_node_mappings[c, :]]
            local_ys = phi * local_coeffs

            append!(x_vals, local_xs)
            append!(y_vals, local_ys)
        end

        #plots[i] = Plots.plot(x_vals, y_vals, label="Function $i")
        xs[i] = x_vals
        ys[i] = y_vals
        labels[i] = "Function $i"
    end

    #Plots.plot(plots..., layout=(length(fss), 1), legend=true)
    #i = 1
    #Plots.plot(xs[i], ys[i]; label=labels[i])
    #for i in 2:length(fss)
    #    Plots.plot!(xs[i], ys[i]; label=labels[i])
    #end
    Plots.plot(xs, ys; label=labels, legend=true)
end

using LinearAlgebra, SparseArrays
export errornorm
"""
    errornorm(f1::FunctionSpace1D, f2::FunctionSpace1D)

    Compute the L^2 norm of the difference between f1 and f2.
"""
function errornorm(f1::FunctionSpace1D, f2::FunctionSpace1D)
    fes1 = f1.fes
    fes2 = f2.fes

    fe1 = fes1.fe
    fe2 = fes2.fe
    mesh = fes1.mesh

    # Create a quadrature rule which is exact for (f1-f2)**2.
    Q = gauss_legendre_quadrature(fe1.cell, 2*max(fe1.degree, fe2.degree))

    # Evaluate the local basis functions at the quadrature points.
    phi = evaluate_basis(fe1, Q.points)
    psi = evaluate_basis(fe2, Q.points)

    norm = 0.
    for c in 1:mesh.entity_counts[end] #for each cell in the mesh
        # Find the appropriate global node numbers for this cell.
        nodes1 = fes1.cell_node_mappings[c, :]
        nodes2 = fes2.cell_node_mappings[c, :]

        # Compute the change of coordinates.
        J = jacobian(mesh, c)

        # Compute the actual cell quadrature.
        norm += dot((phi*f1.coefficients[nodes1] -
                        psi*f2.coefficients[nodes2]).^2,
                       Q.weights) * J
    end

    return sqrt(norm)
end

include("Helmholtz.jl")

end # module FiniteElementSpace
