using SparseArrays

"""
    boundary_nodes(fes::FiniteElementSpace1D)

    Returns the list of nodes that are on the boundary of fes.mesh.
"""
function boundary_nodes(fes::FiniteElementSpace1D)
    a = minimum(fes.mesh.vertices)
    b = maximum(fes.mesh.vertices)
    
    on_boundary = x->begin
        if abs(x-a)<eps() || abs(x-b)<eps()
            return 1.
        else
            return 0.
        end
    end
    fct = FunctionSpace1D(fes, on_boundary)
    boundary_node_indices = findall(x->x!=0.0, fct.coefficients)
    return boundary_node_indices
end

"""
assemble_rhs(fes::FiniteElementSpace1D, f::FunctionSpace1D)
Assembles the right-hand side vector for the 1D Helmholtz equation.

Returns: Vector{Float64}
"""
function assemble_rhs(fes::FiniteElementSpace1D, f::FunctionSpace1D)
    rhs = zeros(Float64, fes.node_count)
    quad = gauss_legendre_quadrature(fes.fe.cell, 2*fes.fe.degree)
    phi = evaluate_basis(fes.fe, quad.points)
    
    nb_node_per_elem = size(phi, 2) # = fes.fe.degree + 1
    
    for c in 1:fes.mesh.entity_counts[end] #For each cell
        for i in 1:nb_node_per_elem
            N = fes.cell_node_mappings[c, i]
            coeff_f_in_c = f.coefficients[fes.cell_node_mappings[c, :]]
            #TODO : complete here
        end
    end

    return rhs
end

"""
assemble_lhs(fes::FiniteElementSpace1D, f::FunctionSpace1D)
Assembles the left-hand side matrix for the 1D Helmholtz equation.

Returns: SparseMatrixCSC{Int, Float64}
"""
function assemble_lhs(fes::FiniteElementSpace1D)
    lhs = spzeros(Float64, fes.node_count, fes.node_count)
    quad = gauss_legendre_quadrature(fes.fe.cell, 2*fes.fe.degree)
    phi = evaluate_basis(fes.fe, quad.points)
    dphi = evaluate_basis(fes.fe, quad.points, true)

    nb_node_per_elem = size(phi, 2) # = fes.fe.degree + 1

    #TODO

    for c in 1:fes.mesh.entity_counts[end]
        for i in 1:nb_node_per_elem
            global_i = fes.cell_node_mappings[c, i]
            for j in 1:nb_node_per_elem
                global_j = fes.cell_node_mappings[c, j]
                #TODO
            end
        end
    end

    return lhs
end

"""
assemble_lhs_rhs(fes::FiniteElementSpace1D, f::FunctionSpace1D)
Assembles the left-hand side matrix and the right-hand side vector for the 1D Helmholtz equation.

Returns: SparseMatrixCSC{Int, Float64}, Vector{Float64}
"""
function assemble_lhs_rhs(fes::FiniteElementSpace1D, f::FunctionSpace1D)
    return assemble_lhs(fes), assemble_rhs(fes, f)
end

"""
apply_dirichlet_elimination(fes::FiniteElementSpace1D, lhs::SparseMatrixCSC{Float64, Int}, rhs::Vector{Float64})
Modifies lhs and rhs to implement Dirichlet boundary conditions using line/column elimination.
"""
function apply_dirichlet_elimination(fes::FiniteElementSpace1D, lhs::SparseMatrixCSC{Float64, Int}, rhs::Vector{Float64})
    #TODO
end

"""
apply_dirichlet_penalization(fes::FiniteElementSpace1D, lhs::SparseMatrixCSC{Float64, Int}, penal::Real = 1e20)
Modifies lhs to implement Dirichlet boundary conditions using diagonal penalization.
"""
function apply_dirichlet_penalization(fes::FiniteElementSpace1D, lhs::SparseMatrixCSC{Float64, Int}, penal::Real = 1e20)
    #TODO
end

"""
solve_helmholtz(fes::FiniteElementSpace1D, f::FunctionSpace1D, dirichlet_bc::Bool=false)
Solves Helmholtz equation in 1D, with right-hand-side function `f`.
If dirichlet_bc=false, homogeneous Neumann  boundary conditions are used.
Otherwise, homogeneous Dirichlet boundary conditions are used.
"""
function solve_helmholtz(fes::FiniteElementSpace1D, f::FunctionSpace1D, dirichlet_bc::Bool=false)
    lhs, rhs = assemble_lhs_rhs(fes, f)

    if dirichlet_bc
        apply_dirichlet_elimination(fes, lhs, rhs)
        #apply_dirichlet_penalization(fes, lhs)
    end

    # Solve the linear system
    coeffs = lhs \ rhs
    return FunctionSpace1D(fes, coeffs)
end