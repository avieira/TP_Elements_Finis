using Pkg
Pkg.activate(".")

using FiniteElementSpace
using LinearAlgebra, SparseArrays

let
    umesh = generate_mesh_structured(0.0, 1.0, 20)
    fe = LagrangeElement1D(ReferenceInterval(), 2)
    fes = FiniteElementSpace1D(umesh, fe)
    fct = FunctionSpace1D(fes, x->(2-x)*(x+1)) #solution is x*(1-x)
    u = solve_helmholtz(fes, fct, true) 

    analytic = FunctionSpace1D(fes, x->x*(1-x))
    FiniteElementSpace.plot([u, analytic]) 
end