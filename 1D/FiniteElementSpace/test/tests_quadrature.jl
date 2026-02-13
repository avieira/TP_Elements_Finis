function test_quadrature1()
    cell = ReferenceInterval(0.0,1.0)
    quadrature = gauss_legendre_quadrature(cell, 2)
    result = integrate(quadrature, x->x^2)
    return isapprox(result, 1.0/3.0, atol=1e-6)
end

function test_quadrature2()
    cell = ReferenceInterval()
    quadrature = gauss_legendre_quadrature(cell, 2)
    result = integrate(quadrature, x->x^3)
    return isapprox(result, 0.0, atol=1e-6)
end

function test_quadrature3()
    cell = ReferenceInterval()
    quadrature = gauss_legendre_quadrature(cell, 1)
    return (quadrature.degree == 1) && (length(quadrature.points) == 1) && (length(quadrature.weights) == 1) &&
             (isapprox(quadrature.points[1], 0.0)) && (isapprox(quadrature.weights[1], 2.0))
end

@testset "Tests Quadrature" begin
    @test test_quadrature1()
    @test test_quadrature2()
    @test test_quadrature3()
end