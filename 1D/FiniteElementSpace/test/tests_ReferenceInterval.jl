function test_ReferenceInterval1()
    ref_interval = ReferenceInterval()
    test1 = ref_interval.a == -1.0
    test2 = ref_interval.b == 1.0
    return test1 && test2
end

function test_ReferenceInterval2()
    ref_interval = ReferenceInterval(-2.0, 3.0)
    test1 = ref_interval.a == -2.0
    test2 = ref_interval.b == 3.0
    return test1 && test2
end


@testset "Tests ReferenceInterval" begin
    @test test_ReferenceInterval1()
    @test test_ReferenceInterval2()
end