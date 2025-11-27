@testset "cheb_lobatto_nodes" begin
    for n in (4, 10)
        nodes = cheb_lobatto_nodes(n)
        @test length(nodes) == n + 1
        @test isapprox(first(nodes), 1.0; atol=1e-12)
        @test isapprox(last(nodes), -1.0; atol=1e-12)
        @test all(isapprox.(nodes, -reverse(nodes); atol=1e-12))
    end
    @test_throws ArgumentError cheb_lobatto_nodes(0)
end

@testset "cheb_lobatto_weights solves moments" begin
    n = 6
    weights = cheb_lobatto_weights(n)
    nodes = cheb_lobatto_nodes(n)
    for deg in 0:n
        lhs = sum(weights .* (nodes .^ deg))
        rhs = iseven(deg) ? 2.0 / (deg + 1) : 0.0
        @test isapprox(lhs, rhs; atol=1e-12, rtol=1e-12)
    end
    @test_throws ArgumentError cheb_lobatto_weights(0)
end

@testset "Cheb-Lobatto quadrature exactness on polynomials" begin
    f1(x) = 1.0
    f2(x) = x^2
    f3(x) = (2x^2 - 1)^2
    I1, I2, I3 = 2.0, 2.0 / 3.0, 14.0 / 15.0
    for n in (4, 6, 8, 10, 12)
        @test isapprox(cheb_lobatto_quadrature(f1, n), I1; atol=1e-12, rtol=1e-12)
        @test isapprox(cheb_lobatto_quadrature(f2, n), I2; atol=1e-12, rtol=1e-12)
        @test isapprox(cheb_lobatto_quadrature(f3, n), I3; atol=1e-12, rtol=1e-12)
    end
end

@testset "Cheb-Lobatto quadrature on sqrt(1-x^2)" begin
    f4(x) = sqrt(1 - x^2)
    I4 = pi / 2
    for n in (10, 20, 40, 80)
        @test isapprox(cheb_lobatto_quadrature(f4, n), I4; atol=1e-2, rtol=1e-2)
    end
end
