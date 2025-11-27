module ChebyshevLobattoQuadrature

export cheb_lobatto_nodes, cheb_lobatto_weights, cheb_lobatto_quadrature

@inline function _checked_order(n::Integer)
    n >= 1 || throw(ArgumentError("n must be a positive integer, got $n"))
    return Int(n)
end

"""
    cheb_lobatto_nodes(n::Integer) -> Vector{Float64}

Return the Chebyshev–Lobatto nodes ``x_k = \\cos(k\\pi/n)`` for ``k = 0,\\ldots,n``.

# Examples
```julia
julia> cheb_lobatto_nodes(4)
5-element Vector{Float64}:
  1.0
  0.7071067811865476
  6.123233995736766e-17
 -0.7071067811865475
 -1.0
```
"""
@inline function cheb_lobatto_nodes(n::Integer)
    N = _checked_order(n)
    nodes = Vector{Float64}(undef, N + 1)
    invN = 1 / N
    @inbounds for k in 0:N
        nodes[k + 1] = cospi(k * invN)
    end
    return nodes
end

"""
    cheb_lobatto_weights(n::Integer) -> Vector{Float64}

Compute Chebyshev–Lobatto weights for integrals on ``[-1, 1]`` with respect to
the standard Lebesgue measure. The weights solve the moment conditions
``\\int_{-1}^{1} x^i\\,dx = \\sum_{k=0}^n w_k x_k^i`` for ``i = 0,\\ldots,n``.
"""
@inline function cheb_lobatto_weights(n::Integer)
    N = _checked_order(n)
    nodes = cheb_lobatto_nodes(N)
    A = Matrix{Float64}(undef, N + 1, N + 1)
    @inbounds for i in 0:N
        A[i + 1, :] = nodes .^ i
    end
    b = Vector{Float64}(undef, N + 1)
    @inbounds for i in 0:N
        b[i + 1] = iseven(i) ? 2.0 / (i + 1) : 0.0
    end
    return A \ b
end

"""
    cheb_lobatto_quadrature(f::Function, n::Integer) -> Float64

Approximate the integral of ``f`` on ``[-1, 1]`` using Chebyshev–Lobatto nodes
and weights obtained from the moment conditions:

``I \\approx \\sum_{k=0}^n w_k f(x_k)``.

# Examples
```julia
julia> cheb_lobatto_quadrature(x -> 1, 32)
2.0
```
"""
@inline function cheb_lobatto_quadrature(f::Function, n::Integer)
    nodes = cheb_lobatto_nodes(n)
    weights = cheb_lobatto_weights(n)
    values = f.(nodes)
    return sum(weights .* values)
end

end # module ChebyshevLobattoQuadrature
