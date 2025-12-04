module Collocation

using LinearAlgebra

using ..Generic: Generalized_Diff_Mat
using ..Chebyshev: ChebyNodes
using ..Legendre: Legendre_Lobatto_Basis

export SpectralGrid, build_grid

"""
    SpectralGrid

Container with nodes and differentiation matrices for a Lobatto-type grid.

- `x_ref`: nodes in the reference interval [-1, 1].
- `x`: nodes mapped to the physical `domain`.
- `D1`, `D2`: first- and second-order differentiation matrices in physical space.
- `basis`: `Symbol` identifying the family (`:chebyshev` or `:legendre`).
- `domain`: tuple `(a, b)` describing the physical interval.
"""
struct SpectralGrid{T}
    x_ref::Vector{T}
    x::Vector{T}
    D1::Matrix{T}
    D2::Matrix{T}
    basis::Symbol
    domain::Tuple{T, T}
end

"""
    build_grid(N; basis=:chebyshev, domain=(-1.0, 1.0))

Build a Lobatto grid with `N` collocation points.  The grid is returned as a
[`SpectralGrid`](@ref) and reused by the BVP/PDE solvers.
"""
function build_grid(N::Integer; basis::Symbol = :chebyshev, domain::Tuple = (-1.0, 1.0))
    N < 2 && error("The collocation grid requires at least two points.")
    a, b = domain
    a >= b && error("Invalid domain ($a, $b). Ensure a < b.")

    ξ = begin
        if basis === :chebyshev
            nodes, _ = ChebyNodes(N - 1, 2)
            collect(float.(nodes))
        elseif basis === :legendre
            _, _, nodes, _ = Legendre_Lobatto_Basis(N - 1)
            collect(float.(nodes))
        else
            error("Unsupported basis $basis. Use :chebyshev or :legendre.")
        end
    end

    scale = 2.0 / (b - a)
    x_phys = @. ((b - a) * ξ + (b + a)) / 2

    Dξ = Generalized_Diff_Mat(ξ)
    D1 = scale .* Dξ
    D2 = (scale^2) .* (Dξ * Dξ)

    return SpectralGrid(ξ, x_phys, D1, D2, basis, (float(a), float(b)))
end

end
