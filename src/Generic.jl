module Generic

using LinearAlgebra

export Bary_Interp, Generalized_Diff_Mat, Poly_Roots

# -- Internal utilities ----------------------------------------------------

"""
    barycentric_weights(xs)

Return barycentric weights for the node vector `xs`.
"""
function barycentric_weights(xs::AbstractVector)
    x = vec(float.(xs))
    N = length(x)
    N >= 2 || error("At least two nodes are required to form weights.")
    w = ones(eltype(x), N)
    @inbounds for j in 1:N
        for k in 1:N
            k == j && continue
            w[j] *= (x[j] - x[k])
        end
    end
    return 1.0 ./ w
end

# -- Public API ------------------------------------------------------------

"""
    Bary_Interp(xk, fk, xnew)

Evaluate the barycentric interpolant that passes through `(xk, fk)` at the
query points `xnew`.  Returns a tuple `(values, P)` where `P` is the
interpolation matrix.
"""
function Bary_Interp(xk, fk, xnew)
    xnodes = vec(float.(xk))
    fvals = vec(float.(fk))
    xquery = vec(float.(xnew))

    length(xnodes) == length(fvals) || error("xk and fk must have the same length.")
    N = length(xnodes)

    w = barycentric_weights(xnodes)
    values = zeros(eltype(xnodes), length(xquery))
    P = zeros(eltype(xnodes), length(xquery), N)

    atol = 10 * eps(eltype(xnodes))
    for (i, xq) in enumerate(xquery)
        idx = findfirst(t -> isapprox(xq, t; atol = atol), xnodes)
        if idx !== nothing
            values[i] = fvals[idx]
            P[i, :] .= 0
            P[i, idx] = 1
            continue
        end
        diff = xq .- xnodes
        numer = w ./ diff
        denom = sum(numer)
        P[i, :] .= numer ./ denom
        values[i] = dot(P[i, :], fvals)
    end

    return values, P
end

"""
    Generalized_Diff_Mat(xs)

Construct the first-order differentiation matrix via barycentric weights for an
arbitrary node distribution `xs`.
"""
function Generalized_Diff_Mat(xs)
    x = vec(float.(xs))
    N = length(x)
    N >= 2 || error("At least two nodes are required to build a differentiation matrix.")
    w = barycentric_weights(x)
    D = zeros(eltype(x), N, N)
    @inbounds for i in 1:N
        for j in 1:N
            i == j && continue
            D[i, j] = (w[j] / w[i]) / (x[i] - x[j])
        end
        D[i, i] = -sum(D[i, k] for k in 1:N if k != i)
    end
    return D
end

"""
    Poly_Roots(vv)

Return the (complex) roots of the polynomial with coefficients `vv` ordered by
increasing real/imaginary parts.
"""
function Poly_Roots(vv)
    coeffs = vec(float.(vv))
    n = length(coeffs) - 1
    n < 1 && return ComplexF64[]

    coeffs_norm = coeffs ./ coeffs[end]
    M = zeros(eltype(coeffs_norm), n, n)
    for k in 1:(n-1)
        M[k, k+1] = 1.0
    end
    M[end, :] .= -coeffs_norm[1:end-1]
    roots = eigvals(M)
    sort!(roots, by = x -> (real(x), imag(x)))
    return roots
end

end
