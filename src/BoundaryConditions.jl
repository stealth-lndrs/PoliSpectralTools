module BoundaryConditions

using LinearAlgebra

export SideSpec, normalize_1d_bc, normalize_2d_bc, eval_bc_value, eval_bc_time_derivative

"""
Side specification for boundary conditions.

`kind` is one of `:dirichlet`, `:neumann`, or `:robin`.
`value` stores the prescribed data (number or callable).
`alpha`, `beta` encode Robin coefficients (`alpha * u + beta * u' = value`).
"""
struct SideSpec{V}
    kind::Symbol
    value::V
    alpha::Float64
    beta::Float64
end

SideSpec(kind::Symbol, value) = SideSpec(kind, value, kind === :dirichlet ? 1.0 : 0.0,
                                         kind === :neumann ? 1.0 : 0.0)

function parse_side(side)
    if side isa Tuple
        kind = side[1]
        if kind === :dirichlet
            return SideSpec(:dirichlet, side[2], 1.0, 0.0)
        elseif kind === :neumann
            return SideSpec(:neumann, side[2], 0.0, 1.0)
        elseif kind === :robin
            length(side) < 4 && error("Robin BC requires (:robin, alpha, beta, value).")
            return SideSpec(:robin, side[4], float(side[2]), float(side[3]))
        else
            error("Unknown boundary type $kind.")
        end
    elseif side isa Number || side isa Function
        return SideSpec(:dirichlet, side, 1.0, 0.0)
    else
        error("Unsupported boundary specification: $side")
    end
end

"""
    normalize_1d_bc(bc)

Ensure both `left` and `right` boundaries are specified.
"""
function normalize_1d_bc(bc::NamedTuple)
    left = haskey(bc, :left) ? parse_side(bc.left) : SideSpec(:dirichlet, 0.0, 1.0, 0.0)
    right = haskey(bc, :right) ? parse_side(bc.right) : SideSpec(:dirichlet, 0.0, 1.0, 0.0)
    return (left = left, right = right)
end

"""
    normalize_2d_bc(bc)

Ensure `left`, `right`, `bottom`, and `top` entries exist (Dirichlet zero by default).
"""
function normalize_2d_bc(bc::NamedTuple)
    defaults = Dict(:left => SideSpec(:dirichlet, 0.0, 1.0, 0.0),
                    :right => SideSpec(:dirichlet, 0.0, 1.0, 0.0),
                    :bottom => SideSpec(:dirichlet, 0.0, 1.0, 0.0),
                    :top => SideSpec(:dirichlet, 0.0, 1.0, 0.0))
    return (left = parse_side(get(bc, :left, defaults[:left])),
            right = parse_side(get(bc, :right, defaults[:right])),
            bottom = parse_side(get(bc, :bottom, defaults[:bottom])),
            top = parse_side(get(bc, :top, defaults[:top])))
end

function call_value(val, x, t)
    if val isa Number
        return float(val)
    elseif val isa Function
        try
            return val(x, t)
        catch err
            if err isa MethodError
                try
                    return val(x)
                catch err2
                    if err2 isa MethodError
                        try
                            return val(t)
                        catch err3
                            if err3 isa MethodError
                                return val()
                            else
                                rethrow(err3)
                            end
                        end
                    else
                        rethrow(err2)
                    end
                end
            else
                rethrow(err)
            end
        end
    else
        error("Unsupported boundary value type $(typeof(val)).")
    end
end

eval_bc_value(spec::SideSpec, x, t) = call_value(spec.value, x, t)

function eval_bc_time_derivative(spec::SideSpec, x, t; ε = 1e-6)
    spec.kind === :neumann || spec.kind === :dirichlet || return 0.0
    if spec.value isa Number
        return 0.0
    end
    return (call_value(spec.value, x, t + ε) - call_value(spec.value, x, t - ε)) / (2ε)
end

end
