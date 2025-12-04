
module SpectralTools

include("Generic.jl")
include("Chebyshev.jl")
include("Legendre.jl")
include("Fourier.jl")
include("BoundaryConditions.jl")
include("Collocation.jl")
include("BVP.jl")
include("PDE.jl")

using .Generic
using .Chebyshev
using .Legendre
using .Fourier
using .BoundaryConditions
using .Collocation
using .BVP
using .PDE

export Generic, Chebyshev, Legendre, Fourier,
       BoundaryConditions, Collocation,
       solve_linear_bvp, solve_nonlinear_bvp,
        solve_diffusion_1d, solve_wave_1d, solve_poisson_2d

end
