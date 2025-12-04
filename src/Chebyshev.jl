module Chebyshev  # Início do módulo

using LinearAlgebra
using FFTW

export D_Cheby, ChebyQuadrature, ChebyshevTransform, Base_Cheby_n,
        M_Prod_Cheby, ChebyNodes, JS_Cheb



# ===================================================================
# MÉTODO 1: O "MOTOR" (SÓ ACEITA VETORES)
# ===================================================================
"""
    ChebyshevTransform(y_in::AbstractVector, type::Int, direction::Int)

Função "motor" unificada. Aplica a transformada de Chebyshev a um VETOR.

# Argumentos
- `y_in`: Vetor de entrada (Físico ou Espectral).
- `type`: Tipo de nó (1=Gauss, 2=Lobatto, 3=Radau-Esquerdo, 4=Radau-Direito).
- `direction`: 1 (Física -> Espectral) ou 2 (Espectral -> Física).
"""
function ChebyshevTransform(y_in::AbstractVector, type::Int, direction::Int)
    
    y = Float64.(y_in)
    N_pts = length(y) # Número de pontos

    if N_pts == 1
        return y
    end

    # --- TIPO 1: GAUSS (DCT-II) ---
    if type == 1
        N = N_pts # n=N pontos (ordem n-1)
        if direction == 1 # Física -> Espectral (DCT-II)
            yp = reverse(y)
            B = FFTW.r2r(yp, FFTW.REDFT10)
            B ./= N
            B[1] /= 2.0
            return B
        elseif direction == 2 # Espectral -> Física (IDCT-II)
            y_scaled = copy(y)
            y_scaled[1] *= 2.0
            y_out = FFTW.r2r(y_scaled, FFTW.REDFT01)
            y_out = reverse(y_out)
            y_out ./= 2.0 
            return y_out
        end

    # --- TIPO 2: LOBATTO (DCT-I) ---
    elseif type == 2
        N = N_pts # n+1 pontos (ordem n)
        n = N - 1
        if direction == 1 # Física -> Espectral
            yp = reverse(y)
            yh = FFTW.r2r(yp, FFTW.REDFT00)
            yh[2:n] ./= n 
            yh[[1, N]] ./= (2 * n)
            return yh
        elseif direction == 2 # Espectral -> Física
            y_scaled = copy(y)
            y_scaled[[1, N]] .*= 2.0
            y_out = FFTW.r2r(y_scaled, FFTW.REDFT00)
            y_out = reverse(y_out)
            y_out ./= 2.0 
            return y_out
        end

    # --- TIPO 3: RADAU-ESQUERDO (FFT) ---
    elseif type == 3
        n = N_pts - 1 # n+1 pontos (ordem n)
        M = 2*n + 1
        if direction == 1 # Física -> Espectral
            y_prime = copy(y)
            y_prime[1] /= 2.0 # Nó em k=0 (x=-1)
            g = zeros(ComplexF64, M)
            g[1] = y_prime[1] 
            g[2:n+1] = y_prime[2:n+1] ./ 2.0
            g[n+2:end] = g[n+1:-1:2]
            YH_fft = fft(g)
            yh_mod = real(YH_fft[1:N_pts]) .* (4.0 / M)
            yh_mod[1] /= 2.0
            yh = yh_mod .* (-1) .^ (0:n)
            return yh
        elseif direction == 2 # Espectral -> Física
            j_vec = (0:n)
            yh_mod = y .* (-1) .^ j_vec
            g = zeros(ComplexF64, M)
            g[1] = yh_mod[1]
            g[2:n+1] = yh_mod[2:end] ./ 2.0
            g[n+2:end] = g[n+1:-1:2]
            Y_fft = fft(g)
            y_out = real(Y_fft[1:N_pts])
            return y_out
        end

    # --- TIPO 4: RADAU-DIREITO (FFT) ---

    elseif type == 4
        n = N_pts - 1 # n+1 pontos (ordem n)
        M = 2*n + 1
        if direction == 1 # Física -> Espectral
            y_prime = copy(y)
            y_prime[end] /= 2.0 # Nó em k=0 (x=1), que é y[end]
            g = zeros(ComplexF64, M)
            g[1] = y_prime[end] 
            g[2:n+1] = y_prime[n:-1:1] ./ 2.0
            g[n+2:end] = g[n+1:-1:2]
            YH_fft = fft(g)
            yh_unscaled = real(YH_fft[1:N_pts])
            yh = yh_unscaled .* (4.0 / M)
            yh[1] /= 2.0
            return yh
        elseif direction == 2 # Espectral -> Física
            g = zeros(ComplexF64, M)
            g[1] = y[1]
            g[2:n+1] = y[2:end] ./ 2.0
            g[n+2:end] = g[n+1:-1:2]
            Y_fft = fft(g)
            y_natural_order = real(Y_fft[1:N_pts])
            y_out = reverse(y_natural_order)
            return y_out
        end
    
    else
        error("Tipo 'type' inválido. Use 1 (Gauss), 2 (Lobatto), 3 (Radau-L) ou 4 (Radau-R).")
    end
end


# ===================================================================
# MÉTODO 2: O "INVÓLUCRO" (ACEITA MATRIZES)
# ===================================================================
"""
    ChebyshevTransform(y_in::AbstractMatrix, type::Int, direction::Int)

Aplica a transformação coluna a coluna, evitando `mapslices`.
"""
function ChebyshevTransform(y_in::AbstractMatrix, type::Int, direction::Int)
    out = Array{Float64}(undef, size(y_in))
    for (j, col) in enumerate(eachcol(y_in))
        out[:, j] = ChebyshevTransform(col, type, direction)
    end
    return out
end


# ====================================================================
# 2. Matriz de Diferenciação Modal de CHEBYSHEV I (D_Cheby)
#    Aplica-se aos coeficientes a[0] a a[N-1]
#    Dimensão: N x N
# ====================================================================
"""
    D_Cheby(N::Int)

Cria a matriz de diferenciação modal de Chebyshev I de dimensão N x N.
Os elementos D[i, j] transformam a[j] em a'[i].
Indices (i, j) vão de 1 a N, correspondendo aos índices modais 0 a N-1.
"""
function D_Cheby(n::Int)
    N = n+1;
    D = zeros(Float64, N, N)

    # i e j são os índices baseados em 1 (1 a N).
    # O índice modal (ordem do polinômio) é i_modal = i - 1 e j_modal = j - 1.
    for i in 1:N
        i_modal = i - 1
        for j in (i + 1):N
            j_modal = j - 1
            
            # Condição de paridade: i_modal + j_modal deve ser ímpar
            if isodd(i_modal + j_modal) 
                
                # Caso i_modal = 0 (Primeira linha da matriz - i=1)
                if i_modal == 0
                    # Condição adicional: j_modal deve ser ímpar.
                    if isodd(j_modal)
                        # Elemento é j_modal
                        D[i, j] = j_modal
                    end
                # Caso 0 < i_modal < j_modal (Linhas seguintes - i > 1)
                else
                    # Elemento é 2 * j_modal
                    D[i, j] = 2 * j_modal
                end
            end
        end
    end
    return D
end


"""
    ChebyQuadrature(a,b,n::Int, tipo::Int = 1) - Clenshawn style

Calcula os nós (z) e pesos (w) para integração espectral de f(x)
no intervalo [-1, 1] (ou seja, ∫ f(x) dx) usando `n` pontos.
tipos 
1. Cheby-Gauss         2. Cheby-Lobatto
3. Radau [-1,1)        4. Radau (-1,1]


"""
function ChebyQuadrature(a::Real, b::Real,n::Int, tipo::Int = 1)
    
    local z, w

    if n == 0
        error("Número de pontos 'n' deve ser >= 1.")
    elseif n == 1
        w = [2.0]
        if tipo == 1
            z = [0.0]
        elseif tipo == 2
            z = [0.0]
        elseif tipo == 3
            z = [-1.0] # Nó de Radau padrão
        end
        return (z, w)
    end


    if tipo == 1
        # --- 1. Quadratura de Fejér (Tipo 1) ---
        # Gauss-Cheby
        k = 1:n
        z = @. cos((2*k - 1) * pi / (2*n))
        theta = @. (2*k - 1) * pi / (2*n)
        
        w = zeros(n)
        for k_idx in 1:n
            theta_k = theta[k_idx]
            soma_j = 0.0
            for j in 1:floor(Int, (n-1)/2)
                soma_j += cos(2*j * theta_k) / (4*j^2 - 1)
            end
            w[k_idx] = (2/n) * (1 - 2 * soma_j)
        end

    elseif tipo == 2
        # --- 2. Quadratura de Clenshaw-Curtis (Trefethen O(N^2)) ---
        # Lobatto-Cheby
        N = n - 1 # Ordem do Polinômio
        
        k = 0:N
        theta = @. pi * k / N
        z = - cos.(theta)
        
        w = zeros(n)
        ii = 2:n-1
        v = ones(n - 2)
        
        theta_ii = theta[ii]
        
        if iseven(N)
            w[1] = 1.0 / (N^2 - 1)
            w[end] = w[1]
            for k_ = 1:(N/2 - 1)
                v .-= 2 .* cos.(2 * k_ .* theta_ii) ./ (4*k_^2 - 1)
            end
            v .-= cos.(N .* theta_ii) ./ (N^2 - 1)
        else
            w[1] = 1.0 / N^2
            w[end] = w[1]
            for k_ = 1:((N-1)/2)
                v .-= 2 .* cos.(2 * k_ .* theta_ii) ./ (4*k_^2 - 1)
            end
        end
        w[ii] .= 2 .* v ./ N
        
    elseif tipo == 3
        # --- 3. Quadratura de Fejér-Radau (Nó em -1) x ∈ [-1,1)  ---
        # (Implementação O(N^2) correta, baseada no seu código MATLAB)
        
        N_poly = n - 1  # Ordem do polinômio (n pontos)
        M = 2*N_poly + 1 # = 2*n - 1
        
        # 1. Nós
        k = 0:N_poly
        theta = (2 * pi / M) .* k
        z = -cos.(theta) # z[1] é -1
        
        # 2. Pesos
        w = zeros(n)
        
        # Coeficientes para j par: 2, 4, ..., N_poly
        j_even = 2:2:N_poly
        coeffs = 2.0 ./ (1.0 .- j_even.^2)
        
        # Loop sobre cada nó `k`
        for k_idx in 1:n
            theta_k = theta[k_idx]
            
            # Soma interna (só sobre j par)
            inner_sum = sum(coeffs .* cos.(j_even .* theta_k))
            
            # Termo constante + soma
            w[k_idx] = 1.0 + inner_sum
        end
        
        # 3. Normalização final
        w .*= (4.0 / M)
        
        # 4. Ajuste para o ponto x = -1
        w[1] /= 2.0
        
        elseif tipo == 4
            # --- 4. Quadratura de Fejér-Radau (Nó em -1) x ∈ (-1,1]  ---
            # Faz a tipo 3 reverte.
        # (Implementação O(N^2) correta, baseada no seu código MATLAB)
        
        N_poly = n - 1  # Ordem do polinômio (n pontos)
        M = 2*N_poly + 1 # = 2*n - 1
        
        # 1. Nós
        k = 0:N_poly
        theta = (2 * pi / M) .* k
        z = -cos.(theta) # z[1] é -1
        
        # 2. Pesos
        w = zeros(n)
        
        # Coeficientes para j par: 2, 4, ..., N_poly
        j_even = 2:2:N_poly
        coeffs = 2.0 ./ (1.0 .- j_even.^2)
        
        # Loop sobre cada nó `k`
        for k_idx in 1:n
            theta_k = theta[k_idx]
            
            # Soma interna (só sobre j par)
            inner_sum = sum(coeffs .* cos.(j_even .* theta_k))
            
            # Termo constante + soma
            w[k_idx] = 1.0 + inner_sum
        end
        
        # 3. Normalização final
        w .*= (4.0 / M)
        
        # 4. Ajuste para o ponto x = -1
        w[1] /= 2.0
        w = reverse(w)
        z = -reverse(z)

    else
        error("Tipo de quadratura inválido: $tipo. Use 1, 2, ou 3.")
    end
    
    z_scal = @. ((b - a) * z + b + a) / 2
    w_scal = @. w * (b - a) / 2

    return (z_scal, w_scal)
end

"""
    Base_Cheby_n(n::Int, x::AbstractVector)

Constrói a matriz de Vandermonde de Chebyshev até a ordem `n` para o vetor `x`.
Retorna uma matriz onde a coluna `j+1` contém os valores de `T_j(x)`.

# Argumentos
- `n`: Ordem máxima do polinómio (a matriz terá n+1 colunas).
- `x`: Vetor de coordenadas.

# Retorno
- `Bn`: Matriz (length(x) x n+1).
"""
function Base_Cheby_n(n::Int, x::AbstractVector)
    # Garante que os valores sejam float (para evitar erros se x for Int)
    vals = float.(x)
    nL = length(vals)
    
    # Pré-aloca a matriz completa com zeros
    # eltype(vals) garante que se x for Float32, a matriz também será
    Bn = zeros(eltype(vals), nL, n + 1)

    # T_0(x) = 1 (Primeira coluna)
    Bn[:, 1] .= 1.0

    # Se n=0, terminamos aqui.
    if n >= 1
        # T_1(x) = x (Segunda coluna)
        Bn[:, 2] .= vals
    end

    # Loop da recorrência: T_k = 2x * T_{k-1} - T_{k-2}
    # Começa de k=2 (calculando a coluna 3, que é T_2)
    for k in 2:n
        # O macro @. aplica a operação elemento a elemento (broadcasting)
        # Coluna k+1 refere-se a T_k
        # Coluna k   refere-se a T_{k-1}
        # Coluna k-1 refere-se a T_{k-2}
        @. Bn[:, k+1] = 2 * vals * Bn[:, k] - Bn[:, k-1]
    end

    return Bn
end


"""
    M_Prod_Cheby(ah_in::AbstractVector, truncate::Bool=false)

Constrói a matriz de multiplicação M(ah) de Townsend.
"""

function M_Prod_Cheby(ah_in::AbstractVector; truncate::Bool=false)
    
    Nb = length(ah_in)

    ah = float.(vec(ah_in))
    T = eltype(ah) 
    
    Na = length(ah)
    N_full = Na + Nb - 1
    nrows = truncate ? Nb : N_full
    
    # Vetores auxiliares
    vec_toep = [2 * ah[1]; ah[2:end]; zeros(T, Nb - 1)]
    vec_hank = [ah[2:end]; zeros(T, Nb)]
    
    Ma = zeros(T, nrows, Nb)
    
    for j in 1:Nb
        for i in 1:nrows
            # --- Parte Toeplitz (B1) ---
            idx_t = abs(i - j) + 1
            val_t = (idx_t <= length(vec_toep)) ? vec_toep[idx_t] : zero(T)
            
            # --- Parte Hankel (B2) ---
            # O Hankel original começa na linha 2 da matriz final (devido ao shift).
            # Portanto, se i=1, a contribuição é ZERO.
            val_h = zero(T)
            
            if i > 1 # <--- CORREÇÃO CRUCIAL: Só calcula Hankel se i > 1
                idx_h = i + j - 2
                if idx_h >= 1 && idx_h <= length(vec_hank)
                    val_h = vec_hank[idx_h]
                end
            end
            
            Ma[i, j] = 0.5 * (val_t + val_h)
        end
    end
    
    return Ma
end

"""
    ChebyNodes(n::Int, opt::Int)

Calcula n+1 nós de Chebyshev e os seus pesos de quadratura.
"""
function ChebyNodes(n::Int, opt::Int)
    k = 0:n
    local x, w

    if opt == 1
        # 1. Gauss (zeros de T_{n+1}, n+1 pontos)
        x = -cos.((k .+ 0.5) .* π / (n + 1))
        w = ones(n + 1) .* (π / (n + 1))
        
    elseif opt == 2
        # 2. Lobatto (extremos de T_n)
        x = -cos.(k .* π / n)
        w = ones(n + 1) .* (π / n)
        w[1] /= 2
        w[end] /= 2 
        
    elseif opt == 3
        # 3. Radau-Esquerdo (nó em x=-1)
        x = -cos.((2 .* k .* π) / (2n + 1)) 
        w = ones(n + 1) .* (2 * π / (2n + 1))
        w[1] /= 2
        
    elseif opt == 4
        # 4. Radau-Direito (nó em x=1)
        x = cos.(2 .* (n .- k) .* π / (2n + 1)) # k=(n:-1:0)
        w = ones(n + 1) .* (2 * π / (2n + 1))
        w[end] /= 2
    else
        error("Opção 'opt' inválida. Use 1, 2, 3 ou 4.")
    end
    
    return x, w # Retorna os nós e os pesos da QUADRATURA (não os pesos da Fejér)
end

"""
    JS_Cheb(N::Int)

Constrói a matriz operacional de integração para polinômios de Chebyshev.

# Argumentos
- `N`: Ordem da expansão de entrada. A matriz será (N+2) x (N+1).
"""
function JS_Cheb(N::Int)
    # A matriz mapeia N+1 coefs -> N+2 coefs
    J = zeros(N + 2, N + 1)

    # --- Caso k=0 (Coluna 1) ---
    # Integral de T_0(x)=1 é x = T_1(x).
    # Ajuste da constante depois.
    J[2, 1] = 1.0 

    # --- Caso k=1 (Coluna 2) ---
    # Integral de T_1(x)=x é x^2/2 = T_2(x)/4 + T_0(x)/4.
    # O termo T_2 (Linha 3) recebe 1/4.
    J[3, 2] = 0.25

    # --- Loop para k >= 2 (Colunas 3 em diante) ---
    for k in 2:N
        # Coluna j = k+1
        col = k + 1
        
        # Termo T_{k+1} (Linha k+2): Coeficiente 1 / (2*(k+1))
        J[k+2, col] = 1.0 / (2 * (k + 1))
        
        # Termo T_{k-1} (Linha k): Coeficiente -1 / (2*(k-1))
        J[k, col] = -1.0 / (2 * (k - 1))
    end

    # --- Ajuste da Constante de Integração (Linha 1 / Coef de T_0) ---
    # Condição F(-1) = 0.
    # T_n(-1) = (-1)^n.
    # Queremos sum(c_i * (-1)^(i-1)) = 0
    # c_1 = - sum_{i=2}^{N+2} c_i * (-1)^(i-1)
    
    for j in 1:(N+1)
        sum_alt = 0.0
        for i in 2:(N+2)
            # i-1 é o grau do polinômio T associado à linha i
            sum_alt += J[i, j] * (-1)^(i-1) 
        end
        J[1, j] = -sum_alt
    end

    return J
end



end  # Módulo
