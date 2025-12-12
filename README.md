# SpectralTools – Plano do Projeto Final

Este documento descreve o plano de desenvolvimento do pacote **SpectralTools.jl**,
utilizado como projeto final da disciplina **Introdução aos Métodos Espectrais**.

---

## Início Rápido

```julia
julia --project=.
using Pkg
Pkg.instantiate()
using SpectralTools
```

O pacote exporta **submódulos** em vez de colocar todas as funções no escopo
principal. Após `using SpectralTools`, importe cada ferramenta de forma
qualificada, por exemplo:

```julia
using SpectralTools.Chebyshev: ChebyshevTransform, ChebyQuadrature
using SpectralTools.BVP: solve_linear_bvp
```

---

## Slides Recomendados por Dupla

Os PDFs estão em `class_slides/`. Cada dupla deve usar os seguintes arquivos como
referência principal:

1. **Dupla 1 – Chebyshev + Difusão**  
   `PE_Aula_05_N.pdf` (pp. 2–5) para colocation em Chebyshev e `PE_Aula_06_N.pdf`
   (pp. 2–3) para Newton e MOL + RK4.
2. **Dupla 2 – Legendre + BVP não linear**  
   `PE_Aula_07_N.pdf`, `PE_Aula_08_N.pdf` (bases Legendre) e `PE_Aula_06_N.pdf`
   (fluxo de Newton) para lidar com não linearidades.
3. **Dupla 3 – Onda 1D**  
   `PE_Aula_10_N.pdf` (pp. 2–3) detalha o esquema leapfrog e as restrições CFL.
4. **Dupla 4 – Poisson 2D + mapeamentos**  
   `PE_Aula_09_N.pdf` (pp. 13–19) cobre Kronecker/Sylvester e o “Programa 16”; as
   pp. 8–10 trazem mapas conformes opcionais.

---

## 1. Estrutura do Pacote

```
SpectralTools/
├── src/
│   ├── SpectralTools.jl       # módulo principal e reexportações
│   ├── Chebyshev.jl           # ferramentas base Chebyshev
│   ├── Legendre.jl            # ferramentas base Legendre
│   ├── Fourier.jl             # ferramentas base Fourier
│   ├── Generic.jl             # utilidades genéricas
│   ├── BoundaryConditions.jl  # normalização de BCs
│   ├── Collocation.jl         # tipo SpectralGrid + grades Lobatto
│   ├── BVP.jl                 # solvers de BVP (novo)
│   └── PDE.jl                 # solvers de EDP (novo)
├── test/                     # suíte de testes
├── docs/                     # notas/Documenter
├── examples/                 # scripts de experimento
└── class_slides/             # PDFs das aulas
```

O módulo `SpectralTools.jl` expõe todos os submódulos e reexporta os solvers de
alto nível (`solve_linear_bvp`, `solve_diffusion_1d`, etc.).

---

## 2. Solvers e Utilidades

- `Collocation.build_grid(N; basis, domain)` – constrói um `SpectralGrid` com nós
  Lobatto e matrizes `D₁`, `D₂` já escaladas.
- `BVP.solve_linear_bvp` – resolve `a(x) y'' + b(x) y' + c(x) y = f(x)` com BCs
  Dirichlet/Neumann/Robin.
- `BVP.solve_nonlinear_bvp` – Newton amortecido para `y'' = g(x, y, y')` com apoio
  opcional de derivadas analíticas `dg_dy`, `dg_dyp`.
- `PDE.solve_diffusion_1d` – metodologia das linhas + RK4 com enforcement de BCs.
- `PDE.solve_wave_1d` – leapfrog clássico com suporte a Dirichlet e Neumann.
- `PDE.solve_poisson_2d` – solução de Poisson 2D via `kron(I, Dy²) + kron(Dx², I)`.

Scripts em `examples/` demonstram o uso das rotinas (BVP linear, difusão etc.).

---

## 3. Código Base (não reimplementar)

- `src/Generic.jl`: `Bary_Interp`, `Generalized_Diff_Mat`, `Poly_Roots`.
- `src/Chebyshev.jl`: `ChebyshevTransform`, `D_Cheby`, `ChebyQuadrature`, `Base_Cheby_n`,
  `M_Prod_Cheby`, `ChebyNodes`, `JS_Cheb`.
- `src/Legendre.jl`: `eval_legendre`, `eval_legendre_p`, bases Gauss/Lobatto/Radau,
  `GaussQuadTypes`, `DS_Legendre`, `Base_Legendre_n`, `M_Prod_Legendre`, `JS_Leg`.
- `src/Fourier.jl`: `DS_Fourier`, `fourdif`, `Multi_Diff_Mat`, `Bary_Trig_Mat`,
  `JS_Fourier`, `fourier_quad`, `real_fft`, `real_ifft`.

---

## 4. Novos Solvers

### `src/BVP.jl`
1. Construir `SpectralGrid` (Chebyshev ou Legendre) com `build_grid`.
2. Montar as matrizes `D₁`, `D₂` e o operador linear `L`.
3. Aplicar condições de contorno via `BoundaryConditions`.
4. Resolver o sistema linear ou iterar Newton para o caso não linear.

### `src/PDE.jl`
- `solve_diffusion_1d`: reproduzir o problema analítico `u = e^{-π² t /4} sin(π(x+1)/2)`.
- `solve_wave_1d`: leapfrog + CFL com acompanhamento de energia.
- `solve_poisson_2d`: problema manufaturado (Programa 16) com Kronecker/Sylvester.

---

## 5. Fluxo de Trabalho

1. Clonar o repositório e criar uma branch `feature/parX-topico`.
2. Implementar as funções atribuídas (ver seção 6).
3. Adicionar scripts em `examples/` para cada experimento exigido.
4. Incluir pelo menos 3 testes novos por dupla.
5. Rodar `julia --project=.\ninclude("test/runtests.jl")` antes do PR.
6. Abrir PR descrevendo funcionalidades, scripts e testes.

---

## 6. Distribuição das Tarefas (arquivos/funções)

- **Dupla 1 – Chebyshev + Difusão** *(base: `PE_Aula_05_N.pdf`, `PE_Aula_06_N.pdf`)*  
  - `src/BVP.jl`: completar `solve_linear_bvp` usando grades de Chebyshev e o fluxo de
    imposição de BCs descrito nos slides da Aula 05.  
  - `src/PDE.jl`: implementar `solve_diffusion_1d` com MOL + RK4 e suporte a
    Dirichlet/Neumann.  
  - `examples/`: manter/estender `bvp_linear.jl` e `diffusion.jl` com análises de erro.

- **Dupla 2 – Legendre + BVP não linear** *(base: `PE_Aula_07_N.pdf`, `PE_Aula_08_N.pdf`, `PE_Aula_06_N.pdf`)*  
  - `src/BVP.jl`: adicionar suporte `basis = :legendre` em `solve_linear_bvp` e
    finalizar `solve_nonlinear_bvp` (Newton amortecido com `dg_dy`, `dg_dyp`).  
  - `test/`: criar casos específicos validando ortogonalidade Legendre e a
    convergência dos BVPs linear e não linear.  
  - `examples/`: script comparando Chebyshev × Legendre conforme Aulas 07/08.

- **Dupla 3 – Onda 1D** *(base: `PE_Aula_10_N.pdf`, pp. 2–3)*  
  - `src/PDE.jl`: preencher `solve_wave_1d` com o esquema leapfrog, inicialização
    do meio passo e opções Dirichlet/Neumann.  
  - `examples/`: scripts de pulso viajante e estudo CFL com medição de energia.  
  - `test/`: validar conservação aproximada de energia e uma solução manufaturada.

- **Dupla 4 – Poisson 2D + mapeamentos** *(base: `PE_Aula_09_N.pdf`, pp. 8–19)*  
  - `src/PDE.jl`: concluir `solve_poisson_2d` com a forma de Sylvester
    `kron(I, Dy^2) + kron(Dx^2, I)` e BCs Dirichlet dos slides.  
  - Opcional: adicionar `src/Mapping.jl` (incluído em `SpectralTools.jl`) com utilitários
    de mapas conformes simples (ex.: `x = sin(ξ)`).  
  - `examples/` + `test/`: reproduzir o Programa 16 e estudar variação de erro.

Cada dupla deve entregar implementações, experimentos reprodutíveis, testes e um
resumo no PR.

---

## 7. Guia de Estilo

- Não adicionar dependências fora da stdlib (LinearAlgebra, FFTW etc.).
- Manter docstrings e exemplos curtos para funções públicas.
- Buscar estabilidade de tipos e reutilizar utilidades existentes.

## 9. Documentação para Usuários

- Consulte `docs/USER_GUIDE.md` para instruções detalhadas de instalação, API e progresso de exemplos/testes.
- A versão web em `docs/web/index.html` (publicável via GitHub Pages) mostra um resumo interativo com figuras.

## 10. Relatório Completo de Testes + Cobertura

Execute o script abaixo para gerar um relatório consolidado (status dos testsets,
falhas detalhadas e cobertura de código). Recomendado rodar com cobertura ativada:

```bash
julia --project --code-coverage=user scripts/run_full_report.jl
```

O relatório indicará quais testsets passaram, quais falharam e exibirá um resumo
de cobertura por arquivo (requer que `.cov` sejam gerados via `--code-coverage`).

---

## 8. Entregáveis

1. Implementações concluídas nas funções designadas.  
2. Scripts de exemplo em `examples/`.  
3. Três testes novos por dupla em `test/`.  
4. PR com descrição, instruções de execução e resumo dos testes.

Seguindo o plano, o pacote resultará em uma biblioteca Julia coesa e testada
para experimentos em métodos espectrais.
