# SpectralTools.jl

Biblioteca de rotinas para métodos espectrais em Julia. O projeto reúne
operações clássicas com bases de Chebyshev, Legendre e Fourier, matrizes de
diferenças, interpolação e, mais recentemente, quadraturas de
Chebyshev–Lobatto para integrais em $[-1,1]$.

## Requisitos

1. Julia 1.10+ (recomendado instalar com
   [`juliaup`](https://github.com/JuliaLang/juliaup)):
   ```bash
   curl -fsSL https://install.julialang.org | sh
   ```
2. `git` para clonar o repositório.

## Instalação

```bash
git clone https://github.com/stealth-lndrs/SpectralTools.git
cd SpectralTools
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Para desenvolver o pacote a partir de outro projeto, abra o REPL em seu projeto
de destino e execute:

```
pkg> dev /caminho/para/SpectralTools
pkg> instantiate
```

## Uso rápido

```julia
julia --project=.
julia> using SpectralTools
julia> x = cheb_lobatto_nodes(16);
julia> w = cheb_lobatto_weights(16);
julia> cheb_lobatto_quadrature(x -> sqrt(1 - x^2), 32)
1.5707602875684286
```

Execute `pkg> st` para verificar que o ambiente está ativo e `pkg> test` para
confirmar que a instalação está íntegra.

## Funcionalidades disponíveis

| Submódulo / Função                       | Descrição resumida                                                    |
| --------------------------------------- | --------------------------------------------------------------------- |
| `Chebyshev` (e.g. `DCT_Chb_Lob`)        | Transformadas discreta de Chebyshev e matrizes de derivadas.          |
| `Legendre`                              | Bases espectrais e matrizes diferenciais associadas aos polinômios.   |
| `Fourier`                               | Ferramentas para problemas periódicos (transformadas e interpolações).|
| `PseudoGen`                             | Funções auxiliares para geração de operadores espectrais genéricos.   |
| `ChebyshevLobattoQuadrature`            | Rotinas `cheb_lobatto_nodes`, `cheb_lobatto_weights` e                |
|                                         | `cheb_lobatto_quadrature` para integrais em $[-1,1]$.                 |

Consulte os docstrings de cada função com `?nome_da_função` no REPL.

## Criando novas funcionalidades

1. Crie um arquivo em `src/` e inclua-o em `src/SpectralTools.jl` via
   `include("MeuModulo.jl")`.
2. Declare um módulo interno e exporte apenas as funções necessárias.
3. Documente cada função com docstrings (`""" ... """`) e exemplos mínimos.
4. Acrescente testes em `test/`, registrando o novo arquivo em `test/runtests.jl`.
5. Execute `julia --project=. -e 'using Pkg; Pkg.test()'` antes de abrir um PR.

### Boas práticas de implementação

- Prefira funções puras e vetorizadas ou loops explícitos com `@inbounds`.
- Evite variáveis globais; use constantes ou argumentos.
- Para novas quadraturas ou operadores, comprove a correção em testes com
  polinômios de baixa ordem ou soluções analíticas.
- Mantenha compatibilidade com pelo menos a versão LTS do Julia.

## Testes e garantia de qualidade

- Execute sempre `pkg> test` depois de modificar o código.
- Para testar apenas um arquivo, use `julia --project=. test/runtests.jl`.
- Se o pacote for desenvolvido via `Pkg.dev`, certifique-se de ativar o
  ambiente correto (`--project=.`) ao rodar experimentos.

## Contribuindo

1. Faça um fork do repositório.
2. Crie um branch descritivo (`feature/nova-quadratura`).
3. Adicione código + testes + documentação.
4. Garanta que `Pkg.test()` passa sem falhas.
5. Abra um pull request descrevendo a motivação e os resultados.

Sugestões, correções e novas funcionalidades são bem-vindas!
