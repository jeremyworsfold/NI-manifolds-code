# Noise-induced manifolds: manuscript code

Julia code accompanying the manuscript *Noise Induced Manifolds: reducing SDEs to PDMPs*.
The repository contains the two models used in the paper:

- `NIManifoldsCode.TK2`: the two-species Togashi--Kaneko model; and
- `NIManifoldsCode.Ants`: foragers switching between two food sources while source A is
  depleted and recovers.

For each model, the code compares the exact microscopic Markov chain with the leading-order
piecewise-deterministic Markov process (PDMP) and its first-order finite-timescale correction.

## Requirements

- Julia 1.11 or 1.12
- Git

Clone the repository and instantiate the package:

```sh
git clone https://github.com/jeremyworsfold/NI-manifolds-code.git
cd NI-manifolds-code
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Dependency versions are bounded by `Project.toml`. Manifests are intentionally not committed,
so Julia can resolve compatible versions on both supported Julia releases.

## Reproducing the results

Instantiate the script environment once:

```sh
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

Quick mode uses smaller exact chains and is intended to check an installation:

```sh
julia --project=scripts scripts/reproduce.jl --mode quick --example all
```

Full mode uses the parameters and resolutions reported in the manuscript:

```sh
julia --project=scripts scripts/reproduce.jl --mode full --example all
```

Replace `all` with `tk2` or `ants` to run one model. Results are written beneath
`output/<mode>/<model>/`. Each run records its configuration and numerical checks in
`metrics.toml`; figure-ready CSV files are written below the corresponding `data/` directory.
Generated output is ignored by Git.

Optional Makie previews can be generated after a reproduction run:

```sh
julia --project=scripts/preview -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts/preview scripts/preview/preview.jl --mode quick --example all
```

## Ant-foraging notation

The code follows the manuscript notation. There are `N` foragers, `n` of which are at source
A, and `u` is the fraction of food remaining at A. Neutral recruitment occurs at rate
`beta*n*(N-n)` in each direction, with `epsilon = 1/beta`. Food changes according to

```text
du/dt = delta * (1 - u - n/N).
```

A lone forager leaves A at `omega0 + omega1*(1-u)^2` and leaves B at
`omega0 + omega1*u^2`. Public functions use the source labels `:A` and `:B`; these correspond
to the manuscript's `+` and `-` branches, respectively. The exact chain discretises the food
level into `food_steps = M` bounded upwind steps.

Example:

```julia
using NIManifoldsCode

p = Ants.ModelParameters(
    foragers = 200,
    epsilon = 0.1,
    delta = 0.22,
    omega0 = 0.02,
    omega1 = 0.2,
)

leading = Ants.stationary_solution(p)
corrected = Ants.stationary_solution(p; correction = :local)
chain = Ants.exact_chain(p; food_steps = 500)

Ants.source_probabilities(corrected)       # (P_A, P_B)
Ants.exact_residence(chain, :A).mean
Ants.wasserstein_distance(corrected, chain)
```

## TK2

TK2 provides `ModelParameters` for the microscopic chain, `ReducedParameters` for the reduced
theory, and leading, local and reset-corrected stationary solutions. For example:

```julia
using NIManifoldsCode

p = TK2.ReducedParameters(kappa = 0.4, beta = 0.0, epsilon = 0.1)
solution = TK2.stationary_solution(p; correction = :reset)

TK2.branch_probabilities(solution)
TK2.stationary_moments(solution)
```

## Tests and continuous integration

Run the package tests with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

GitHub Actions runs the package tests on Julia 1.11 and 1.12 and runs the quick reproduction
and preview workflows on Julia 1.12.

## License and citation

The code is released under the MIT License. If you use it, please cite the associated
manuscript and the software metadata in [`CITATION.cff`](CITATION.cff).
