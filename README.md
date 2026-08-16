# Bayes-Joint-Model-ctDNA-PFS

Stan/mrgsolve code for the joint longitudinal-survival model described in:

> Johnson M, Serra Traynor C, Vishwanathan K, Overend P, Hartmaier R, Markovets A, Chmielecki J, Mugundu GM, Barrett JC, Tomkinson H, Ramalingam SS. **Longitudinal Circulating Tumor DNA Modeling to Predict Disease Progression in First-Line Mutant Epidermal Growth Factor Receptor Non-Small Cell Lung Cancer.** *Clin Pharmacol Ther.* 2024;115(2):349-360. doi: [10.1002/cpt.3113](https://doi.org/10.1002/cpt.3113). PMID: [38010260](https://pubmed.ncbi.nlm.nih.gov/38010260/).

The model jointly describes longitudinal circulating tumor DNA (ctDNA) kinetics and progression-free survival (PFS) in a Bayesian nonlinear mixed-effects framework, linking a biexponential (growth/decay) ctDNA trajectory to the hazard of disease progression.

## Repository structure

```
Models/
  jmODEklog.stan   Joint model, ODE formulation (solved with Torsten's
                    pmx_solve_group_rk45), NONMEM-style dosing/event data.
  jmStein.stan      Joint model, analytical (closed-form) Stein–Fojo
                     biomarker trajectory with Weibull survival. Plain
                     Stan (>= 2.33), no Torsten dependency.

Simulations/
  simjmklog.cpp     mrgsolve simulation model corresponding to the ODE
                     joint model, used for simulation-based validation.
```

Both Stan models share the same biological structure:

- **Biomarker (ctDNA) sub-model**: growth rate `kg` and decay rate `kd`
  govern a biexponential trajectory, `jmODEklog.stan` solves this as an
  ODE while `jmStein.stan` uses the closed-form Stein–Fojo solution.
- **Survival sub-model**: a Weibull hazard for time to disease
  progression, linked to the biomarker sub-model through shared
  individual-level parameters (`kg`, `kd`, baseline) via a log-linear
  predictor `eta`.
- **Between-subject variability**: modeled with a multivariate
  log-normal random-effects structure (Cholesky-factorized correlation).

## Requirements

- [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) / [Stan](https://mc-stan.org/) (>= 2.33 for `jmStein.stan`)
- [Torsten](https://metrumresearchgroup.github.io/Torsten/) (for `jmODEklog.stan`, which calls `pmx_solve_group_rk45`)
- [mrgsolve](https://mrgsolve.org/) (R package, for `Simulations/simjmklog.cpp`)

## Usage

These are model definition files, not a packaged pipeline. To fit a model,
supply the data blocks documented at the top of each `.stan` file (ctDNA
observations, survival times/censoring, treatment covariate, and prior
hyperparameters) via your preferred Stan interface (CmdStan, CmdStanR,
CmdStanPy, etc.), e.g.:

```r
library(StanFlowR)
cmdstan_mkdir("Models/jmStein.stan")
```

`Models/jmODEklog.stan` additionally requires linking against Torsten to
resolve `pmx_solve_group_rk45`.

## License

Released under the [MIT License](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff), or cite the publication above.
