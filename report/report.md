# 1D Heat Equation Solver Report
## Generalized Crank-Nicolson (Theta) Scheme
*Generated: 2025-12-14 03:59:59*

---

## Table of Contents
1. [Introduction](#introduction)
2. [Part A: Von Neumann Stability Analysis](#von-neumann-stability-analysis)
3. [Part B: Time-Marching Computations](#time-marching-computations)
4. [Convergence Results](#convergence-results)
5. [Conclusions](#conclusions)

---

## Introduction

This report presents the analysis and numerical solution of the 1D heat (diffusion) equation using the generalized Crank-Nicolson (theta) scheme.

### The Heat Equation

We solve the partial differential equation:

```
u_t = alpha * u_xx
```

where:
- `u(x,t)` is the temperature distribution
- `alpha` is the thermal diffusivity
- Subscripts denote partial derivatives

### The Theta Scheme

The generalized Crank-Nicolson discretization is:

```
(u_i^(n+1) - u_i^n)/dt = alpha * [
    theta * (u_(i+1)^(n+1) - 2*u_i^(n+1) + u_(i-1)^(n+1)) / dx^2
  + (1-theta) * (u_(i+1)^n - 2*u_i^n + u_(i-1)^n) / dx^2
]
```

Special cases:
- `theta = 0`: Explicit (FTCS) - Forward Time, Central Space
- `theta = 0.5`: Crank-Nicolson - second-order accurate
- `theta = 1`: Implicit (Backward Euler)

Key parameter: `mu = alpha * dt / dx^2` (mesh ratio / CFL number)

---


### Stability Plots

---

## Time-Marching Computations
### Problem Setup
- Domain: [0, 1]
- Grid points: J = 8, giving 9 total points
- Grid spacing: dx = 1/8 = 0.125000
- Diffusion coefficient: alpha = 1.0
- Time step: dt = mu * dx^2 (varies with mu)

### Boundary Conditions
- Left (x=0): u = 1 (Dirichlet)
- Right (x=1): u = 0 (Dirichlet)

### Initial Condition
Step profile:
- u(x) = 1 for x <= 0.5
- u(x) = 0 for x > 0.5

### Key Observations

#### Explicit Scheme (theta = 0)
- Stable for mu in {0.25, 0.5}
- **Unstable** for mu in {1.0, 2.5, 10.0, 100.0, 1000.0}
- Confirms theoretical stability limit: mu <= 0.5
- Fastest per-step computation (no matrix solve)
- Requires many more steps for large dt

#### Crank-Nicolson (theta = 0.5)
- Unconditionally stable for all tested mu values
- Shows overshoot for mu in {100.0, 1000.0}
- Oscillatory behavior for mu in {2.5, 10.0, 100.0, 1000.0}
- Second-order accurate, often the best balance of accuracy and efficiency

#### Implicit Scheme (theta = 1)
- Unconditionally stable for all tested mu values
- No oscillatory behavior (g always positive)
- No overshoot observed - monotonic approach to steady state
- First-order accurate in time, may introduce excessive numerical diffusion

### Behavior Summary

| Behavior | Explicit | Crank-Nicolson | Implicit |
|----------|----------|----------------|----------|
| Stability | mu <= 0.5 | Unconditional | Unconditional |
| Oscillations | No (when stable) | Yes (large mu) | No |
| Accuracy | O(dt, dx^2) | O(dt^2, dx^2) | O(dt, dx^2) |
| Matrix solve | No | Yes (Thomas) | Yes (Thomas) |


### Computation Plots

#### Evolution Explicit
![Evolution Explicit](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/evolution_explicit.png)

#### Evolution Crank-Nicolson
![Evolution Crank-Nicolson](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/evolution_crank-nicolson.png)

#### Evolution Implicit
![Evolution Implicit](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/evolution_implicit.png)

#### Final Solution Mu 0.25
![Final Solution Mu 0.25](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_0.25.png)

#### Final Solution Mu 0.50
![Final Solution Mu 0.50](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_0.50.png)

#### Final Solution Mu 1.00
![Final Solution Mu 1.00](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_1.00.png)

#### Final Solution Mu 2.50
![Final Solution Mu 2.50](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_2.50.png)

#### Final Solution Mu 10.00
![Final Solution Mu 10.00](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_10.00.png)

#### Final Solution Mu 100.00
![Final Solution Mu 100.00](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_100.00.png)

#### Final Solution Mu 1000.00
![Final Solution Mu 1000.00](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/final_solution_mu_1000.00.png)

#### Convergence History
![Convergence History](../../../../../../../var/folders/_f/38hnp529435cdmt7y6j_lw6r0000gn/T/tmp_wtf29d8/convergence_history.png)

---

## Convergence Results

| Scheme | theta | mu | dt | Steps | Status | Overshoot | Oscillatory |
|--------|-------|-----|-----|-------|--------|-----------|-------------|
| Explicit (FTCS) | 0.0 | 0.25 | 3.91e-03 | 153 | OK | No | No |
| Explicit (FTCS) | 0.0 | 0.5 | 7.81e-03 | 125 | OK | No | Yes |
| Explicit (FTCS) | 0.0 | 1.0 | 1.56e-02 | - | UNSTABLE | Yes | No |
| Explicit (FTCS) | 0.0 | 2.5 | 3.91e-02 | - | UNSTABLE | Yes | No |
| Explicit (FTCS) | 0.0 | 10.0 | 1.56e-01 | - | UNSTABLE | Yes | No |
| Explicit (FTCS) | 0.0 | 100.0 | 1.56e+00 | - | UNSTABLE | Yes | No |
| Explicit (FTCS) | 0.0 | 1000.0 | 1.56e+01 | - | UNSTABLE | Yes | No |
| Crank-Nicolson | 0.5 | 0.25 | 3.91e-03 | 156 | OK | No | No |
| Crank-Nicolson | 0.5 | 0.5 | 7.81e-03 | 87 | OK | No | No |
| Crank-Nicolson | 0.5 | 1.0 | 1.56e-02 | 49 | OK | No | No |
| Crank-Nicolson | 0.5 | 2.5 | 3.91e-02 | 25 | OK | No | Yes |
| Crank-Nicolson | 0.5 | 10.0 | 1.56e-01 | 96 | OK | No | Yes |
| Crank-Nicolson | 0.5 | 100.0 | 1.56e+00 | 949 | OK | Yes | Yes |
| Crank-Nicolson | 0.5 | 1000.0 | 1.56e+01 | >1000 | NO CONV | Yes | Yes |
| Implicit (BE) | 1.0 | 0.25 | 3.91e-03 | 158 | OK | No | No |
| Implicit (BE) | 1.0 | 0.5 | 7.81e-03 | 90 | OK | No | No |
| Implicit (BE) | 1.0 | 1.0 | 1.56e-02 | 52 | OK | No | No |
| Implicit (BE) | 1.0 | 2.5 | 3.91e-02 | 26 | OK | No | No |
| Implicit (BE) | 1.0 | 10.0 | 1.56e-01 | 11 | OK | No | No |
| Implicit (BE) | 1.0 | 100.0 | 1.56e+00 | 5 | OK | No | No |
| Implicit (BE) | 1.0 | 1000.0 | 1.56e+01 | 3 | OK | No | No |

---

## Conclusions

### Summary of Findings

1. **Stability Analysis Confirms Theory**
   - Explicit scheme is conditionally stable with limit `mu <= 0.5`
   - Crank-Nicolson and Implicit are unconditionally stable
   - CN can exhibit oscillations for large mu despite stability

2. **Computational Efficiency Trade-offs**
   - Explicit: cheapest per step, but limited by stability
   - Implicit: allows large time steps, but first-order accurate
   - CN: best accuracy per computational cost for smooth solutions

3. **Solution Quality Considerations**
   - CN oscillations are bounded but may be undesirable
   - Implicit scheme provides monotonic convergence
   - For stiff problems (large mu), implicit methods are preferred

### Recommendations

- For high accuracy: use Crank-Nicolson with moderate mu
- For robustness: use Implicit with any mu
- For efficiency with accuracy: use Explicit with mu near 0.5
- Avoid CN with very large mu if oscillations are problematic

---

*End of Report*
