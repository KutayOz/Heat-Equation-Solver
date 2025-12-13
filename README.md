# 1D Heat Equation Solver

A complete numerical methods program for solving the 1D heat (diffusion) equation using the generalized Crank-Nicolson (theta) scheme.

## Features

- **Part A: Von Neumann Stability Analysis**
  - Amplification factor derivation and visualization
  - Stability comparison across schemes (Explicit, Crank-Nicolson, Implicit)
  - Automatic stability report generation

- **Part B: Time-Marching Computations**
  - Three numerical schemes: Explicit (FTCS), Crank-Nicolson, Implicit (Backward Euler)
  - Custom Thomas algorithm (TDMA) for tridiagonal systems - no full matrices
  - Convergence to steady state with configurable tolerance
  - Behavioral diagnostics (overshoot, oscillation detection)

## Installation

### Requirements

- Python 3.11+
- NumPy
- Matplotlib

### Setup

```bash
cd heat_equation_solver
pip install -r requirements.txt
```

## Usage

### Interactive Menu (No Command-Line Arguments)

```bash
python3 -m src.menu
```

This launches an interactive menu where you can choose options without typing any arguments.

### Command-Line Interface

#### Run Everything

```bash
python3 -m src.main --run_all
```

### Run Specific Parts

```bash
# Stability analysis only
python -m src.main --part analysis_only

# Time-marching computations only
python -m src.main --part compute_only
```

### Custom Parameters

```bash
# Custom mu values
python -m src.main --mu_list 0.25,0.5,1.0,2.5,10,100,1000

# Custom convergence settings
python -m src.main --max_steps 100000 --tol 1e-6

# Custom grid resolution
python -m src.main --J 64
```

### All Options

```
--run_all           Run both stability analysis and computations
--part              Choose: analysis_only, compute_only, both
--mu_list           Comma-separated mu values for computations
--mu_analysis       Comma-separated mu values for stability plots
--max_steps         Maximum time steps (default: 200000)
--tol               Convergence tolerance (default: 1e-5)
--J                 Grid intervals (default: 32, giving 33 points)
--output_dir        Output directory (default: outputs)
--report_dir        Report directory (default: report)
--show_plots        Display plots interactively
```

## Project Structure

```
heat_equation_solver/
├── src/
│   ├── __init__.py           # Package initialization
│   ├── main.py               # CLI entry point
│   ├── stability_analysis.py # Von Neumann analysis (Part A)
│   ├── time_marching.py      # Time stepping (Part B)
│   ├── thomas_solver.py      # Tridiagonal solver (Thomas algorithm)
│   └── report_generator.py   # Automatic report generation
├── outputs/                   # Generated plots and CSV files
├── report/                    # Generated markdown report
├── README.md
└── requirements.txt
```

## Mathematical Background

### The Heat Equation

```
u_t = alpha * u_xx
```

where `u(x,t)` is temperature, `alpha` is thermal diffusivity.

### Theta Scheme Discretization

```
(u_i^(n+1) - u_i^n)/dt = alpha * [
    theta * (u_(i+1)^(n+1) - 2*u_i^(n+1) + u_(i-1)^(n+1)) / dx^2
  + (1-theta) * (u_(i+1)^n - 2*u_i^n + u_(i-1)^n) / dx^2
]
```

- `theta = 0`: Explicit (FTCS) - Forward Time, Central Space
- `theta = 0.5`: Crank-Nicolson - second-order accurate in time and space
- `theta = 1`: Implicit (Backward Euler) - first-order in time

### Key Parameter

```
mu = alpha * dt / dx^2    (mesh ratio / CFL number)
```

### Amplification Factor (Von Neumann Analysis)

```
s2 = sin^2(beta/2)
g(theta, mu, beta) = (1 - 4*(1-theta)*mu*s2) / (1 + 4*theta*mu*s2)
```

Stability requires `|g| <= 1` for all `beta` in `[0, pi]`.

## Implementation Details

### Grid Indexing

- Grid has `J+1` points: indices `j = 0, 1, ..., J`
- `j = 0` is left boundary (x = 0)
- `j = J` is right boundary (x = 1)
- Interior nodes: `j = 1, ..., J-1`
- Grid spacing: `dx = 1.0 / J`

### Time Step Derivation

From the definition of mu:
```
mu = alpha * dt / dx^2
dt = mu * dx^2 / alpha
```

With `alpha = 1.0` (default), `dt = mu * dx^2`.

### Boundary Conditions

Dirichlet conditions enforced at each time step:
- `u[0] = 1` (left boundary, fixed temperature)
- `u[J] = 0` (right boundary, fixed temperature)

### Thomas Algorithm

The tridiagonal system for implicit/CN schemes:
```
-lambda*mu * u[j-1]^(n+1) + (1 + 2*lambda*mu) * u[j]^(n+1) - lambda*mu * u[j+1]^(n+1) = RHS
```

Solved using O(n) Thomas algorithm with:
1. Forward elimination (remove lower diagonal)
2. Back substitution (solve from last equation)

**No full matrices are formed** - only 1D arrays for the three diagonals and RHS.

### Stability Behavior

| Scheme | Condition | Behavior |
|--------|-----------|----------|
| Explicit (theta=0) | mu <= 0.5 | Conditionally stable |
| Explicit (theta=0) | mu > 0.5 | UNSTABLE (exponential blowup) |
| Crank-Nicolson (theta=0.5) | any mu | Unconditionally stable, may oscillate |
| Implicit (theta=1) | any mu | Unconditionally stable, strongly damping |

### Why Explicit Becomes Unstable

For `theta = 0`, the amplification factor is:
```
g = 1 - 4*mu*sin^2(beta/2)
```

At `beta = pi` (highest frequency): `g = 1 - 4*mu`

For `|g| <= 1`: need `-1 <= 1 - 4*mu`, which gives `mu <= 0.5`.

For `mu > 0.5`, high-frequency modes have `|g| > 1` and grow exponentially.

### Why CN Can Be Stable But Oscillatory

For `theta = 0.5`:
```
g = (1 - 2*mu*s2) / (1 + 2*mu*s2)
```

- Always `|g| < 1` (stable)
- But `g` becomes negative when `1 - 2*mu*s2 < 0`
- Negative `g` means solution alternates sign each step (oscillatory)

## Design Decisions

1. **No full matrices**: All linear algebra uses 1D arrays only. The Thomas algorithm operates in O(n) time and space.

2. **Convergence criterion**: L2 norm (RMS) of change between time steps over interior nodes. Threshold: 1e-5 by default.

3. **Stability detection**: Solution marked unstable if `|u| > 1e6` or NaN/Inf detected.

4. **Overshoot detection**: Flagged if solution exceeds physical bounds [0, 1] at any time step.

5. **Oscillation detection**: Tracked via sign changes in incremental updates at probe points.

## Output Files

After running, you will find:

### In `outputs/`
- `amplification_mu_*.png` - Amplification factor plots for Part A
- `stability_regions.png` - Stability region visualization
- `evolution_*.png` - Solution evolution plots
- `final_solution_mu_*.png` - Final solution comparisons
- `convergence_history.png` - Residual history
- `convergence_results.csv` - Full results table

### In `report/`
- `report.md` - Complete markdown report with all analysis

## License

MIT License - Free for educational and research use.

## Author

Generated by numerical methods analysis tool.
