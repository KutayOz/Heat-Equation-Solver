"""
Time-Marching Computations for 1D Heat Equation
================================================

This module implements time-marching to steady state for the 1D heat equation
using the generalized Crank-Nicolson (theta) scheme.

Problem Setup:
    Domain: 0 <= x <= 1
    Grid: J = 32, so J+1 = 33 points (indices j = 0, 1, ..., J)
    dx = 1.0 / J
    
    Boundary conditions (Dirichlet):
        u[0]^n = 1  (left boundary)
        u[J]^n = 0  (right boundary)
    
    Initial condition (step profile):
        u[j]^0 = 1 for j <= J/2
        u[j]^0 = 0 for j > J/2

Time stepping:
    dt = mu * dx^2 / alpha
    where mu is the mesh ratio (user specified)

Convergence criterion:
    L2_norm(u^(n+1) - u^n) <= tol
    where L2_norm is RMS over interior nodes

Three schemes implemented:
    - Simple explicit (theta = 0.0): FTCS, conditionally stable
    - Crank-Nicolson (theta = 0.5): second-order accurate, unconditionally stable
    - Simple implicit (theta = 1.0): backward Euler, unconditionally stable
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
import os
import csv
import time as time_module

from .thomas_solver import thomas_solve, build_tridiagonal_coefficients, compute_rhs


@dataclass
class SimulationResult:
    """Container for simulation results and diagnostics."""
    scheme_name: str
    lam: float  # theta parameter
    mu: float
    dx: float
    dt: float
    J: int
    
    # Convergence info
    converged: bool
    n_steps: int
    final_residual: float
    
    # Stability info
    stable: bool
    failure_reason: str = ""
    
    # Behavior diagnostics
    overshoot: bool = False  # True if min(u) < 0 or max(u) > 1 at any time
    min_u_ever: float = 0.0
    max_u_ever: float = 1.0
    oscillatory: bool = False  # Sign changes in incremental updates
    
    # Solution data
    u_final: np.ndarray = field(default_factory=lambda: np.array([]))
    x: np.ndarray = field(default_factory=lambda: np.array([]))
    
    # Early time history for plotting
    u_history: List[np.ndarray] = field(default_factory=list)
    history_steps: List[int] = field(default_factory=list)
    
    # Convergence history
    residual_history: List[float] = field(default_factory=list)
    
    # Runtime
    runtime_seconds: float = 0.0


def initialize_grid(J: int = 32) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Initialize the spatial grid.
    
    Parameters
    ----------
    J : int
        Number of grid intervals. Total points = J + 1.
        Grid indices are j = 0, 1, ..., J.
    
    Returns
    -------
    x : np.ndarray
        Grid point coordinates. x[j] = j * dx.
    u : np.ndarray
        Initial condition array (step profile).
    dx : float
        Grid spacing.
    
    Grid Details
    ------------
    - j = 0 is the left boundary (x = 0)
    - j = J is the right boundary (x = 1)
    - Interior nodes are j = 1, ..., J-1
    - dx = 1.0 / J
    """
    # Grid spacing
    # Domain is [0, 1] with J intervals, so dx = 1/J
    dx = 1.0 / J
    
    # Grid point coordinates
    # x[j] = j * dx for j = 0, 1, ..., J
    x = np.linspace(0.0, 1.0, J + 1)
    
    # Initial condition: step profile
    # u[j] = 1 for j <= J/2
    # u[j] = 0 for j > J/2
    # With J = 32, J/2 = 16, so:
    #   j = 0..16 -> u = 1
    #   j = 17..32 -> u = 0
    u = np.zeros(J + 1)
    midpoint = J // 2  # Integer division gives 16 for J=32
    u[:midpoint + 1] = 1.0  # j = 0, 1, ..., 16
    u[midpoint + 1:] = 0.0  # j = 17, 18, ..., 32
    
    return x, u, dx


def apply_boundary_conditions(u: np.ndarray, u_left: float = 1.0, u_right: float = 0.0):
    """
    Apply Dirichlet boundary conditions in-place.
    
    Parameters
    ----------
    u : np.ndarray
        Solution array (modified in-place).
    u_left : float
        Left boundary value (u[0]).
    u_right : float
        Right boundary value (u[J]).
    """
    u[0] = u_left
    u[-1] = u_right


def explicit_step(u: np.ndarray, mu: float) -> np.ndarray:
    """
    Perform one explicit (FTCS) time step.
    
    For theta = 0, the update is:
        u[j]^(n+1) = u[j]^n + mu * (u[j+1]^n - 2*u[j]^n + u[j-1]^n)
    
    This is applied to interior nodes j = 1, ..., J-1.
    Boundary values are preserved.
    
    Parameters
    ----------
    u : np.ndarray
        Solution at time level n.
    mu : float
        Mesh ratio: alpha * dt / dx^2
    
    Returns
    -------
    u_new : np.ndarray
        Solution at time level n+1.
    
    Notes
    -----
    No matrix solve is needed for the explicit scheme.
    This is a direct update formula.
    
    Stability: requires mu <= 0.5 for stability.
    For mu > 0.5, high-frequency modes will grow exponentially.
    """
    J = len(u) - 1
    u_new = u.copy()
    
    # Update interior nodes: j = 1, 2, ..., J-1
    for j in range(1, J):
        # Second-order central difference for u_xx
        # (u[j+1] - 2*u[j] + u[j-1]) / dx^2
        # Multiplied by alpha * dt = mu * dx^2, giving mu * (...)
        laplacian = u[j + 1] - 2.0 * u[j] + u[j - 1]
        u_new[j] = u[j] + mu * laplacian
    
    # Boundary conditions are preserved (not updated)
    return u_new


def implicit_step(u: np.ndarray, mu: float, lam: float,
                  u_left: float = 1.0, u_right: float = 0.0) -> np.ndarray:
    """
    Perform one implicit or Crank-Nicolson time step using Thomas algorithm.
    
    The discretized system for interior nodes is:
        -lam*mu * u[j-1]^(n+1) + (1 + 2*lam*mu) * u[j]^(n+1) - lam*mu * u[j+1]^(n+1)
            = u[j]^n + (1-lam)*mu * (u[j+1]^n - 2*u[j]^n + u[j-1]^n)
    
    Parameters
    ----------
    u : np.ndarray
        Solution at time level n.
    mu : float
        Mesh ratio.
    lam : float
        Theta parameter. 0.5 = Crank-Nicolson, 1.0 = fully implicit.
    u_left : float
        Left boundary condition.
    u_right : float
        Right boundary condition.
    
    Returns
    -------
    u_new : np.ndarray
        Solution at time level n+1.
    
    Notes
    -----
    Uses the Thomas algorithm (TDMA) to solve the tridiagonal system.
    No full matrices are formed - only 1D arrays for diagonals and RHS.
    
    For Crank-Nicolson (lam=0.5):
        - Second-order accurate in time and space
        - Unconditionally stable (|g| <= 1)
        - Can produce oscillations for large mu (g becomes negative)
    
    For implicit (lam=1.0):
        - First-order accurate in time, second-order in space
        - Unconditionally stable
        - Strongly damping, no oscillations
    """
    J = len(u) - 1
    n_interior = J - 1  # Number of interior unknowns
    
    # Build tridiagonal coefficient arrays
    # a: lower diagonal, b: main diagonal, c: upper diagonal
    a, b, c = build_tridiagonal_coefficients(n_interior, lam, mu)
    
    # Compute right-hand side vector
    rhs = compute_rhs(u, lam, mu, u_left, u_right)
    
    # Solve tridiagonal system using Thomas algorithm
    # This does NOT form any full matrices
    u_interior = thomas_solve(a, b, c, rhs)
    
    # Assemble full solution with boundary conditions
    u_new = np.zeros(J + 1)
    u_new[0] = u_left
    u_new[1:J] = u_interior
    u_new[J] = u_right
    
    return u_new


def compute_residual(u_new: np.ndarray, u_old: np.ndarray) -> float:
    """
    Compute L2 norm (RMS) of the change between time steps.
    
    Uses interior nodes only (j = 1, ..., J-1) for consistency.
    
    Parameters
    ----------
    u_new : np.ndarray
        Solution at time level n+1.
    u_old : np.ndarray
        Solution at time level n.
    
    Returns
    -------
    residual : float
        RMS of (u_new - u_old) over interior nodes.
    """
    # Interior nodes: indices 1 to J-1 (excluding boundaries)
    delta = u_new[1:-1] - u_old[1:-1]
    
    # RMS (Root Mean Square)
    residual = np.sqrt(np.mean(delta ** 2))
    
    return residual


def check_stability(u: np.ndarray, threshold: float = 1e6) -> Tuple[bool, str]:
    """
    Check if solution has become unstable.
    
    Parameters
    ----------
    u : np.ndarray
        Current solution.
    threshold : float
        Maximum allowed magnitude before declaring instability.
    
    Returns
    -------
    stable : bool
        True if solution is stable.
    reason : str
        Reason for instability (empty if stable).
    """
    if np.any(np.isnan(u)):
        return False, "NaN detected"
    if np.any(np.isinf(u)):
        return False, "Inf detected"
    if np.max(np.abs(u)) > threshold:
        return False, f"Solution exceeded threshold {threshold}"
    
    return True, ""


def run_simulation(lam: float, mu: float, J: int = 32, alpha: float = 1.0,
                   tol: float = 1e-5, max_steps: int = 200000,
                   history_steps: List[int] = None) -> SimulationResult:
    """
    Run time-marching simulation to steady state.
    
    Parameters
    ----------
    lam : float
        Theta parameter (0=explicit, 0.5=CN, 1=implicit).
    mu : float
        Mesh ratio.
    J : int
        Number of grid intervals.
    alpha : float
        Diffusion coefficient (default 1.0).
    tol : float
        Convergence tolerance for L2 norm of change.
    max_steps : int
        Maximum number of time steps.
    history_steps : list of int
        Time steps at which to save solution snapshots.
    
    Returns
    -------
    result : SimulationResult
        Container with all simulation results and diagnostics.
    """
    # Default history steps for plotting
    if history_steps is None:
        history_steps = [0, 1, 2, 5, 10, 20, 50, 100]
    
    # Scheme name for reporting
    if lam == 0.0:
        scheme_name = "Explicit (FTCS)"
    elif lam == 0.5:
        scheme_name = "Crank-Nicolson"
    elif lam == 1.0:
        scheme_name = "Implicit (BE)"
    else:
        scheme_name = f"Theta={lam}"
    
    # Initialize grid and solution
    x, u, dx = initialize_grid(J)
    
    # Compute time step from mu
    # mu = alpha * dt / dx^2  =>  dt = mu * dx^2 / alpha
    dt = mu * dx ** 2 / alpha
    
    # Initialize result container
    result = SimulationResult(
        scheme_name=scheme_name,
        lam=lam,
        mu=mu,
        dx=dx,
        dt=dt,
        J=J,
        converged=False,
        n_steps=0,
        final_residual=np.inf,
        stable=True,
        x=x.copy()
    )
    
    # Save initial condition
    if 0 in history_steps:
        result.u_history.append(u.copy())
        result.history_steps.append(0)
    
    # Track min/max for overshoot detection
    min_u = np.min(u)
    max_u = np.max(u)
    result.min_u_ever = min_u
    result.max_u_ever = max_u
    
    # Track sign changes for oscillation detection
    prev_delta_sign = None
    sign_changes = 0
    
    # Boundary conditions
    u_left = 1.0
    u_right = 0.0
    
    # Start timing
    start_time = time_module.time()
    
    # Time marching loop
    for n in range(max_steps):
        # Perform one time step
        if lam == 0.0:
            # Explicit scheme - no matrix solve needed
            u_new = explicit_step(u, mu)
        else:
            # Implicit or Crank-Nicolson - use Thomas algorithm
            u_new = implicit_step(u, mu, lam, u_left, u_right)
        
        # Ensure boundary conditions (should already be correct)
        apply_boundary_conditions(u_new, u_left, u_right)
        
        # Check stability
        stable, reason = check_stability(u_new)
        if not stable:
            result.stable = False
            result.failure_reason = reason
            result.n_steps = n + 1
            result.u_final = u.copy()  # Save last stable solution
            result.runtime_seconds = time_module.time() - start_time
            return result
        
        # Compute residual for convergence check
        residual = compute_residual(u_new, u)
        result.residual_history.append(residual)
        
        # Track min/max for overshoot detection
        current_min = np.min(u_new)
        current_max = np.max(u_new)
        if current_min < result.min_u_ever:
            result.min_u_ever = current_min
        if current_max > result.max_u_ever:
            result.max_u_ever = current_max
        
        # Check for overshoot (solution outside [0, 1])
        if current_min < -1e-10 or current_max > 1.0 + 1e-10:
            result.overshoot = True
        
        # Track oscillation (sign changes in delta at a probe point)
        # Use midpoint as probe
        delta_mid = u_new[J // 2] - u[J // 2]
        current_sign = np.sign(delta_mid) if abs(delta_mid) > 1e-15 else 0
        if prev_delta_sign is not None and current_sign != 0:
            if current_sign != prev_delta_sign and prev_delta_sign != 0:
                sign_changes += 1
        prev_delta_sign = current_sign if current_sign != 0 else prev_delta_sign
        
        # Save history at specified steps
        if (n + 1) in history_steps:
            result.u_history.append(u_new.copy())
            result.history_steps.append(n + 1)
        
        # Update solution
        u = u_new
        
        # Check convergence
        if residual <= tol:
            result.converged = True
            result.n_steps = n + 1
            result.final_residual = residual
            result.u_final = u.copy()
            result.oscillatory = sign_changes > 2
            result.runtime_seconds = time_module.time() - start_time
            return result
    
    # Did not converge within max_steps
    result.n_steps = max_steps
    result.final_residual = residual
    result.u_final = u.copy()
    result.oscillatory = sign_changes > 2
    result.failure_reason = f"Did not converge in {max_steps} steps"
    result.runtime_seconds = time_module.time() - start_time
    
    return result


def run_all_simulations(mu_values: List[float] = None,
                        J: int = 32,
                        tol: float = 1e-5,
                        max_steps: int = 200000) -> List[SimulationResult]:
    """
    Run simulations for all schemes and mu values.
    
    Parameters
    ----------
    mu_values : list of float
        Mesh ratio values to test.
    J : int
        Number of grid intervals.
    tol : float
        Convergence tolerance.
    max_steps : int
        Maximum time steps.
    
    Returns
    -------
    results : list of SimulationResult
        Results for all (scheme, mu) combinations.
    """
    if mu_values is None:
        mu_values = [0.25, 0.5, 1.0, 2.5, 10, 100, 1000]
    
    lam_values = [0.0, 0.5, 1.0]
    
    results = []
    
    print("=" * 70)
    print("PART B: Time-Marching Computations")
    print("=" * 70)
    print(f"\nGrid: J = {J}, {J+1} points, dx = {1.0/J:.6f}")
    print(f"Convergence tolerance: {tol}")
    print(f"Maximum steps: {max_steps}")
    print()
    
    for lam in lam_values:
        for mu in mu_values:
            print(f"Running: theta={lam}, mu={mu}...", end=" ", flush=True)
            
            result = run_simulation(
                lam=lam,
                mu=mu,
                J=J,
                tol=tol,
                max_steps=max_steps
            )
            
            results.append(result)
            
            # Print summary
            if result.stable and result.converged:
                print(f"Converged in {result.n_steps} steps "
                      f"(dt={result.dt:.2e}, time={result.runtime_seconds:.2f}s)")
            elif not result.stable:
                print(f"UNSTABLE: {result.failure_reason}")
            else:
                print(f"Did not converge: residual={result.final_residual:.2e}")
    
    print()
    return results


def generate_results_table(results: List[SimulationResult], 
                           output_dir: str = "outputs") -> str:
    """
    Generate summary table and save as CSV.
    
    Parameters
    ----------
    results : list of SimulationResult
        Simulation results.
    output_dir : str
        Directory to save CSV file.
    
    Returns
    -------
    table_md : str
        Markdown formatted table.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # CSV output
    csv_path = os.path.join(output_dir, "convergence_results.csv")
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'scheme', 'theta', 'mu', 'dt', 'steps', 'status',
            'final_residual', 'runtime_s', 'overshoot', 'oscillatory',
            'min_u', 'max_u'
        ])
        
        for r in results:
            if not r.stable:
                status = f"UNSTABLE: {r.failure_reason}"
            elif r.converged:
                status = "CONVERGED"
            else:
                status = "NOT_CONVERGED"
            
            writer.writerow([
                r.scheme_name, r.lam, r.mu, f"{r.dt:.6e}", r.n_steps, status,
                f"{r.final_residual:.6e}", f"{r.runtime_seconds:.3f}",
                r.overshoot, r.oscillatory,
                f"{r.min_u_ever:.6f}", f"{r.max_u_ever:.6f}"
            ])
    
    print(f"Saved convergence results to: {csv_path}")
    
    # Markdown table
    lines = []
    lines.append("## Convergence Results\n")
    lines.append("| Scheme | theta | mu | dt | Steps | Status | Overshoot | Oscillatory |")
    lines.append("|--------|-------|-----|-----|-------|--------|-----------|-------------|")
    
    for r in results:
        if not r.stable:
            status = "UNSTABLE"
            steps = "-"
        elif r.converged:
            status = "OK"
            steps = str(r.n_steps)
        else:
            status = "NO CONV"
            steps = f">{r.n_steps}"
        
        overshoot = "Yes" if r.overshoot else "No"
        oscillatory = "Yes" if r.oscillatory else "No"
        
        lines.append(f"| {r.scheme_name} | {r.lam} | {r.mu} | {r.dt:.2e} | "
                    f"{steps} | {status} | {overshoot} | {oscillatory} |")
    
    lines.append("")
    return "\n".join(lines)


def plot_solutions(results: List[SimulationResult],
                   output_dir: str = "outputs",
                   show_plots: bool = False) -> List[str]:
    """
    Generate solution plots for all simulations.
    
    Parameters
    ----------
    results : list of SimulationResult
        Simulation results.
    output_dir : str
        Output directory for plots.
    show_plots : bool
        Whether to display plots interactively.
    
    Returns
    -------
    plot_files : list of str
        Paths to generated plot files.
    """
    os.makedirs(output_dir, exist_ok=True)
    plot_files = []
    
    # Group results by scheme
    schemes = {0.0: [], 0.5: [], 1.0: []}
    for r in results:
        schemes[r.lam].append(r)
    
    scheme_names = {0.0: "Explicit", 0.5: "Crank-Nicolson", 1.0: "Implicit"}
    
    # Plot 1: Early time evolution for each scheme (select representative mu values)
    for lam, scheme_results in schemes.items():
        # Find a successful result with small mu and one with large mu
        small_mu_result = None
        large_mu_result = None
        
        for r in scheme_results:
            if r.stable and len(r.u_history) > 1:
                if r.mu <= 0.5 and small_mu_result is None:
                    small_mu_result = r
                elif r.mu >= 2.5 and large_mu_result is None:
                    large_mu_result = r
        
        # Plot early time evolution
        results_to_plot = [r for r in [small_mu_result, large_mu_result] if r is not None]
        
        if results_to_plot:
            fig, axes = plt.subplots(1, len(results_to_plot), figsize=(7*len(results_to_plot), 5))
            if len(results_to_plot) == 1:
                axes = [axes]
            
            for ax, r in zip(axes, results_to_plot):
                for u_snap, step in zip(r.u_history, r.history_steps):
                    alpha_val = 0.3 + 0.7 * (step / max(r.history_steps) if max(r.history_steps) > 0 else 1)
                    ax.plot(r.x, u_snap, label=f'n={step}', alpha=min(1.0, alpha_val))
                
                ax.set_xlabel('x', fontsize=12)
                ax.set_ylabel('u(x)', fontsize=12)
                ax.set_title(f'{scheme_names[lam]}, mu={r.mu}', fontsize=13)
                ax.legend(loc='best', fontsize=9)
                ax.grid(True, alpha=0.3)
                ax.set_xlim([0, 1])
                ax.set_ylim([-0.2, 1.2])
            
            plt.tight_layout()
            filename = os.path.join(output_dir, f"evolution_{scheme_names[lam].lower().replace(' ', '_')}.png")
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            plot_files.append(filename)
            
            if show_plots:
                plt.show()
            plt.close()
    
    # Plot 2: Final solutions comparison for each mu
    mu_values = sorted(set(r.mu for r in results))
    
    for mu in mu_values:  # Generate plots for all mu values
        fig, ax = plt.subplots(figsize=(8, 5))
        
        colors = {0.0: 'red', 0.5: 'blue', 1.0: 'green'}
        
        for lam in [0.0, 0.5, 1.0]:
            matching = [r for r in results if r.lam == lam and r.mu == mu]
            if matching:
                r = matching[0]
                if r.stable and len(r.u_final) > 0:
                    linestyle = '-' if r.converged else '--'
                    ax.plot(r.x, r.u_final, color=colors[lam], linestyle=linestyle,
                           linewidth=2, label=f'{r.scheme_name}')
        
        # Analytical steady state (linear profile)
        x_analytical = np.linspace(0, 1, 100)
        u_analytical = 1.0 - x_analytical
        ax.plot(x_analytical, u_analytical, 'k--', linewidth=1.5, 
               alpha=0.5, label='Steady state')
        
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('u(x)', fontsize=12)
        ax.set_title(f'Final Solutions (mu = {mu})', fontsize=13)
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, 1])
        
        plt.tight_layout()
        filename = os.path.join(output_dir, f"final_solution_mu_{mu:.2f}.png")
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plot_files.append(filename)
        
        if show_plots:
            plt.show()
        plt.close()
    
    # Plot 3: Convergence history
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for idx, lam in enumerate([0.0, 0.5, 1.0]):
        ax = axes[idx]
        
        for r in schemes[lam]:
            if r.stable and len(r.residual_history) > 10:
                # Sample for readability
                n_points = min(1000, len(r.residual_history))
                indices = np.linspace(0, len(r.residual_history)-1, n_points, dtype=int)
                steps = indices + 1
                residuals = [r.residual_history[i] for i in indices]
                
                ax.semilogy(steps, residuals, label=f'mu={r.mu}', alpha=0.8)
        
        ax.axhline(y=1e-5, color='gray', linestyle=':', label='Tolerance')
        ax.set_xlabel('Time step', fontsize=11)
        ax.set_ylabel('L2 residual', fontsize=11)
        ax.set_title(f'{scheme_names[lam]}', fontsize=12)
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filename = os.path.join(output_dir, "convergence_history.png")
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plot_files.append(filename)
    
    if show_plots:
        plt.show()
    plt.close()
    
    print(f"\nGenerated {len(plot_files)} plots in {output_dir}/")
    for pf in plot_files:
        print(f"  {os.path.basename(pf)}")
    
    return plot_files


def generate_computation_report(results: List[SimulationResult]) -> str:
    """
    Generate markdown report section for computations.
    
    Parameters
    ----------
    results : list of SimulationResult
        Simulation results.
    
    Returns
    -------
    report : str
        Markdown formatted report section.
    """
    lines = []
    lines.append("## Time-Marching Computations\n")
    
    lines.append("### Problem Setup\n")
    if results:
        r = results[0]
        lines.append(f"- Domain: [0, 1]\n")
        lines.append(f"- Grid points: J = {r.J}, giving {r.J + 1} total points\n")
        lines.append(f"- Grid spacing: dx = 1/{r.J} = {r.dx:.6f}\n")
        lines.append(f"- Diffusion coefficient: alpha = 1.0\n")
        lines.append(f"- Time step: dt = mu * dx^2 (varies with mu)\n")
    
    lines.append("\n### Boundary Conditions\n")
    lines.append("- Left (x=0): u = 1 (Dirichlet)\n")
    lines.append("- Right (x=1): u = 0 (Dirichlet)\n")
    
    lines.append("\n### Initial Condition\n")
    lines.append("Step profile:\n")
    lines.append("- u(x) = 1 for x <= 0.5\n")
    lines.append("- u(x) = 0 for x > 0.5\n")
    
    lines.append("\n### Key Observations\n\n")
    
    # Analyze results by scheme
    scheme_summary = {0.0: [], 0.5: [], 1.0: []}
    for r in results:
        scheme_summary[r.lam].append(r)
    
    lines.append("#### Explicit Scheme (theta = 0)\n")
    explicit_results = scheme_summary[0.0]
    stable_mu = [r.mu for r in explicit_results if r.stable]
    unstable_mu = [r.mu for r in explicit_results if not r.stable]
    
    if stable_mu:
        lines.append(f"- Stable for mu in {{{', '.join(map(str, sorted(stable_mu)))}}}\n")
    if unstable_mu:
        lines.append(f"- **Unstable** for mu in {{{', '.join(map(str, sorted(unstable_mu)))}}}\n")
    lines.append("- Confirms theoretical stability limit: mu <= 0.5\n")
    lines.append("- Fastest per-step computation (no matrix solve)\n")
    lines.append("- Requires many more steps for large dt\n\n")
    
    lines.append("#### Crank-Nicolson (theta = 0.5)\n")
    cn_results = scheme_summary[0.5]
    cn_overshoot = [r.mu for r in cn_results if r.overshoot and r.stable]
    cn_oscillatory = [r.mu for r in cn_results if r.oscillatory and r.stable]
    
    lines.append("- Unconditionally stable for all tested mu values\n")
    if cn_overshoot:
        lines.append(f"- Shows overshoot for mu in {{{', '.join(map(str, sorted(cn_overshoot)))}}}\n")
    if cn_oscillatory:
        lines.append(f"- Oscillatory behavior for mu in {{{', '.join(map(str, sorted(cn_oscillatory)))}}}\n")
    lines.append("- Second-order accurate, often the best balance of accuracy and efficiency\n\n")
    
    lines.append("#### Implicit Scheme (theta = 1)\n")
    implicit_results = scheme_summary[1.0]
    impl_overshoot = [r.mu for r in implicit_results if r.overshoot and r.stable]
    
    lines.append("- Unconditionally stable for all tested mu values\n")
    lines.append("- No oscillatory behavior (g always positive)\n")
    if not impl_overshoot:
        lines.append("- No overshoot observed - monotonic approach to steady state\n")
    lines.append("- First-order accurate in time, may introduce excessive numerical diffusion\n\n")
    
    lines.append("### Behavior Summary\n\n")
    lines.append("| Behavior | Explicit | Crank-Nicolson | Implicit |\n")
    lines.append("|----------|----------|----------------|----------|\n")
    lines.append("| Stability | mu <= 0.5 | Unconditional | Unconditional |\n")
    lines.append("| Oscillations | No (when stable) | Yes (large mu) | No |\n")
    lines.append("| Accuracy | O(dt, dx^2) | O(dt^2, dx^2) | O(dt, dx^2) |\n")
    lines.append("| Matrix solve | No | Yes (Thomas) | Yes (Thomas) |\n\n")
    
    return "".join(lines)
