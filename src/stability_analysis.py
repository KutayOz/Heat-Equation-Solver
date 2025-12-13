"""
Von Neumann Stability Analysis for 1D Heat Equation
====================================================

This module implements the Von Neumann (Fourier) stability analysis
for the generalized Crank-Nicolson (theta) scheme applied to:

    u_t = alpha * u_xx    (1D heat/diffusion equation)

Discretization (theta scheme):
    (u_i^(n+1) - u_i^n)/dt = alpha * [
        theta * (u_(i+1)^(n+1) - 2*u_i^(n+1) + u_(i-1)^(n+1)) / dx^2
      + (1-theta) * (u_(i+1)^n - 2*u_i^n + u_(i-1)^n) / dx^2
    ]

Key parameters:
    mu = alpha * dt / dx^2   (mesh ratio / CFL number)
    beta = wavenumber * dx   (nondimensional wavenumber, beta in [0, pi])
    theta (lambda) = implicitness parameter

Amplification Factor Derivation
-------------------------------
Assume Fourier mode: u_j^n = G^n * exp(i * j * beta)
where G is the amplification factor per time step.

Substituting into the scheme and simplifying:
    G - 1 = mu * [theta * G * (exp(i*beta) - 2 + exp(-i*beta))
                + (1-theta) * (exp(i*beta) - 2 + exp(-i*beta))]
    
Using: exp(i*beta) - 2 + exp(-i*beta) = 2*cos(beta) - 2 = -4*sin^2(beta/2)

Let s2 = sin^2(beta/2), then:
    G - 1 = -4*mu*s2 * [theta*G + (1-theta)]
    G + 4*mu*theta*s2*G = 1 - 4*mu*(1-theta)*s2
    G * (1 + 4*mu*theta*s2) = 1 - 4*mu*(1-theta)*s2

Therefore:
    g(theta, mu, beta) = (1 - 4*(1-theta)*mu*s2) / (1 + 4*theta*mu*s2)

where s2 = sin^2(beta/2).

Stability criterion: |g| <= 1 for all beta in [0, pi].
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional
import os


def amplification_factor(lam: float, mu: float, beta: np.ndarray) -> np.ndarray:
    """
    Compute the discrete amplification factor g for the theta scheme.
    
    Formula:
        s2 = sin^2(beta/2)
        g = (1 - 4*(1-lam)*mu*s2) / (1 + 4*lam*mu*s2)
    
    Parameters
    ----------
    lam : float
        Theta parameter (lambda).
        0 = explicit (FTCS), 0.5 = Crank-Nicolson, 1 = implicit (backward Euler)
    mu : float
        Mesh ratio: alpha * dt / dx^2
    beta : np.ndarray
        Nondimensional wavenumber array, beta in [0, pi].
    
    Returns
    -------
    g : np.ndarray
        Amplification factor (real-valued for this scheme).
    """
    # s2 = sin^2(beta/2)
    s2 = np.sin(beta / 2.0) ** 2
    
    # Numerator: 1 - 4*(1-lam)*mu*s2
    numerator = 1.0 - 4.0 * (1.0 - lam) * mu * s2
    
    # Denominator: 1 + 4*lam*mu*s2
    denominator = 1.0 + 4.0 * lam * mu * s2
    
    # Amplification factor
    g = numerator / denominator
    
    return g


def pde_amplification_factor(mu: float, beta: np.ndarray) -> np.ndarray:
    """
    Compute the exact (PDE) amplification factor over one time step.
    
    For the continuous heat equation u_t = alpha * u_xx,
    a Fourier mode exp(i*k*x) evolves as exp(-alpha*k^2*t).
    
    Over one time step dt, the amplification is:
        g_pde = exp(-alpha * k^2 * dt)
    
    With k*dx = beta and mu = alpha*dt/dx^2:
        g_pde = exp(-mu * beta^2)
    
    Note: This uses beta^2 directly, which is the correct form for
    the continuous PDE. The discrete scheme uses sin^2(beta/2) which
    approaches beta^2/4 for small beta.
    
    Parameters
    ----------
    mu : float
        Mesh ratio: alpha * dt / dx^2
    beta : np.ndarray
        Nondimensional wavenumber array.
    
    Returns
    -------
    g_pde : np.ndarray
        Exact PDE amplification factor.
    """
    # For the actual PDE, the mode k decays as exp(-alpha*k^2*dt)
    # With beta = k*dx, we have k = beta/dx
    # alpha*k^2*dt = alpha*beta^2*dt/dx^2 = mu*beta^2/dx^2 * dx^2 = mu*(beta/dx*dx)^2
    # Actually: if beta = k*dx, then alpha*k^2*dt = alpha*(beta/dx)^2*dt = mu*beta^2
    # Wait, let me reconsider...
    # mu = alpha*dt/dx^2
    # For mode e^{ikx}, decay is e^{-alpha*k^2*dt}
    # With nondimensional beta = k*dx, k = beta/dx
    # alpha*k^2*dt = alpha*(beta/dx)^2*dt = alpha*beta^2*dt/dx^2 = mu*beta^2
    # But this would make g_pde = exp(-mu*beta^2), which for beta=pi would be
    # exp(-mu*pi^2), a very strong decay.
    # 
    # The standard convention in Von Neumann analysis uses:
    # g_pde = exp(-4*mu*sin^2(beta/2)) for comparison to discrete schemes
    # But the assignment asks for g_pde = exp(-mu*beta^2).
    # 
    # Using the assignment's formula:
    return np.exp(-mu * beta ** 2)


def analyze_stability(lam: float, mu: float, 
                      n_beta: int = 200) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Perform Von Neumann stability analysis for given theta and mu.
    
    Parameters
    ----------
    lam : float
        Theta parameter.
    mu : float
        Mesh ratio.
    n_beta : int
        Number of beta points in [0, pi].
    
    Returns
    -------
    beta : np.ndarray
        Wavenumber array.
    g : np.ndarray
        Amplification factor array.
    g_min : float
        Minimum value of g (can be negative for oscillatory behavior).
    g_max : float
        Maximum value of |g| (for stability check).
    """
    beta = np.linspace(0, np.pi, n_beta)
    g = amplification_factor(lam, mu, beta)
    
    g_min = np.min(g)
    g_max = np.max(np.abs(g))
    
    return beta, g, g_min, g_max


def is_stable(lam: float, mu: float, tol: float = 1e-10) -> Tuple[bool, float]:
    """
    Check if the scheme is stable for given theta and mu.
    
    Stability requires |g| <= 1 for all beta in [0, pi].
    
    Parameters
    ----------
    lam : float
        Theta parameter.
    mu : float
        Mesh ratio.
    tol : float
        Tolerance for stability check.
    
    Returns
    -------
    stable : bool
        True if |g| <= 1 + tol for all beta.
    max_abs_g : float
        Maximum |g| over all beta.
    """
    _, g, _, max_abs_g = analyze_stability(lam, mu, n_beta=1000)
    stable = max_abs_g <= 1.0 + tol
    return stable, max_abs_g


def find_stability_limit_explicit(tol: float = 1e-6) -> float:
    """
    Find the stability limit mu_max for the explicit scheme (theta=0).
    
    For theta=0:
        g = 1 - 4*mu*s2
    
    The worst case is s2 = 1 (beta = pi), giving:
        g = 1 - 4*mu
    
    For stability, we need |g| <= 1:
        -1 <= 1 - 4*mu <= 1
        -2 <= -4*mu <= 0
        0 <= mu <= 0.5
    
    Returns
    -------
    mu_max : float
        Maximum stable mu for explicit scheme.
    """
    # Binary search for stability limit
    mu_low, mu_high = 0.0, 1.0
    
    while mu_high - mu_low > tol:
        mu_mid = (mu_low + mu_high) / 2.0
        stable, _ = is_stable(0.0, mu_mid)
        if stable:
            mu_low = mu_mid
        else:
            mu_high = mu_mid
    
    return mu_low


def plot_amplification_factors(mu_values: List[float],
                               lam_values: List[float] = [0.0, 0.5, 1.0],
                               output_dir: str = "outputs",
                               show_plots: bool = False) -> List[str]:
    """
    Generate plots of amplification factor vs beta for various mu and theta.
    
    Parameters
    ----------
    mu_values : list of float
        Mesh ratio values to plot.
    lam_values : list of float
        Theta values to plot.
    output_dir : str
        Directory to save plots.
    show_plots : bool
        Whether to display plots interactively.
    
    Returns
    -------
    plot_files : list of str
        Paths to generated plot files.
    """
    os.makedirs(output_dir, exist_ok=True)
    plot_files = []
    
    # Color and line style mappings
    scheme_names = {0.0: "Explicit (FTCS)", 0.5: "Crank-Nicolson", 1.0: "Implicit (BE)"}
    colors = {0.0: 'red', 0.5: 'blue', 1.0: 'green'}
    
    beta = np.linspace(0, np.pi, 200)
    
    # Plot 1: g vs beta for each mu (all schemes on same plot)
    for mu in mu_values:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Left plot: g (signed)
        ax1 = axes[0]
        for lam in lam_values:
            g = amplification_factor(lam, mu, beta)
            ax1.plot(beta, g, color=colors[lam], linewidth=2,
                    label=f"{scheme_names[lam]}")
        
        # Plot PDE amplification factor
        g_pde = pde_amplification_factor(mu, beta)
        ax1.plot(beta, g_pde, 'k--', linewidth=1.5, label='PDE exact')
        
        # Stability bounds
        ax1.axhline(y=1, color='gray', linestyle=':', alpha=0.7)
        ax1.axhline(y=-1, color='gray', linestyle=':', alpha=0.7)
        ax1.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
        
        ax1.set_xlabel('beta (wavenumber)', fontsize=12)
        ax1.set_ylabel('g (amplification factor)', fontsize=12)
        ax1.set_title(f'Amplification Factor vs Wavenumber (mu = {mu})', fontsize=13)
        ax1.set_xlim([0, np.pi])
        ax1.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax1.set_xticklabels(['0', 'pi/4', 'pi/2', '3pi/4', 'pi'])
        ax1.legend(loc='best', fontsize=10)
        ax1.grid(True, alpha=0.3)
        
        # Right plot: |g| (magnitude)
        ax2 = axes[1]
        for lam in lam_values:
            g = amplification_factor(lam, mu, beta)
            ax2.plot(beta, np.abs(g), color=colors[lam], linewidth=2,
                    label=f"{scheme_names[lam]}")
        
        # PDE
        ax2.plot(beta, np.abs(g_pde), 'k--', linewidth=1.5, label='PDE exact')
        
        # Stability bound
        ax2.axhline(y=1, color='gray', linestyle=':', alpha=0.7, label='Stability limit')
        
        ax2.set_xlabel('beta (wavenumber)', fontsize=12)
        ax2.set_ylabel('|g| (magnitude)', fontsize=12)
        ax2.set_title(f'Amplification Magnitude vs Wavenumber (mu = {mu})', fontsize=13)
        ax2.set_xlim([0, np.pi])
        ax2.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax2.set_xticklabels(['0', 'pi/4', 'pi/2', '3pi/4', 'pi'])
        ax2.legend(loc='best', fontsize=10)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        filename = os.path.join(output_dir, f"amplification_mu_{mu:.2f}.png")
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plot_files.append(filename)
        
        if show_plots:
            plt.show()
        plt.close()
    
    # Plot 2: Combined plot showing stability regions
    fig, ax = plt.subplots(figsize=(10, 6))
    
    mu_range = np.linspace(0.01, 3.0, 100)
    
    for lam in lam_values:
        max_g_values = []
        min_g_values = []
        for mu in mu_range:
            _, g, g_min, g_max = analyze_stability(lam, mu)
            max_g_values.append(g_max)
            min_g_values.append(g_min)
        
        ax.plot(mu_range, max_g_values, color=colors[lam], linewidth=2,
               label=f"{scheme_names[lam]} - max|g|")
        ax.plot(mu_range, min_g_values, color=colors[lam], linewidth=2,
               linestyle='--', alpha=0.7, label=f"{scheme_names[lam]} - min(g)")
    
    ax.axhline(y=1, color='gray', linestyle=':', alpha=0.7, label='Stability limit')
    ax.axhline(y=-1, color='gray', linestyle=':', alpha=0.7)
    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    ax.axvline(x=0.5, color='red', linestyle=':', alpha=0.5, 
              label='Explicit stability limit (mu=0.5)')
    
    ax.set_xlabel('mu (mesh ratio)', fontsize=12)
    ax.set_ylabel('Amplification factor bounds', fontsize=12)
    ax.set_title('Stability Analysis: Amplification Factor Bounds vs mu', fontsize=13)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 3])
    ax.set_ylim([-2, 2])
    
    plt.tight_layout()
    filename = os.path.join(output_dir, "stability_regions.png")
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plot_files.append(filename)
    
    if show_plots:
        plt.show()
    plt.close()
    
    return plot_files


def generate_stability_report(mu_values: List[float],
                              lam_values: List[float] = [0.0, 0.5, 1.0]) -> str:
    """
    Generate a markdown report section on stability analysis.
    
    Parameters
    ----------
    mu_values : list of float
        Mesh ratio values analyzed.
    lam_values : list of float
        Theta values analyzed.
    
    Returns
    -------
    report : str
        Markdown formatted stability report.
    """
    scheme_names = {0.0: "Explicit (FTCS)", 0.5: "Crank-Nicolson", 1.0: "Implicit (BE)"}
    
    report = []
    report.append("## Von Neumann Stability Analysis\n")
    
    report.append("### Amplification Factor Derivation\n")
    report.append("For the 1D heat equation `u_t = alpha * u_xx` discretized with the ")
    report.append("generalized Crank-Nicolson (theta) scheme:\n\n")
    report.append("```\n")
    report.append("(u_i^(n+1) - u_i^n)/dt = alpha * [\n")
    report.append("    theta * (u_(i+1)^(n+1) - 2*u_i^(n+1) + u_(i-1)^(n+1)) / dx^2\n")
    report.append("  + (1-theta) * (u_(i+1)^n - 2*u_i^n + u_(i-1)^n) / dx^2\n")
    report.append("]\n")
    report.append("```\n\n")
    
    report.append("Using Von Neumann analysis with Fourier mode `u_j^n = G^n * exp(i*j*beta)`,\n")
    report.append("the amplification factor is derived as:\n\n")
    report.append("```\n")
    report.append("s2 = sin^2(beta/2)\n")
    report.append("g(theta, mu, beta) = (1 - 4*(1-theta)*mu*s2) / (1 + 4*theta*mu*s2)\n")
    report.append("```\n\n")
    report.append("where:\n")
    report.append("- `mu = alpha * dt / dx^2` (mesh ratio)\n")
    report.append("- `beta` is the nondimensional wavenumber in `[0, pi]`\n")
    report.append("- `theta` (lambda) is the implicitness parameter\n\n")
    
    report.append("### Stability Results by Scheme\n\n")
    
    # Table header
    report.append("| Scheme | theta | mu | max|g| | min(g) | Stable? | Oscillatory? |\n")
    report.append("|--------|-------|-----|--------|--------|---------|---------------|\n")
    
    for lam in lam_values:
        for mu in mu_values:
            _, g, g_min, g_max = analyze_stability(lam, mu)
            stable = g_max <= 1.0 + 1e-10
            oscillatory = g_min < -1e-10  # g can become negative
            
            report.append(f"| {scheme_names[lam]} | {lam} | {mu} | {g_max:.4f} | ")
            report.append(f"{g_min:.4f} | {'Yes' if stable else 'NO'} | ")
            report.append(f"{'Yes' if oscillatory else 'No'} |\n")
    
    report.append("\n### Stability Discussion\n\n")
    
    report.append("#### Effect of theta (lambda) on Stability\n\n")
    report.append("1. **Explicit scheme (theta = 0.0)**: Conditionally stable.\n")
    report.append("   - The amplification factor simplifies to: `g = 1 - 4*mu*sin^2(beta/2)`\n")
    report.append("   - For stability (`|g| <= 1`), we need `mu <= 0.5`\n")
    report.append("   - The worst case occurs at `beta = pi` (highest frequency mode)\n")
    report.append("   - For `mu > 0.5`, high-frequency modes grow unboundedly\n\n")
    
    report.append("2. **Crank-Nicolson (theta = 0.5)**: Unconditionally stable in magnitude.\n")
    report.append("   - `g = (1 - 2*mu*s2) / (1 + 2*mu*s2)` where `s2 = sin^2(beta/2)`\n")
    report.append("   - For any `mu > 0`, we have `|g| < 1` (strictly damping)\n")
    report.append("   - **However**, for large `mu`, `g` becomes negative for high frequencies\n")
    report.append("   - This causes oscillatory (non-monotonic) behavior in time\n")
    report.append("   - The sign flip occurs when `1 - 2*mu*s2 < 0`, i.e., `mu > 0.5/s2`\n\n")
    
    report.append("3. **Implicit scheme (theta = 1.0)**: Unconditionally stable, strongly damping.\n")
    report.append("   - `g = 1 / (1 + 4*mu*s2)`\n")
    report.append("   - Since denominator `> 1`, we always have `0 < g < 1`\n")
    report.append("   - No oscillatory behavior (g is always positive)\n")
    report.append("   - High frequencies are strongly damped for large `mu`\n\n")
    
    report.append("#### Effect of mu on Stability\n\n")
    report.append("- **Small mu (mu < 0.5)**: All schemes are stable\n")
    report.append("  - Explicit and CN behave similarly with moderate damping\n")
    report.append("  - All schemes preserve solution monotonicity\n\n")
    
    report.append("- **Moderate mu (0.5 < mu < 1)**: Explicit becomes unstable\n")
    report.append("  - CN and implicit remain stable\n")
    report.append("  - CN may show slight oscillations for high-frequency components\n\n")
    
    report.append("- **Large mu (mu >> 1)**: Only implicit methods remain practical\n")
    report.append("  - Explicit: violently unstable (exponential blowup)\n")
    report.append("  - CN: stable but increasingly oscillatory\n")
    report.append("  - Implicit: stable and monotonic, but excessive numerical diffusion\n\n")
    
    report.append("### Key Observations\n\n")
    report.append("- The stability limit `mu <= 0.5` for explicit schemes is exact\n")
    report.append("- CN is often preferred for accuracy (second-order in time and space)\n")
    report.append("- For very large time steps, implicit methods avoid oscillations\n")
    report.append("- The choice of scheme involves trade-offs between:\n")
    report.append("  - Computational cost (explicit is cheapest per step)\n")
    report.append("  - Stability (implicit is most robust)\n")
    report.append("  - Accuracy (CN is most accurate for smooth solutions)\n")
    report.append("  - Solution quality (implicit avoids oscillations)\n\n")
    
    return "".join(report)


def run_stability_analysis(mu_values: List[float] = [0.25, 0.5, 1.0, 2.5],
                           output_dir: str = "outputs",
                           show_plots: bool = False) -> Tuple[List[str], str]:
    """
    Run complete Von Neumann stability analysis and generate outputs.
    
    Parameters
    ----------
    mu_values : list of float
        Mesh ratio values to analyze.
    output_dir : str
        Output directory for plots.
    show_plots : bool
        Whether to display plots interactively.
    
    Returns
    -------
    plot_files : list of str
        Paths to generated plot files.
    report : str
        Markdown stability report.
    """
    print("=" * 60)
    print("PART A: Von Neumann Stability Analysis")
    print("=" * 60)
    
    lam_values = [0.0, 0.5, 1.0]
    scheme_names = {0.0: "Explicit", 0.5: "Crank-Nicolson", 1.0: "Implicit"}
    
    # Print stability summary
    print("\nStability Summary:")
    print("-" * 60)
    print(f"{'Scheme':<20} {'theta':<8} {'mu':<8} {'max|g|':<10} {'Stable?':<10}")
    print("-" * 60)
    
    for lam in lam_values:
        for mu in mu_values:
            _, g, g_min, g_max = analyze_stability(lam, mu)
            stable = "Yes" if g_max <= 1.0 + 1e-10 else "NO"
            print(f"{scheme_names[lam]:<20} {lam:<8.1f} {mu:<8.2f} {g_max:<10.4f} {stable:<10}")
    
    print("-" * 60)
    
    # Find explicit stability limit
    mu_limit = find_stability_limit_explicit()
    print(f"\nExplicit scheme stability limit: mu <= {mu_limit:.4f}")
    
    # Generate plots
    print("\nGenerating amplification factor plots...")
    plot_files = plot_amplification_factors(mu_values, lam_values, output_dir, show_plots)
    
    for pf in plot_files:
        print(f"  Saved: {pf}")
    
    # Generate report
    print("\nGenerating stability report...")
    report = generate_stability_report(mu_values, lam_values)
    
    print("Stability analysis complete.\n")
    
    return plot_files, report
