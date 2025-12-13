"""
Report Generator for Heat Equation Solver
==========================================

Automatically generates a markdown report combining results from:
- Part A: Von Neumann stability analysis
- Part B: Time-marching computations

The report includes:
- Theoretical derivations and explanations
- Generated plots (referenced by path)
- Summary tables
- Key observations and discussion
"""

import os
from typing import List, Optional
from datetime import datetime


def generate_full_report(stability_report: str,
                        computation_report: str,
                        results_table: str,
                        stability_plots: List[str],
                        computation_plots: List[str],
                        output_dir: str = "report") -> str:
    """
    Generate the complete markdown report.
    
    Parameters
    ----------
    stability_report : str
        Markdown content for Part A (stability analysis).
    computation_report : str
        Markdown content for Part B (time marching).
    results_table : str
        Markdown table of convergence results.
    stability_plots : list of str
        Paths to stability analysis plots.
    computation_plots : list of str
        Paths to computation plots.
    output_dir : str
        Directory to save the report.
    
    Returns
    -------
    report_path : str
        Path to the generated report file.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    report_lines = []
    
    # Title and header
    report_lines.append("# 1D Heat Equation Solver Report\n")
    report_lines.append("## Generalized Crank-Nicolson (Theta) Scheme\n")
    report_lines.append(f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
    
    report_lines.append("---\n\n")
    
    # Table of contents
    report_lines.append("## Table of Contents\n")
    report_lines.append("1. [Introduction](#introduction)\n")
    report_lines.append("2. [Part A: Von Neumann Stability Analysis](#von-neumann-stability-analysis)\n")
    report_lines.append("3. [Part B: Time-Marching Computations](#time-marching-computations)\n")
    report_lines.append("4. [Convergence Results](#convergence-results)\n")
    report_lines.append("5. [Conclusions](#conclusions)\n\n")
    
    report_lines.append("---\n\n")
    
    # Introduction
    report_lines.append("## Introduction\n\n")
    report_lines.append("This report presents the analysis and numerical solution of the ")
    report_lines.append("1D heat (diffusion) equation using the generalized Crank-Nicolson ")
    report_lines.append("(theta) scheme.\n\n")
    
    report_lines.append("### The Heat Equation\n\n")
    report_lines.append("We solve the partial differential equation:\n\n")
    report_lines.append("```\n")
    report_lines.append("u_t = alpha * u_xx\n")
    report_lines.append("```\n\n")
    report_lines.append("where:\n")
    report_lines.append("- `u(x,t)` is the temperature distribution\n")
    report_lines.append("- `alpha` is the thermal diffusivity\n")
    report_lines.append("- Subscripts denote partial derivatives\n\n")
    
    report_lines.append("### The Theta Scheme\n\n")
    report_lines.append("The generalized Crank-Nicolson discretization is:\n\n")
    report_lines.append("```\n")
    report_lines.append("(u_i^(n+1) - u_i^n)/dt = alpha * [\n")
    report_lines.append("    theta * (u_(i+1)^(n+1) - 2*u_i^(n+1) + u_(i-1)^(n+1)) / dx^2\n")
    report_lines.append("  + (1-theta) * (u_(i+1)^n - 2*u_i^n + u_(i-1)^n) / dx^2\n")
    report_lines.append("]\n")
    report_lines.append("```\n\n")
    
    report_lines.append("Special cases:\n")
    report_lines.append("- `theta = 0`: Explicit (FTCS) - Forward Time, Central Space\n")
    report_lines.append("- `theta = 0.5`: Crank-Nicolson - second-order accurate\n")
    report_lines.append("- `theta = 1`: Implicit (Backward Euler)\n\n")
    
    report_lines.append("Key parameter: `mu = alpha * dt / dx^2` (mesh ratio / CFL number)\n\n")
    
    report_lines.append("---\n\n")
    
    # Part A: Stability Analysis
    report_lines.append(stability_report)
    report_lines.append("\n")
    
    # Include stability plots
    report_lines.append("### Stability Plots\n\n")
    for plot_path in stability_plots:
        rel_path = os.path.relpath(plot_path, output_dir)
        plot_name = os.path.basename(plot_path).replace('.png', '').replace('_', ' ').title()
        report_lines.append(f"#### {plot_name}\n")
        report_lines.append(f"![{plot_name}]({rel_path})\n\n")
    
    report_lines.append("---\n\n")
    
    # Part B: Computations
    report_lines.append(computation_report)
    report_lines.append("\n")
    
    # Include computation plots
    report_lines.append("### Computation Plots\n\n")
    for plot_path in computation_plots:
        rel_path = os.path.relpath(plot_path, output_dir)
        plot_name = os.path.basename(plot_path).replace('.png', '').replace('_', ' ').title()
        report_lines.append(f"#### {plot_name}\n")
        report_lines.append(f"![{plot_name}]({rel_path})\n\n")
    
    report_lines.append("---\n\n")
    
    # Results table
    report_lines.append(results_table)
    report_lines.append("\n")
    
    report_lines.append("---\n\n")
    
    # Conclusions
    report_lines.append("## Conclusions\n\n")
    
    report_lines.append("### Summary of Findings\n\n")
    
    report_lines.append("1. **Stability Analysis Confirms Theory**\n")
    report_lines.append("   - Explicit scheme is conditionally stable with limit `mu <= 0.5`\n")
    report_lines.append("   - Crank-Nicolson and Implicit are unconditionally stable\n")
    report_lines.append("   - CN can exhibit oscillations for large mu despite stability\n\n")
    
    report_lines.append("2. **Computational Efficiency Trade-offs**\n")
    report_lines.append("   - Explicit: cheapest per step, but limited by stability\n")
    report_lines.append("   - Implicit: allows large time steps, but first-order accurate\n")
    report_lines.append("   - CN: best accuracy per computational cost for smooth solutions\n\n")
    
    report_lines.append("3. **Solution Quality Considerations**\n")
    report_lines.append("   - CN oscillations are bounded but may be undesirable\n")
    report_lines.append("   - Implicit scheme provides monotonic convergence\n")
    report_lines.append("   - For stiff problems (large mu), implicit methods are preferred\n\n")
    
    report_lines.append("### Recommendations\n\n")
    report_lines.append("- For high accuracy: use Crank-Nicolson with moderate mu\n")
    report_lines.append("- For robustness: use Implicit with any mu\n")
    report_lines.append("- For efficiency with accuracy: use Explicit with mu near 0.5\n")
    report_lines.append("- Avoid CN with very large mu if oscillations are problematic\n\n")
    
    report_lines.append("---\n\n")
    report_lines.append("*End of Report*\n")
    
    # Write report
    report_content = "".join(report_lines)
    report_path = os.path.join(output_dir, "report.md")
    
    with open(report_path, 'w') as f:
        f.write(report_content)
    
    print(f"Report saved to: {report_path}")
    
    return report_path
