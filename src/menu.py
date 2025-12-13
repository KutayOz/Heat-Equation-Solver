"""
Interactive Menu for Heat Equation Solver
==========================================

A user-friendly menu interface that doesn't require command-line arguments.
Simply run: python3 -m src.menu

This provides the same functionality as main.py but with interactive prompts.
"""

import sys
import os

# Add parent directory to path for module imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.stability_analysis import run_stability_analysis
from src.time_marching import (
    run_all_simulations,
    generate_results_table,
    plot_solutions,
    generate_computation_report
)
from src.report_generator import generate_full_report


def clear_screen():
    """Clear the terminal screen."""
    os.system('cls' if os.name == 'nt' else 'clear')


def print_header():
    """Print the program header."""
    print("=" * 60)
    print("       1D HEAT EQUATION SOLVER")
    print("   Generalized Crank-Nicolson (Theta) Scheme")
    print("=" * 60)
    print()


def print_menu():
    """Print the main menu options."""
    print("Main Menu:")
    print("-" * 40)
    print("  1. Run ALL (stability analysis + computations)")
    print("  2. Run stability analysis only (Part A)")
    print("  3. Run time-marching computations only (Part B)")
    print("  4. Custom run with parameters")
    print("  5. Quick info about the solver")
    print("  0. Exit")
    print("-" * 40)


def get_choice(prompt: str, valid_choices: list) -> str:
    """Get a valid choice from the user."""
    while True:
        choice = input(prompt).strip()
        if choice in valid_choices:
            return choice
        print(f"Invalid choice. Please enter one of: {', '.join(valid_choices)}")


def get_float_list(prompt: str, default: str) -> list:
    """Get a comma-separated list of floats from the user."""
    print(f"{prompt}")
    print(f"  (Press Enter for default: {default})")
    user_input = input("  > ").strip()
    
    if not user_input:
        user_input = default
    
    try:
        values = [float(x.strip()) for x in user_input.split(',')]
        return values
    except ValueError:
        print("Invalid input. Using default values.")
        return [float(x.strip()) for x in default.split(',')]


def get_int(prompt: str, default: int) -> int:
    """Get an integer from the user."""
    print(f"{prompt}")
    print(f"  (Press Enter for default: {default})")
    user_input = input("  > ").strip()
    
    if not user_input:
        return default
    
    try:
        return int(user_input)
    except ValueError:
        print(f"Invalid input. Using default: {default}")
        return default


def get_float(prompt: str, default: float) -> float:
    """Get a float from the user."""
    print(f"{prompt}")
    print(f"  (Press Enter for default: {default})")
    user_input = input("  > ").strip()
    
    if not user_input:
        return default
    
    try:
        return float(user_input)
    except ValueError:
        print(f"Invalid input. Using default: {default}")
        return default


def get_yes_no(prompt: str, default: bool = False) -> bool:
    """Get a yes/no answer from the user."""
    default_str = "Y/n" if default else "y/N"
    user_input = input(f"{prompt} [{default_str}]: ").strip().lower()
    
    if not user_input:
        return default
    return user_input in ['y', 'yes']


def run_full(output_dir: str = "outputs", report_dir: str = "report",
             mu_analysis: list = None, mu_compute: list = None,
             J: int = 32, tol: float = 1e-5, max_steps: int = 200000):
    """Run both stability analysis and computations."""
    if mu_analysis is None:
        mu_analysis = [0.25, 0.5, 1.0, 2.5]
    if mu_compute is None:
        mu_compute = [0.25, 0.5, 1.0, 2.5, 10, 100, 1000]
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)
    
    # Part A
    print("\n" + "=" * 60)
    stability_plots, stability_report = run_stability_analysis(
        mu_values=mu_analysis,
        output_dir=output_dir,
        show_plots=False
    )
    
    # Part B
    results = run_all_simulations(
        mu_values=mu_compute,
        J=J,
        tol=tol,
        max_steps=max_steps
    )
    
    results_table = generate_results_table(results, output_dir)
    computation_plots = plot_solutions(results, output_dir, show_plots=False)
    computation_report = generate_computation_report(results)
    
    # Generate report
    report_path = generate_full_report(
        stability_report=stability_report,
        computation_report=computation_report,
        results_table=results_table,
        stability_plots=stability_plots,
        computation_plots=computation_plots,
        output_dir=report_dir
    )
    
    print("\n" + "=" * 60)
    print("COMPLETE!")
    print("=" * 60)
    print(f"Output files saved to: {output_dir}/")
    print(f"Report saved to: {report_path}")


def run_analysis_only(output_dir: str = "outputs", 
                      mu_analysis: list = None):
    """Run stability analysis only."""
    if mu_analysis is None:
        mu_analysis = [0.25, 0.5, 1.0, 2.5]
    
    os.makedirs(output_dir, exist_ok=True)
    
    stability_plots, stability_report = run_stability_analysis(
        mu_values=mu_analysis,
        output_dir=output_dir,
        show_plots=False
    )
    
    print("\n" + "=" * 60)
    print("STABILITY ANALYSIS COMPLETE!")
    print("=" * 60)
    print(f"Plots saved to: {output_dir}/")


def run_compute_only(output_dir: str = "outputs", report_dir: str = "report",
                     mu_compute: list = None, J: int = 32, 
                     tol: float = 1e-5, max_steps: int = 200000):
    """Run time-marching computations only."""
    if mu_compute is None:
        mu_compute = [0.25, 0.5, 1.0, 2.5, 10, 100, 1000]
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)
    
    results = run_all_simulations(
        mu_values=mu_compute,
        J=J,
        tol=tol,
        max_steps=max_steps
    )
    
    results_table = generate_results_table(results, output_dir)
    computation_plots = plot_solutions(results, output_dir, show_plots=False)
    
    print("\n" + "=" * 60)
    print("COMPUTATIONS COMPLETE!")
    print("=" * 60)
    print(f"Results saved to: {output_dir}/")


def custom_run_menu():
    """Interactive menu for custom parameters."""
    print("\n" + "-" * 40)
    print("Custom Run Configuration")
    print("-" * 40)
    
    # What to run
    print("\nWhat would you like to run?")
    print("  1. Both (stability + computations)")
    print("  2. Stability analysis only")
    print("  3. Computations only")
    run_choice = get_choice("Choice [1/2/3]: ", ['1', '2', '3'])
    
    # Get parameters based on choice
    mu_analysis = None
    mu_compute = None
    J = 32
    tol = 1e-5
    max_steps = 200000
    
    if run_choice in ['1', '2']:
        print("\n--- Stability Analysis Parameters ---")
        mu_analysis = get_float_list(
            "mu values for stability plots:",
            "0.25,0.5,1.0,2.5"
        )
    
    if run_choice in ['1', '3']:
        print("\n--- Computation Parameters ---")
        mu_compute = get_float_list(
            "mu values for computations:",
            "0.25,0.5,1.0,2.5,10,100,1000"
        )
        J = get_int("Grid intervals (J):", 32)
        tol = get_float("Convergence tolerance:", 1e-5)
        max_steps = get_int("Maximum time steps:", 200000)
    
    # Output directories
    print("\n--- Output Directories ---")
    output_dir = input("Output directory (Enter for 'outputs'): ").strip()
    if not output_dir:
        output_dir = "outputs"
    
    report_dir = input("Report directory (Enter for 'report'): ").strip()
    if not report_dir:
        report_dir = "report"
    
    # Confirm and run
    print("\n" + "-" * 40)
    print("Configuration Summary:")
    if mu_analysis:
        print(f"  Stability mu values: {mu_analysis}")
    if mu_compute:
        print(f"  Computation mu values: {mu_compute}")
        print(f"  Grid intervals: {J}")
        print(f"  Tolerance: {tol}")
        print(f"  Max steps: {max_steps}")
    print(f"  Output dir: {output_dir}")
    print(f"  Report dir: {report_dir}")
    print("-" * 40)
    
    if get_yes_no("Proceed with this configuration?", default=True):
        if run_choice == '1':
            run_full(output_dir, report_dir, mu_analysis, mu_compute, J, tol, max_steps)
        elif run_choice == '2':
            run_analysis_only(output_dir, mu_analysis)
        else:
            run_compute_only(output_dir, report_dir, mu_compute, J, tol, max_steps)
    else:
        print("Cancelled.")


def show_info():
    """Display information about the solver."""
    print("\n" + "=" * 60)
    print("ABOUT THIS SOLVER")
    print("=" * 60)
    print("""
This program solves the 1D heat (diffusion) equation:

    u_t = alpha * u_xx

using the generalized Crank-Nicolson (theta) scheme.

SCHEMES AVAILABLE:
  - Explicit (theta=0): Fast but conditionally stable (mu <= 0.5)
  - Crank-Nicolson (theta=0.5): Second-order accurate, unconditionally stable
  - Implicit (theta=1): First-order, unconditionally stable, no oscillations

KEY PARAMETER:
  mu = alpha * dt / dx^2  (mesh ratio)

PART A - STABILITY ANALYSIS:
  Von Neumann (Fourier) analysis of the amplification factor:
  g = (1 - 4*(1-theta)*mu*s2) / (1 + 4*theta*mu*s2)
  where s2 = sin^2(beta/2)

PART B - TIME MARCHING:
  Solve to steady state on [0,1] with:
  - Left BC: u(0) = 1
  - Right BC: u(1) = 0
  - Initial: step function at x = 0.5

TRIDIAGONAL SOLVER:
  Uses Thomas algorithm (TDMA) - O(n) time, no full matrices.
""")
    input("\nPress Enter to return to menu...")


def main():
    """Main menu loop."""
    while True:
        clear_screen()
        print_header()
        print_menu()
        
        choice = get_choice("\nEnter your choice: ", ['0', '1', '2', '3', '4', '5'])
        
        if choice == '0':
            print("\nGoodbye!")
            break
        elif choice == '1':
            print("\nRunning full analysis with default parameters...")
            print("(mu_analysis: 0.25, 0.5, 1.0, 2.5)")
            print("(mu_compute: 0.25, 0.5, 1.0, 2.5, 10, 100, 1000)")
            print("(J=32, tol=1e-5, max_steps=200000)")
            if get_yes_no("\nProceed?", default=True):
                run_full()
            input("\nPress Enter to continue...")
        elif choice == '2':
            print("\nRunning stability analysis only...")
            if get_yes_no("Proceed with default mu values (0.25, 0.5, 1.0, 2.5)?", default=True):
                run_analysis_only()
            input("\nPress Enter to continue...")
        elif choice == '3':
            print("\nRunning time-marching computations only...")
            if get_yes_no("Proceed with default parameters?", default=True):
                run_compute_only()
            input("\nPress Enter to continue...")
        elif choice == '4':
            custom_run_menu()
            input("\nPress Enter to continue...")
        elif choice == '5':
            show_info()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
