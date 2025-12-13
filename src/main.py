"""
Main Entry Point for Heat Equation Solver
==========================================

Command-line interface for running:
- Part A: Von Neumann stability analysis
- Part B: Time-marching computations

Usage:
    python -m src.main --run_all
    python -m src.main --part analysis_only
    python -m src.main --part compute_only
    python -m src.main --mu_list 0.25,0.5,1.0,2.5
    python -m src.main --max_steps 100000 --tol 1e-6
"""

import argparse
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


def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns
    -------
    args : argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="1D Heat Equation Solver with Crank-Nicolson Theta Scheme",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m src.main --run_all
    python -m src.main --part analysis_only
    python -m src.main --part compute_only
    python -m src.main --mu_list 0.25,0.5,1.0 --max_steps 50000
        """
    )
    
    parser.add_argument(
        '--run_all',
        action='store_true',
        help='Run both stability analysis and computations (default behavior)'
    )
    
    parser.add_argument(
        '--part',
        type=str,
        choices=['analysis_only', 'compute_only', 'both'],
        default='both',
        help='Which part to run: analysis_only, compute_only, or both'
    )
    
    parser.add_argument(
        '--mu_list',
        type=str,
        default='0.25,0.5,1.0,2.5,10,100,1000',
        help='Comma-separated list of mu values for computations'
    )
    
    parser.add_argument(
        '--mu_analysis',
        type=str,
        default='0.25,0.5,1.0,2.5',
        help='Comma-separated list of mu values for stability analysis plots'
    )
    
    parser.add_argument(
        '--max_steps',
        type=int,
        default=200000,
        help='Maximum number of time steps for convergence'
    )
    
    parser.add_argument(
        '--tol',
        type=float,
        default=1e-5,
        help='Convergence tolerance (L2 norm of change)'
    )
    
    parser.add_argument(
        '--J',
        type=int,
        default=32,
        help='Number of grid intervals (J+1 total points)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='outputs',
        help='Directory for output files'
    )
    
    parser.add_argument(
        '--report_dir',
        type=str,
        default='report',
        help='Directory for report files'
    )
    
    parser.add_argument(
        '--show_plots',
        action='store_true',
        help='Display plots interactively (requires display)'
    )
    
    return parser.parse_args()


def main():
    """
    Main entry point for the heat equation solver.
    
    Runs stability analysis and/or time-marching computations
    based on command-line arguments.
    """
    args = parse_arguments()
    
    # Parse mu values
    mu_compute = [float(x.strip()) for x in args.mu_list.split(',')]
    mu_analysis = [float(x.strip()) for x in args.mu_analysis.split(',')]
    
    # Determine what to run
    run_analysis = args.run_all or args.part in ['analysis_only', 'both']
    run_compute = args.run_all or args.part in ['compute_only', 'both']
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.report_dir, exist_ok=True)
    
    print("=" * 70)
    print("1D HEAT EQUATION SOLVER")
    print("Generalized Crank-Nicolson (Theta) Scheme")
    print("=" * 70)
    print()
    
    # Initialize report components
    stability_report = ""
    computation_report = ""
    results_table = ""
    stability_plots = []
    computation_plots = []
    
    # Part A: Stability Analysis
    if run_analysis:
        stability_plots, stability_report = run_stability_analysis(
            mu_values=mu_analysis,
            output_dir=args.output_dir,
            show_plots=args.show_plots
        )
    
    # Part B: Time-Marching Computations
    if run_compute:
        # Run simulations
        results = run_all_simulations(
            mu_values=mu_compute,
            J=args.J,
            tol=args.tol,
            max_steps=args.max_steps
        )
        
        # Generate results table and CSV
        results_table = generate_results_table(results, args.output_dir)
        
        # Generate solution plots
        computation_plots = plot_solutions(
            results,
            args.output_dir,
            args.show_plots
        )
        
        # Generate computation report section
        computation_report = generate_computation_report(results)
    
    # Generate full report
    if run_analysis or run_compute:
        report_path = generate_full_report(
            stability_report=stability_report,
            computation_report=computation_report,
            results_table=results_table,
            stability_plots=stability_plots,
            computation_plots=computation_plots,
            output_dir=args.report_dir
        )
        
        print()
        print("=" * 70)
        print("COMPLETE")
        print("=" * 70)
        print(f"Output files: {args.output_dir}/")
        print(f"Report: {report_path}")
        print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
