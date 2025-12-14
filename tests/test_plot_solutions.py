"""
Integration test for verifying PNG generation for all mu values including mu=10, 100, 1000.

This test runs the main program with a subset of mu values and verifies output files.
"""
import os
import sys
import subprocess
import tempfile
import glob


def test_plot_solutions_generates_all_mu_pngs():
    """
    Test that running the heat equation solver generates PNG files for ALL mu values,
    including mu=10, mu=100, and mu=1000.
    """
    mu_values = [0.25, 0.5, 1.0, 2.5, 10, 100, 1000]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run the main program with compute_only to generate the plots
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        venv_python = os.path.join(project_root, 'venv311', 'bin', 'python')
        
        cmd = [
            venv_python, '-m', 'src.main',
            '--part', 'compute_only',
            '--mu_list', ','.join(str(m) for m in mu_values),
            '--output_dir', tmpdir,
            '--max_steps', '1000',  # Low steps for faster test
            '--J', '8'  # Smaller grid for faster test
        ]
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=project_root, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"Main program failed with return code {result.returncode}")
        
        # Check that PNG files exist for ALL mu values
        missing = []
        found = []
        for mu in mu_values:
            expected_file = os.path.join(tmpdir, f"final_solution_mu_{mu:.2f}.png")
            if os.path.exists(expected_file):
                found.append(mu)
            else:
                missing.append(mu)
        
        # List all generated files
        all_pngs = glob.glob(os.path.join(tmpdir, "*.png"))
        print(f"\nGenerated PNG files ({len(all_pngs)}):")
        for f in sorted(all_pngs):
            print(f"  - {os.path.basename(f)}")
        
        if missing:
            print(f"\nMISSING PNGs for mu values: {missing}")
            raise AssertionError(f"Missing PNG files for mu values: {missing}")
        
        print(f"\nSUCCESS: All {len(mu_values)} mu value PNGs generated correctly")
        print(f"Found PNGs for mu: {found}")
    
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("Integration Test: PNG generation for all mu values")
    print("=" * 60)
    print()
    
    try:
        test_plot_solutions_generates_all_mu_pngs()
        print()
        print("=" * 60)
        print("ALL TESTS PASSED")
        print("=" * 60)
    except AssertionError as e:
        print(f"\nTEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
