"""
Output Manager for Heat Equation Solver
=======================================

Manages output folder structure with timestamps for organized results.

Folder Structure:
    outputs/
    ├── stability_analysis_2025-12-13_23-45-30/
    │   ├── amplification_mu_0.25.png
    │   ├── amplification_mu_0.50.png
    │   └── stability_regions.png
    ├── time_marching_2025-12-13_23-46-15/
    │   ├── evolution_explicit.png
    │   ├── evolution_crank-nicolson.png
    │   └── convergence_results.csv
    └── full_run_2025-12-13_23-50-00/
        ├── stability/
        │   └── ...
        ├── computations/
        │   └── ...
        └── report/
            └── report.md
"""

import os
from datetime import datetime
from typing import Tuple, Optional


def get_timestamp() -> str:
    """
    Get current timestamp formatted for folder names.
    
    Returns
    -------
    timestamp : str
        Formatted as YYYY-MM-DD_HH-MM-SS
    """
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def create_run_folder(base_dir: str, run_type: str, 
                      timestamp: Optional[str] = None) -> str:
    """
    Create a timestamped folder for a specific run type.
    
    Parameters
    ----------
    base_dir : str
        Base output directory (e.g., "outputs")
    run_type : str
        Type of run: "stability_analysis", "time_marching", or "full_run"
    timestamp : str, optional
        Custom timestamp. If None, current time is used.
    
    Returns
    -------
    folder_path : str
        Path to the created folder.
    """
    if timestamp is None:
        timestamp = get_timestamp()
    
    folder_name = f"{run_type}_{timestamp}"
    folder_path = os.path.join(base_dir, folder_name)
    
    os.makedirs(folder_path, exist_ok=True)
    
    return folder_path


def create_full_run_structure(base_dir: str, 
                              timestamp: Optional[str] = None) -> Tuple[str, str, str, str]:
    """
    Create complete folder structure for a full run.
    
    Parameters
    ----------
    base_dir : str
        Base output directory.
    timestamp : str, optional
        Custom timestamp.
    
    Returns
    -------
    run_folder : str
        Main run folder path.
    stability_folder : str
        Folder for stability analysis outputs.
    computation_folder : str
        Folder for computation outputs.
    report_folder : str
        Folder for report files.
    """
    if timestamp is None:
        timestamp = get_timestamp()
    
    run_folder = os.path.join(base_dir, f"full_run_{timestamp}")
    stability_folder = os.path.join(run_folder, "stability_analysis")
    computation_folder = os.path.join(run_folder, "time_marching")
    report_folder = os.path.join(run_folder, "report")
    
    os.makedirs(stability_folder, exist_ok=True)
    os.makedirs(computation_folder, exist_ok=True)
    os.makedirs(report_folder, exist_ok=True)
    
    return run_folder, stability_folder, computation_folder, report_folder


def create_stability_folder(base_dir: str, 
                           timestamp: Optional[str] = None) -> str:
    """
    Create folder for stability analysis outputs.
    
    Parameters
    ----------
    base_dir : str
        Base output directory.
    timestamp : str, optional
        Custom timestamp.
    
    Returns
    -------
    folder_path : str
        Path to stability analysis folder.
    """
    return create_run_folder(base_dir, "stability_analysis", timestamp)


def create_computation_folder(base_dir: str,
                             timestamp: Optional[str] = None) -> str:
    """
    Create folder for time-marching computation outputs.
    
    Parameters
    ----------
    base_dir : str
        Base output directory.
    timestamp : str, optional
        Custom timestamp.
    
    Returns
    -------
    folder_path : str
        Path to computation folder.
    """
    return create_run_folder(base_dir, "time_marching", timestamp)


def get_latest_run_folder(base_dir: str, run_type: str) -> Optional[str]:
    """
    Get the most recent run folder of a given type.
    
    Parameters
    ----------
    base_dir : str
        Base output directory.
    run_type : str
        Type of run to search for.
    
    Returns
    -------
    folder_path : str or None
        Path to the latest folder, or None if not found.
    """
    if not os.path.exists(base_dir):
        return None
    
    matching_folders = []
    for item in os.listdir(base_dir):
        if item.startswith(run_type) and os.path.isdir(os.path.join(base_dir, item)):
            matching_folders.append(item)
    
    if not matching_folders:
        return None
    
    # Sort by name (timestamp in name makes this work)
    matching_folders.sort(reverse=True)
    
    return os.path.join(base_dir, matching_folders[0])


def print_output_summary(output_folder: str):
    """
    Print a summary of files in the output folder.
    
    Parameters
    ----------
    output_folder : str
        Path to output folder.
    """
    print(f"\nOutput folder: {output_folder}")
    print("-" * 50)
    
    for root, dirs, files in os.walk(output_folder):
        level = root.replace(output_folder, '').count(os.sep)
        indent = '  ' * level
        folder_name = os.path.basename(root)
        
        if level == 0:
            print(f"{indent}{folder_name}/")
        else:
            print(f"{indent}{folder_name}/")
        
        sub_indent = '  ' * (level + 1)
        for file in sorted(files):
            file_path = os.path.join(root, file)
            size_kb = os.path.getsize(file_path) / 1024
            print(f"{sub_indent}{file} ({size_kb:.1f} KB)")
