"""
Thomas Algorithm (Tridiagonal Matrix Algorithm - TDMA)
======================================================

Solves a tridiagonal system of linear equations:
    a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1] = d[i]

where:
    - a is the lower diagonal (a[0] is not used)
    - b is the main diagonal
    - c is the upper diagonal (c[n-1] is not used)
    - d is the right-hand side vector

This implementation does NOT form any full NxN matrices.
All operations use 1D arrays only.

The Thomas algorithm consists of two phases:
1. Forward elimination: eliminate the lower diagonal
2. Back substitution: solve for unknowns from bottom to top
"""

import numpy as np


def thomas_solve(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """
    Solve a tridiagonal system using the Thomas algorithm.
    
    The system is:
        a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1] = d[i]
    
    Parameters
    ----------
    a : np.ndarray
        Lower diagonal coefficients. a[0] is not used.
        Length n.
    b : np.ndarray
        Main diagonal coefficients. Length n.
    c : np.ndarray
        Upper diagonal coefficients. c[n-1] is not used.
        Length n.
    d : np.ndarray
        Right-hand side vector. Length n.
    
    Returns
    -------
    x : np.ndarray
        Solution vector. Length n.
    
    Notes
    -----
    - This function modifies copies of the input arrays to avoid side effects.
    - The algorithm is O(n) in both time and space.
    - No full matrices are formed; only 1D arrays are used.
    
    Algorithm Details
    -----------------
    Forward Elimination (i = 1 to n-1):
        w = a[i] / b[i-1]
        b[i] = b[i] - w * c[i-1]
        d[i] = d[i] - w * d[i-1]
    
    Back Substitution (i = n-2 down to 0):
        x[n-1] = d[n-1] / b[n-1]
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    """
    n = len(d)
    
    # Create copies to avoid modifying input arrays
    # These are 1D arrays only - no matrices
    b_mod = b.copy()
    d_mod = d.copy()
    
    # Allocate solution array
    x = np.zeros(n)
    
    # Phase 1: Forward Elimination
    # ----------------------------
    # Eliminate the lower diagonal by subtracting a scaled version
    # of the previous equation from the current equation.
    # After this phase, the system becomes upper bidiagonal.
    for i in range(1, n):
        # Compute the multiplier for elimination
        # w = a[i] / b[i-1] is the factor needed to eliminate a[i]
        w = a[i] / b_mod[i - 1]
        
        # Update the main diagonal: b[i] = b[i] - w * c[i-1]
        b_mod[i] = b_mod[i] - w * c[i - 1]
        
        # Update the RHS: d[i] = d[i] - w * d[i-1]
        d_mod[i] = d_mod[i] - w * d_mod[i - 1]
    
    # Phase 2: Back Substitution
    # --------------------------
    # Starting from the last equation (which now has only one unknown),
    # solve for each unknown moving backward.
    
    # Last unknown: x[n-1] = d[n-1] / b[n-1]
    x[n - 1] = d_mod[n - 1] / b_mod[n - 1]
    
    # Remaining unknowns: x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    for i in range(n - 2, -1, -1):
        x[i] = (d_mod[i] - c[i] * x[i + 1]) / b_mod[i]
    
    return x


def build_tridiagonal_coefficients(n: int, lam: float, mu: float) -> tuple:
    """
    Build the tridiagonal coefficient arrays for the generalized
    Crank-Nicolson scheme.
    
    The discretized system for interior nodes j = 1, ..., J-1 is:
    
    -lam*mu * u[j-1]^(n+1) + (1 + 2*lam*mu) * u[j]^(n+1) - lam*mu * u[j+1]^(n+1)
        = RHS (computed separately)
    
    Parameters
    ----------
    n : int
        Number of interior unknowns (J-1 for grid with J+1 total points).
    lam : float
        Theta parameter (lambda). 
        0 = explicit, 0.5 = Crank-Nicolson, 1 = implicit.
    mu : float
        Mesh ratio: alpha * dt / dx^2
    
    Returns
    -------
    a : np.ndarray
        Lower diagonal: -lam * mu (except a[0] unused)
    b : np.ndarray
        Main diagonal: 1 + 2 * lam * mu
    c : np.ndarray
        Upper diagonal: -lam * mu (except c[n-1] unused)
    
    Notes
    -----
    For the explicit scheme (lam=0), the LHS matrix is identity,
    so no Thomas solve is needed. This function is still valid
    but the solve becomes trivial.
    """
    # Lower diagonal coefficient
    # a[i] multiplies u[j-1]^(n+1) in equation for u[j]
    a = np.full(n, -lam * mu)
    
    # Main diagonal coefficient
    # b[i] multiplies u[j]^(n+1)
    b = np.full(n, 1.0 + 2.0 * lam * mu)
    
    # Upper diagonal coefficient
    # c[i] multiplies u[j+1]^(n+1)
    c = np.full(n, -lam * mu)
    
    return a, b, c


def compute_rhs(u: np.ndarray, lam: float, mu: float, 
                u_left: float, u_right: float) -> np.ndarray:
    """
    Compute the right-hand side vector for the tridiagonal system.
    
    The RHS for interior node j is:
        RHS[j] = u[j]^n + (1-lam)*mu * (u[j+1]^n - 2*u[j]^n + u[j-1]^n)
    
    Plus boundary contributions:
        - For j=1 (first interior): add lam*mu * u_left to RHS
        - For j=J-1 (last interior): add lam*mu * u_right to RHS
    
    Parameters
    ----------
    u : np.ndarray
        Full solution array at time level n, including boundaries.
        u[0] = left BC, u[J] = right BC, u[1..J-1] = interior.
    lam : float
        Theta parameter.
    mu : float
        Mesh ratio.
    u_left : float
        Left boundary value at time n+1 (Dirichlet BC).
    u_right : float
        Right boundary value at time n+1 (Dirichlet BC).
    
    Returns
    -------
    rhs : np.ndarray
        Right-hand side vector for interior nodes.
        Length = len(u) - 2 (number of interior nodes).
    """
    J = len(u) - 1  # u has J+1 points, indices 0..J
    n_interior = J - 1  # Interior nodes are j = 1..J-1
    
    rhs = np.zeros(n_interior)
    
    # Coefficient for explicit part
    coef_explicit = (1.0 - lam) * mu
    
    # Compute RHS for each interior node
    for i in range(n_interior):
        j = i + 1  # Map array index i to grid index j
        
        # Base term: u[j]^n
        rhs[i] = u[j]
        
        # Explicit part: (1-lam)*mu * (u[j+1]^n - 2*u[j]^n + u[j-1]^n)
        rhs[i] += coef_explicit * (u[j + 1] - 2.0 * u[j] + u[j - 1])
    
    # Boundary contributions (from moving known BC terms to RHS)
    # For j=1 (i=0): the term -lam*mu * u[j-1]^(n+1) = -lam*mu * u_left
    # moves to RHS as +lam*mu * u_left
    rhs[0] += lam * mu * u_left
    
    # For j=J-1 (i=n_interior-1): the term -lam*mu * u[j+1]^(n+1) = -lam*mu * u_right
    # moves to RHS as +lam*mu * u_right
    rhs[n_interior - 1] += lam * mu * u_right
    
    return rhs
