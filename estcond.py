"""
    Lehrenfeld, Christoph; Voulis, Igor; van Beeck, Tim, 2023, "estcond.py",
    Replication Data for: Analysis of divergence-preserving unfitted finite
    element methods for the mixed Poisson problem,
    https://doi.org/10.25625/LYR2CG/FMMIFA, GRO.data, V2

    Licence: CC0 1.0 Universal
"""


from ngsolve import Norm


def EstimateConditionNumber(mat, invmat, tol=5e-3, output=True):
    """
    Crude estimation of condition number based on a power iteration
    and an inverse power iteration (use only for "well-behaved" (s.p.d)
    matrices).
    Input: matrix and its inverse and a starting vector z.
    Output: condition number (estimate) as a scalar.
    """
    z = mat.CreateColVector()
    z.FV().NumPy()[:] = 1.0
    x = z.CreateVector()
    y = z.CreateVector()
    vmin, vmax = z.CreateVector(), z.CreateVector()
    y.data = z
    lmax = -1
    for i in range(10000):
        x.data = 1.0 / Norm(y) * y
        y.data = mat * x
        lmax_old = lmax
        lmax = Norm(y)
        print("\r iteration for largest EV: ", i, end="")
        if (abs(lmax - lmax_old) / abs(lmax_old) < tol):
            break
    print("\r                                    ", end="")
    vmax.data = x
    y.data = z
    lmin = -1
    for j in range(10000):
        x.data = 1.0 / Norm(y) * y
        y.data = invmat * x
        lmin_old = lmin
        lmin = 1.0 / Norm(y)
        print("\r iteration for smallest EV: ", j, end="")
        if (abs(lmin - lmin_old) / abs(lmin_old) < tol):
            break
    print("\r                                    \r", end="")
    vmin.data = x
    if output:
        print(f"needed {i} iterations for power iteration and {j} iterations for inverse power iteration.")
    return lmin, lmax, lmax / lmin, vmin, vmax
