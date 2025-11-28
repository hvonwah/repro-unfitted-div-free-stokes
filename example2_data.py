from netgen.occ import WorkPlane, OCCGeometry
from ngsolve import *

__all__ = ['nu', 'lset', 'uex', 'duex', 'ubnd', 'pex', 'force', 'geo']

rr = x**2 + y**2
phi = sin(pi * rr / 2) / pi
lset = sqrt(x**2 + y**2) - (0.7 + 0.2 * cos(6 * atan2(y, x)))

uex = uex = CF((0, 0))
duex = CF((0, 0, 0, 0), dims=(2, 2))
pex = x**5 + y**5

ubnd = CF((0, 0))
nu = 1

force = CF((pex.Diff(x), pex.Diff(y)))

geo = OCCGeometry(WorkPlane().MoveTo(-1, -1).Rectangle(2, 2).Face(), dim=2)
