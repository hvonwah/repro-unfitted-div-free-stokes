from netgen.occ import WorkPlane, OCCGeometry
from ngsolve import *

__all__ = ['nu', 'lset', 'uex', 'duex', 'ubnd', 'pex', 'force', 'geo']

rr = x**2 + y**2
phi = sin(pi * rr / 2) / pi
lset = sqrt(rr) - 1

uex = CF((-phi.Diff(y), phi.Diff(x)))
duex = CF((uex.Diff(x), uex.Diff(y)), dims=(2, 2)).trans
pex = x**3 * y**2

ubnd = CF((0, 0))
nu = 1

force = CF((- nu * duex[0, 0].Diff(x) - nu * duex[0, 1].Diff(y) + pex.Diff(x),
            - nu * duex[1, 0].Diff(x) - nu * duex[1, 1].Diff(y) + pex.Diff(y)))

geo = OCCGeometry(WorkPlane().MoveTo(-1.5, -1.5).Rectangle(3, 3).Face(), dim=2)
