from netgen.occ import WorkPlane, OCCGeometry
from ngsolve import *

__all__ = ['nu', 'lset', 'uex', 'duex', 'ubnd', 'pex', 'force', 'geo']

rr = x**4 + y**4
phi = sin(2 * pi * rr) / (2 * pi)
lset = rr - 1 / 4

uex = CF((-phi.Diff(y), phi.Diff(x)))
ubnd = CF((0, 0))
duex = CF((uex.Diff(x), uex.Diff(y)), dims=(2, 2)).trans
pex = sin(pi * x * y) + x**3 + y**3

ubnd = CF((0, 0))
nu = 1

force = CF((- nu * duex[0, 0].Diff(x) - nu * duex[0, 1].Diff(y) + pex.Diff(x),
            - nu * duex[1, 0].Diff(x) - nu * duex[1, 1].Diff(y) + pex.Diff(y)))

geo = OCCGeometry(WorkPlane().MoveTo(-1, -1).Rectangle(2, 2).Face(), dim=2)
