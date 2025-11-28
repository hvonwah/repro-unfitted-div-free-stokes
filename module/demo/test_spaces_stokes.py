from netgen.occ import *
from ngsolve import *
from ngsolve.krylovspace import RichardsonSolver
import piola_lagrange_space as addon
from math import log2
SetNumThreads(4)


def solve_stokes(mesh, force, uex, duex, pex, piola=True, inverse="pardiso", lagr=False):
    if piola:
        V = addon.P2VecPiola(mesh, dirichlet=".*")
    else:
        V = VectorH1(mesh, order=2, dirichlet=".*")
    Q = L2(mesh, order=1)
    if lagr is True:
        N = NumberSpace(mesh)
        X = V * Q * N
        (u, p, lam), (v, q, mu) = X.TnT()
    else:
        X = V * Q
        (u, p), (v, q) = X.TnT()
    gfu = GridFunction(X)

    stokes = (grad(u) | grad(v)) - div(v) * p  - div(u) * q
    if lagr:
        stokes += mu * p + lam * q
    else:
        stokes += 1e-10 * p * q

    a = BilinearForm(X)
    a += stokes * dx
    a.Assemble()

    f = LinearForm(X)
    f += (force * v) * dx
    f.Assemble()

    pre = a.mat.Inverse(X.FreeDofs(), inverse=inverse)
    inv = RichardsonSolver(mat=a.mat, pre=pre, tol=1e-12, maxiter=4,
                           printrates=True)
    gfu.vec.data = inv * f.vec

    vel, pre = gfu.components[0], gfu.components[1]
    err_l2u = sqrt(Integrate((vel - uex) | (vel - uex), mesh))
    err_h1u = sqrt(Integrate((grad(vel) - duex) | (grad(vel) - duex), mesh))
    err_div = sqrt(Integrate(div(vel)**2, mesh))
    err_l2p = sqrt(Integrate((pre - pex)**2, mesh))
    ndof = sum(X.FreeDofs())

    del X, V, Q, gfu, inv, pre, a, f

    return ndof, err_l2u, err_h1u, err_div, err_l2p


geo = OCCGeometry(Circle(gp_Pnt2d(0, 0), 1).Face(), dim=2)

rr = x**2 + y**2
uex = CF((y * cos(pi * rr / 2), -x * cos(pi * rr / 2)))
pex = x**3 * y**2
duex = CF((uex[0].Diff(x), uex[0].Diff(y), uex[1].Diff(x), uex[1].Diff(y)), dims=(2, 2))
force = CF((-duex[0, 0].Diff(x) + -duex[0, 1].Diff(y) + pex.Diff(x),
            -duex[1, 0].Diff(x) + -duex[1, 1].Diff(y) + pex.Diff(y)))

errs_h1u_ngs, errs_l2u_ngs, errs_l2p_ngs = [], [], []
errs_h1u_piola, errs_l2u_piola, errs_l2p_piola = [], [], []
errs_div_ngs, errs_div_piola = [], []
ndof = []

red_refine = True
h0 = 0.4
inverse = 'pardiso'
lagr = False

for i in range(5):
    if red_refine is False:
        ngmesh = geo.GenerateMesh(maxh=h0 / 2**i)
    else:
        ngmesh = geo.GenerateMesh(maxh=h0)
        for j in range(i):
            ngmesh.Refine()
    ngmesh.SplitAlfeld()
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    with TaskManager():
        print('solve isoparametric')
        out = solve_stokes(mesh, force, uex, duex, pex, False, inverse, lagr)
    errs_l2u_ngs.append(out[1])
    errs_h1u_ngs.append(out[2])
    errs_l2p_ngs.append(out[4])
    errs_div_ngs.append(out[3])

    with TaskManager():
        print('solve piola')
        out = solve_stokes(mesh, force, uex, duex, pex, True, inverse, lagr)
    errs_l2u_piola.append(out[1])
    errs_h1u_piola.append(out[2])
    errs_l2p_piola.append(out[4])
    errs_div_piola.append(out[3])
    ndof.append(out[0])

print('   ndof l2u           h1u           l2p           l2div     '
      + '| l2u           h1u           l2p           l2div')
for i in range(len(ndof)):
    if i == 0:
        l2un, h1un, l2pn, l2up, h1up, l2pp = [float('nan')] * 6
    else:
        l2un = log2(errs_l2u_ngs[i - 1] / errs_l2u_ngs[i])
        h1un = log2(errs_h1u_ngs[i - 1] / errs_h1u_ngs[i])
        l2pn = log2(errs_l2p_ngs[i - 1] / errs_l2p_ngs[i])

        l2up = log2(errs_l2u_piola[i - 1] / errs_l2u_piola[i])
        h1up = log2(errs_h1u_piola[i - 1] / errs_h1u_piola[i])
        l2pp = log2(errs_l2p_piola[i - 1] / errs_l2p_piola[i])

    print(f'{ndof[i]:7d} {errs_l2u_ngs[i]:.3e} {l2un:.1f} {errs_h1u_ngs[i]:.3e} '
          + f'{h1un:.1f} {errs_l2p_ngs[i]:.3e} {l2pn:.1f} {errs_div_ngs[i]:.3e} '
          + f'| {errs_l2u_piola[i]:.3e} {l2up:.1f} {errs_h1u_piola[i]:.3e} '
          + f'{h1up:.1f} {errs_l2p_piola[i]:.3e} {l2pp:.1f} {errs_div_piola[i]:.3e}')

    # del mesh, ngmesh
    # print(i, 'done')
    # input('pause')
