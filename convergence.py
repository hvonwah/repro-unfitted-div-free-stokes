from ngsolve import *
from math import log2
from stokes import cutfem_stokes_sv as solver
import example0_data
import example1_data
import example2_data
import argparse

parser = argparse.ArgumentParser(description='Solve convergence study for Stokes flow using unfitted Piola mapped Scott-Vogelius elements',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--example", help="Example", type=int, default=1, choices=[0, 1, 2])
parser.add_argument("--maxh", help="maxh", type=float, default=0.4)
parser.add_argument("--nref", help="number of refinement levels", type=int, default=6)
parser.add_argument("--gp", help="gamma_gp", type=float, default=0.1)
parser.add_argument("--stabil", help="facet choicee", type=str, choices=['global', 'patch'], default='global')
parser.add_argument("--sigma", help="nitsche", type=float, default=40)
parser.add_argument("--gl", help="gamma_l", type=float, default=0.1)
parser.add_argument("--hogeo", help="Piola mapping", type=int, choices=[0, 1], default=1)
parser.add_argument("--threshold", help="limit mappig", type=float, default=1)
parser.add_argument("--scalepatch", help="Doenscale GP patch", type=float, default=0.25)
parser.add_argument("--lm", help="Order of Laggrang muliplier space", type=int, default=1)
parser.add_argument("--nthreads", help="number of threads", type=int, default=4)
parser.add_argument("--vtk", help="write vtks", type=int, choices=[0, 1], default=0)
args = parser.parse_args()
options = vars(args)
print(options)

# ------------------------------------ DATA ----------------------------------
example = options['example']
h0 = options['maxh']
nref = options['nref']
gamma_gp = options["gp"]
stabil = options["stabil"]
sigma = options["sigma"]
gamma_l0 = options["gl"]
hogeo = bool(options["hogeo"])
threshold = options["threshold"]
scalepatch = options["scalepatch"]
order_lm = options['lm']
do_vtk = bool(options['vtk'])

filename = f'UnfittedStokesPiolaMappedSV_Example{example}_h{h0}'
filename += f'hogeo{int(hogeo)}{stabil}gp{gamma_gp}nitsche{sigma}'
filename += f'lagr{gamma_l0}threshold{threshold}scalepatch{scalepatch}'
filename += f'lm{order_lm}'

if example == 0:
    print('loading example 0')
    data = example0_data
elif example == 1:
    print('loading example 1')
    data = example1_data
else:  # example == 2:
    print('loading example 2')
    data = example2_data

nu = data.nu
uex, duex, ubnd, pex = data.uex, data.duex, data.ubnd, data.pex
force = data.force
lset = data.lset
geo = data.geo

# ------------------------------ CONVERGENCE STUDY ----------------------------
SetNumThreads(options['nthreads'])
err_l2v, err_h1v, err_l2d, err_l2p, err_l2pa, dofs = [[] for _ in range(6)]
with open(f'{filename}.data', 'w') as fid:
    fid.write('lvl l2u h1u l2d l2p l2p* dofs\n')

for i in range(nref):
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=h0))
        for _ in range(i):
            mesh.ngmesh.Refine()
        out = solver(mesh=mesh, lset=lset, nu=nu, rhs=force, ubnd=ubnd,
                     uex=uex, duex=duex, pex=pex, hogeo=hogeo,
                     gamma_gp=gamma_gp, sigma=sigma, gamma_l0=gamma_l0,
                     inverse='pardiso', condest=False, threshold=threshold,
                     scale_patch=min(scalepatch * (i + 1), 1), stabil=stabil,
                     order_lm=order_lm)
    err_l2v.append(out[5])
    err_h1v.append(out[6])
    err_l2d.append(out[7])
    err_l2p.append(out[8])
    err_l2pa.append(out[9])
    dofs.append(out[10])
    with open(f'{filename}.data', 'a') as fid:
        fid.write(f'{i} {out[5]} {out[6]} {out[7]} {out[8]} {out[9]} {out[10]}\n')

    if do_vtk is True:
        _vel, _pre, _lagr, _num = out[0].components
        _mesh = _vel.space.mesh
        _mesh.SetDeformation(out[2])
        vtk = VTKOutput(ma=_mesh, coefs=[1], names=['1'],
                        filename=f'vtk/{filename}mesh_lvl{i}', order=2,
                        subdivision=0)
        vtk.Do()
        vtk = VTKOutput(ma=_mesh, order=2, filename=f'vtk/{filename}_lvl{i}',
                        coefs=[_vel, _pre, div(_vel), out[1]],
                        names=['vel', 'pre', 'div', 'lsetp1'],
                        subdivision=1)
        vtk.Do(drawelems=out[3])
        vtk = VTKOutput(ma=_mesh, order=2, filename=f'vtk/{filename}lagr_lvl{i}',
                        coefs=[_lagr], names=['lagr'], subdivision=1)
        vtk.Do(drawelems=out[4])
    del out


# ------------------------------ CONVERGNCE TABLE -----------------------------
print('    dof l2v           h1v           l2p           l2pa          l2div')
for i in range(len(err_l2v)):
    if i == 0:
        l2v, h1v, l2p, l2pa = [float('nan')] * 4
    else:
        l2v = log2(err_l2v[i - 1] / err_l2v[i])
        h1v = log2(err_h1v[i - 1] / err_h1v[i])
        l2p = log2(err_l2p[i - 1] / err_l2p[i])
        l2pa = log2(err_l2pa[i - 1] / err_l2pa[i])

    print(f'{dofs[i]:7d} {err_l2v[i]:.3e} {l2v:.1f} {err_h1v[i]:.3e} '
          + f'{h1v:.1f} {err_l2p[i]:.3e} {l2p:.1f} {err_l2pa[i]:.3e} '
          + f'{l2pa:.1f} {err_l2d[i]:.3e}')
