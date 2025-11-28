from netgen.occ import WorkPlane, OCCGeometry
from ngsolve import *
from stokes import cutfem_stokes_sv as solver
import argparse
SetNumThreads(4)

parser = argparse.ArgumentParser(description='Compute condition number estimates for Stokes system using unfitted Piola mapped Scott-Vogelius elements',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--maxh", help="maxh", type=float, default=0.1)
parser.add_argument("--lm", help="Order of Laggrang muliplier space", type=int, default=1)
parser.add_argument("--steps", help="number of shift steps ", type=int, default=100)
args = parser.parse_args()
options = vars(args)
print(options)

# ------------------------------------ DATA ----------------------------------
example = options['maxh']
h0 = 0.1
gamma_gp = 0.1
sigma = 40
gamma_l0 = 0.1
hogeo = True
threshold = 1
scale_patch = 0.25
order_lm = options['lm']

filename = f'UnfittedStokesPiolaMappedSV_Example{example}_h{h0}hogeo'
filename += f'{int(hogeo)}gp{gamma_gp}nitsche{sigma}lagr{gamma_l0}theshold'
filename += f'{threshold}scale_patch{scale_patch}lm{order_lm}condest.txt'


# ------------------------------ CONVERGENCE STUDY ----------------------------
_x = Parameter(0.0)
rr = (x - _x)**4 + y**4
lset = rr - 1 / 4
nu = 1

uex = CF((0, 0))
ubnd = CF((0, 0))
duex = CF((0, 0, 0, 0), dims=(2, 2)).trans
pex = CF(0)
force = CF((0, 0))

geo = OCCGeometry(WorkPlane().MoveTo(-1, -1).Rectangle(2, 2).Face(), dim=2)

mesh = Mesh(geo.GenerateMesh(maxh=h0))

_n = options['steps']
offset = 0.2
with open(filename, 'w') as fid:
    fid.write('i condest\n')

for i in range(_n + 1):
    _x.Set(-offset + i * 2 * offset / _n)

    with TaskManager():
        out = solver(mesh=mesh, lset=lset, nu=nu, rhs=force, ubnd=ubnd,
                     uex=uex, duex=duex, pex=pex, hogeo=hogeo,
                     gamma_gp=gamma_gp, sigma=sigma, gamma_l0=gamma_l0,
                     inverse='pardiso', condest=True, threshold=threshold,
                     scale_patch=scale_patch, order_lm=order_lm)
    with open(filename, 'a') as fid:
        fid.write(f'{i} {out[-1]:.5e}\n')
