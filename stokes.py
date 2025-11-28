
from ngsolve import *
from ngsolve import dx as dxngs
from ngsolve.krylovspace import RichardsonSolver
from xfem import *
from xfem.lsetcurv import *
from piola_lagrange_space import P2VecPiola
from netgen.meshing import Pnt, FaceDescriptor, Element1D, Element2D, MeshPoint
from netgen.meshing import Mesh as ngmesh
from estcond import EstimateConditionNumber

__all__ = ['cutfem_stokes_sv']


def alfeld_split_2d(mesh):
    ngsmesh = mesh
    mesh = mesh.ngmesh
    newmesh = ngmesh()
    newmesh.SetGeometry(mesh.GetGeometry())
    newmesh.dim = mesh.dim

    for p in mesh.Points():
        newmesh.Add(p)

    def average_pnt(lv):
        ll = len(lv)
        xd_coord = [0, 0]
        for d in [0, 1]:
            xd_coord[d] = sum([mesh.Points()[v].p[d] for v in lv]) / ll
        return Pnt(xd_coord[0], xd_coord[1], 0.0)

    newmesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))

    for edge in mesh.Elements1D():
        v1, v2 = edge.vertices
        newmesh.Add(Element1D(vertices=[v1, v2], index=edge.index))

    for i, bnd in enumerate(ngsmesh.GetBoundaries()):
        newmesh.SetBCName(i, bnd)

    for i, el in enumerate(mesh.Elements2D()):
        v1, v2, v3 = el.vertices
        cellcenter = average_pnt(el.vertices)
        v123 = newmesh.Add(MeshPoint(cellcenter))
        newmesh.Add(Element2D(index=el.index, vertices=[v1, v2, v123]))
        newmesh.Add(Element2D(index=el.index, vertices=[v3, v1, v123]))
        newmesh.Add(Element2D(index=el.index, vertices=[v2, v3, v123]))

    newmesh.Update()

    return newmesh


def cutfem_stokes_sv(mesh: ngsolve.Mesh, lset: ngsolve.CoefficientFunction, nu: float = 1.0, rhs: ngsolve.CoefficientFunction = CF((0, 0)), ubnd: ngsolve.CoefficientFunction = CF((0, 0)), uex: ngsolve.CoefficientFunction = CF((0, 0)), duex: ngsolve.CoefficientFunction = CF((0, 0, 0, 0), dims=(2, 2)), pex: ngsolve.CoefficientFunction = CF(0), hogeo: bool = True, gamma_gp: float = 0.1, sigma: float = 10, gamma_l0: float = 0.1, inverse: str = 'pardiso', condest: bool = False, threshold: float = 0.05, scale_patch: float = 1, stabil: str = 'global', order_lm: int = 1):
    # -------------------------- ALFELD SPLIT MESH ----------------------------
    mesh_m = Mesh(mesh.ngmesh.Copy())
    mesh = Mesh(alfeld_split_2d(mesh_m))

    # -------------------------- LEVELSET & CUT-INFO --------------------------
    if hogeo:
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=2,
                                              threshold=threshold,
                                              discontinuous_qn=True,
                                              levelset=lset)
    else:
        lsetmeshadap = NoDeformation(mesh, lset, False)
    lsetp1 = lsetmeshadap.lset_p1
    deform = lsetmeshadap.deform

    ci = CutInfo(mesh, lsetp1)

    lsetp1_m = GridFunction(H1(mesh_m, order=1))
    InterpolateToP1(lset, lsetp1_m)
    ci_m = CutInfo(mesh_m, lsetp1_m)

    # ---------------------------- ELEMENT MARKERS ----------------------------
    els_hasneg = ci.GetElementsOfType(HASNEG)
    els_if = ci.GetElementsOfType(IF)

    els_hasneg_m = ci_m.GetElementsOfType(HASNEG)
    els_if_m = ci_m.GetElementsOfType(IF)

    facets_gp = BitArray(mesh.nfacet)

    if stabil == 'global':
        facets_gp_m = GetFacetsWithNeighborTypes(mesh_m, a=els_hasneg_m,
                                                 b=els_if_m, use_and=True)
        els_gp_m = GetElementsWithNeighborFacets(mesh_m, facets_gp_m)

        els_gp = BitArray(mesh.ne)
        els_gp.Clear()
        for i in range(mesh_m.ne):
            if els_gp_m[i]:
                els_gp[3 * i: 3 * (i + 1)] = True
        facets_gp[:] = GetFacetsWithNeighborTypes(mesh, a=els_gp, b=els_gp,
                                                  use_and=True)
    elif stabil == 'patch':
        els_neg_m = ci_m.GetElementsOfType(NEG)
        ea_m = ElementAggregation(mesh_m, els_neg_m, els_if_m)
        n_patch = max(ea_m.element_to_patch) + 1
        els_patches = [BitArray(mesh.ne) for _ in range(n_patch)]
        for i in range(n_patch):
            els_patches[i].Clear()
        for j in range(mesh_m.ne):
            i = ea_m.element_to_patch[j]
            if i != -1:
                els_patches[i][3 * j: 3 * (j + 1)] = True
        facets_gp.Clear()
        for els_patch in els_patches:
            facets_gp |= GetFacetsWithNeighborTypes(mesh, a=els_patch,
                                                    b=els_patch, use_and=True)
    else:
        raise ValueError('Unknown gp facet choice')

    els_outer = BitArray(mesh.ne)
    els_outer.Clear()
    for i in range(mesh_m.ne):
        if els_hasneg_m[i]:
            els_outer[3 * i: 3 * (i + 1)] = True

    facets_if = GetFacetsWithNeighborTypes(mesh, a=els_if, b=els_if,
                                           use_and=True)
    facets_active = BitArray(mesh.nfacet)
    facets_active[:] = facets_if | facets_gp

    # ------------------------- FINITE ELEMENT SPACE --------------------------
    V = P2VecPiola(mesh)
    Q = L2(mesh, order=1)
    L = H1(mesh, order=order_lm)
    N = NumberSpace(mesh)
    X = FESpace([V, Q, L, N], dgjumps=True)
    (u, p, lam, r), (v, q, mu, s) = X.TnT()
    gfu = GridFunction(X)
    freedofs = CompoundBitArray([GetDofsOfElements(V, els_outer),
                                 GetDofsOfElements(Q, els_outer),
                                 GetDofsOfElements(L, els_if),
                                 N.FreeDofs()])
    ndof = sum(freedofs)

    # ------------------------------ INTEGRATORS ------------------------------
    dx = dCut(levelset=lsetp1, domain_type=NEG, definedonelements=els_hasneg,
              deformation=deform)
    ds = dCut(levelset=lsetp1, domain_type=IF, definedonelements=els_if,
              deformation=deform)
    dw = dFacetPatch(definedonelements=facets_gp, deformation=deform,
                     downscale=scale_patch)
    dX = dxngs(definedonelements=els_outer, deformation=deform)
    dXint = dxngs(definedonelements=els_hasneg, deformation=deform)
    dXG = dxngs(definedonelements=els_if, deformation=deform)

    # --------------------------- (BI)LINEAR FORMS ----------------------------
    h = specialcf.mesh_size
    n_lset = 1.0 / Norm(grad(lsetp1)) * grad(lsetp1)

    stokes = nu * (grad(u) | grad(v))
    div_constraint = - div(v) * p - div(u) * q + p * s + q * r - 1e-11 * p * q
    nitsche = - nu * (grad(u) * n_lset) * v - nu * (grad(v) * n_lset) * u
    nitsche += nu * (sigma * 2**2 / h) * u * v
    nitsche += lam * n_lset * v + mu * n_lset * u

    ghost_penalty = nu * gamma_gp / h**2 * (u - u.Other()) * (v - v.Other())

    stabil_lagr0 = - gamma_l0 * h * (n_lset * Grad(lam)) * (n_lset * Grad(mu))

    forcing = rhs * v
    f_bnd = nu * (grad(v) * n_lset) * ubnd
    f_bnd += nu * (sigma * 2**2 / h) * ubnd * v
    f_bnd += mu * n_lset * ubnd

    # ------------------------------ INTEGRATORS ------------------------------
    a = RestrictedBilinearForm(X, element_restriction=els_outer,
                               facet_restriction=facets_active,
                               check_unused=False)
    a += stokes.Compile() * dx
    a += nitsche.Compile() * ds
    a += ghost_penalty.Compile() * dw
    a += div_constraint.Compile() * dX
    a += stabil_lagr0.Compile() * dXG
    a.Assemble()
    pre = a.mat.Inverse(freedofs=freedofs, inverse=inverse)

    f = LinearForm(X)
    f += forcing.Compile() * dx
    f += f_bnd.Compile() * ds
    f.Assemble()

    # ----------------------------- WRITE TO FILE -----------------------------
    if condest is True:
        lam1, lamn, kappa, vmin, vmax = EstimateConditionNumber(a.mat, pre,
                                                                1e-6, True)
    else:
        kappa = float('Nan')

    # ----------------------------- SOLVE PROBLEM -----------------------------
    inv = RichardsonSolver(mat=a.mat, pre=pre, maxiter=4, tol=1e-20,
                           printrates=True)
    gfu.vec.data = inv * f.vec

    # --------------------------- ERROR COMPUTATION ---------------------------
    vel, pre = gfu.components[0], gfu.components[1]
    l2v = sqrt(Integrate((vel - uex)**2 * dx, mesh))
    h1v = sqrt(Integrate(((grad(vel) - duex) | (grad(vel) - duex)) * dx, mesh))
    l2d = sqrt(Integrate(div(vel)**2 * dX, mesh))
    l2p = sqrt(Integrate((pre - pex)**2 * dXint, mesh))

    # -------------------------- POSTPROCESS PRESSURE -------------------------
    Qast = H1(mesh, order=1, dgjumps=True)
    Xast = Qast * N

    (ps, rs), (qs, ss) = Xast.TnT()
    gfuast = GridFunction(Xast)

    active_dofs = CompoundBitArray([GetDofsOfElements(Qast, els_hasneg),
                                    N.FreeDofs()])

    facets_gp[:] = GetFacetsWithNeighborTypes(mesh, a=els_hasneg,
                                              b=els_if, use_and=True)
    dw = dFacetPatch(definedonelements=facets_gp, deformation=deform)

    t_lset = CF((-n_lset[1], n_lset[0]))
    grad_vel = CF(Grad(vel), dims=(2, 2))
    curl_vel = grad_vel[0, 1] - grad_vel[1, 0]

    poisson = Grad(ps) * Grad(qs) + ps * ss + rs * qs
    ghost_penalty = gamma_gp / h**2 * (ps - ps.Other()) * (qs - qs.Other())

    forcing = rhs * Grad(qs)
    forcing_bnd = nu * curl_vel * t_lset * Grad(qs)

    a_ast = RestrictedBilinearForm(Xast, element_restriction=els_hasneg,
                                   facet_restriction=facets_gp,
                                   check_unused=False)
    a_ast += poisson.Compile() * dx
    a_ast += ghost_penalty.Compile() * dw
    a_ast.Assemble()

    f_ast = LinearForm(Xast)
    f_ast += forcing.Compile() * dx
    f_ast += forcing_bnd.Compile() * ds
    f_ast.Assemble()

    pres = a_ast.mat.Inverse(freedofs=active_dofs, inverse=inverse)
    invs = RichardsonSolver(mat=a_ast.mat, pre=pres, maxiter=4, tol=1e-19,
                            printrates=True)
    gfuast.vec.data = invs * f_ast.vec

    preast = gfuast.components[0]

    Draw(BitArrayCF(els_hasneg), mesh, 'hasneg')
    Draw(pre * IfPos(lsetp1, float('nan'), 1), mesh, 'pre')
    Draw(pex * IfPos(lsetp1, float('nan'), 1), mesh, 'preex')
    Draw(preast * IfPos(lsetp1, float('nan'), 1), mesh, 'prea')
    Draw((preast - pex) * IfPos(lsetp1, float('nan'), 1), mesh, 'prea-err')
    Draw(IfPos(lsetp1, float('nan'), 1) * vel, mesh, 'vel')
    # input('pause')
    l2pa = sqrt(Integrate((preast - pex)**2 * dx, mesh))

    return gfu, lsetp1, deform, els_outer, els_if, l2v, h1v, l2d, l2p, \
        l2pa, ndof, kappa
