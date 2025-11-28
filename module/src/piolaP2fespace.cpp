/*

FESpace for quadratic vector valued elements triangular elements.
Elements are mapped using the Piola transformation and corrected
to be continuous at vertices and edge midpoints.

*/

#include <comp.hpp> // provides FESpace, ...

#include "piolaP2element.hpp"
#include "piolaP2fespace.hpp"
#include "piolaP2diffop.hpp"

namespace ngcomp
{

  P2VecPiola ::P2VecPiola (shared_ptr<MeshAccess> ama, const Flags &flags) : FESpace (ama, flags)
  {
    evaluator[VOL] = make_shared<T_DifferentialOperator<MyVecDiffOpId>> ();
    evaluator[BND] = make_shared<T_DifferentialOperator<MyVecDiffOpId>> ();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<MyVecDiffOpGradient>> ();
  }

  DocInfo P2VecPiola ::GetDocu ()
  {
    auto docu = FESpace::GetDocu ();
    docu.short_docu = "A Piola transformed P2 with cortrection ";
    docu.long_docu =
        R"raw_string(A vector valued second order finite element
space transformed using the Piola tranform and linear combination of
basis function to obtain conformity in Lagrange nodes.

Needed for lowest order Scott-Vogelius on curved meshes.
)raw_string";

    return docu;
  }

  void P2VecPiola ::Update ()
  {
    ne = ma->GetNE ();
    nv = ma->GetNV ();
    ned = ma->GetNEdges ();
    comp_dim = nv + ned;
    SetNDof (2 * comp_dim);

    fine_edge.SetSize (ned);
    fine_edge = 0;
    for (int i = 0; i < ne; i++)
    {
      ElementId ei (VOL, i);
      auto eledges = ma->GetElEdges (ei);
      for (int j = 0; j < eledges.Size (); j++)
      {
        fine_edge[eledges[j]] = 1;
      }
    }
    UpdateCouplingDofArray ();
  }

  void P2VecPiola ::UpdateCouplingDofArray ()
  {
    ctofdof.SetSize (ndof);
    ctofdof = WIREBASKET_DOF;
    for (int edge = 0; edge < ned; edge++)
    {
      ctofdof[nv + edge] = fine_edge[edge] ? INTERFACE_DOF : UNUSED_DOF;
      ctofdof[comp_dim + nv + edge] = fine_edge[edge] ? INTERFACE_DOF : UNUSED_DOF;
    }
  }

  void P2VecPiola ::GetDofNrs (ElementId ei, Array<DofId> &dnums) const
  {
    dnums.SetSize (0);

    for (size_t i : Range (2))
    {
      for (auto v : ma->GetElement (ei).Vertices ())
        dnums.Append (v + i * comp_dim);
      for (auto e : ma->GetElement (ei).Edges ())
        dnums.Append (nv + e + i * comp_dim);
    }
  }

  FiniteElement &P2VecPiola ::GetFE (ElementId ei, Allocator &alloc) const
  {
    switch (ma->GetElement (ei).GetType ())
    {
    case ET_TRIG:
      return *new (alloc) MyQuadraticVecTrig;
    case ET_SEGM:
      return *new (alloc) MyQuadraticVecSegm;
    default:
      throw Exception ("P2VecPiola: Element of type " + ToString (ma->GetElement (ei).GetType ())
                       + " not available\n");
    }
  }

}
