#ifndef FILE_PIOLAP2FESPACE_HPP
#define FILE_PIOLAP2FESPACE_HPP

namespace ngcomp
{
  class P2VecPiola : public FESpace
  {
    int ne, nv, ned, comp_dim;
    Array<bool> fine_edge;

  public:
    // Constructor
    P2VecPiola (shared_ptr<MeshAccess> ama, const Flags &flags);
    // Name for FE space
    string GetClassName () const override { return "P2VecPiola"; }
    // Class documentation
    static DocInfo GetDocu ();
    // Organize FE space
    void Update () override;
    // Dof numbers from element id
    void GetDofNrs (ElementId ei, Array<DofId> &dnums) const override;
    // Generate finite element for given element id
    FiniteElement &GetFE (ElementId ei, Allocator &alloc) const override;
    // Specify wirebasket/interface/unused dofs (needed if mesh refined)
    virtual void UpdateCouplingDofArray () override;
  };

}

#endif
