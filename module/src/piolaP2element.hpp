#ifndef FILE_PIOLAP2ELEMENT_HPP
#define FILE_PIOLAP2ELEMENT_HPP

namespace ngfem
{

  /*
    Finite element base class. The elements have to implement shape functions and their derivatives.
  */

  class MyBaseVecElement : public FiniteElement
  {
  public:
    MyBaseVecElement (int ndof, int order) : FiniteElement (ndof, order) {}

    virtual void CalcShape (const IntegrationPoint &ip, BareSliceMatrix<> shape) const = 0;
    virtual void CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const = 0;
  };

  class MyQuadraticVecTrig : public MyBaseVecElement
  {
  public:
    MyQuadraticVecTrig () : MyBaseVecElement (12, 2) {}
    ELEMENT_TYPE ElementType () const override { return ET_TRIG; }

    void CalcShape (const IntegrationPoint &ip, BareSliceMatrix<> shape) const override;
    void CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const override;
  };

  class MyQuadraticVecSegm : public MyBaseVecElement
  {
  public:
    MyQuadraticVecSegm () : MyBaseVecElement (6, 2) {}
    ELEMENT_TYPE ElementType () const override { return ET_SEGM; }

    void CalcShape (const IntegrationPoint &ip, BareSliceMatrix<> shape) const override;
    void CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const override;
  };

}

#endif // FILE_PIOLAP2ELEMENT_HPP
