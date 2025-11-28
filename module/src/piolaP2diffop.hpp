#ifndef FILE_PIOLAP2DIFFOP_HPP
#define FILE_PIOLAP2DIFFOP_HPP

#include <fem.hpp>
#include "piolaP2element.hpp"

namespace ngfem
{

// Omplementation of the identity operator
  class MyVecDiffOpId : public DiffOp<MyVecDiffOpId>
  {
  public:
    static constexpr int DIM = 1;         // dimension of the input
    static constexpr int DIM_SPACE = 2;   // dimension of the space
    static constexpr int DIM_ELEMENT = 2; // spatial dimension of the element
    static constexpr int DIM_DMAT = 2;    // dimension of the output
    static constexpr int DIFFORDER = 0;   // order of differentiation

    static bool SupportsVB (VorB checkvb) { return true; } // can do VOL and BND terms

    static const MyBaseVecElement &Cast (const FiniteElement &fel)
    {
      return static_cast<const MyBaseVecElement &> (fel);
    }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip, MAT &mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      FlatMatrix<double> shape (fel.GetNDof (), DIM_ELEMENT, lh);
      FlatMatrix<double> piolashape (fel.GetNDof (), DIM_ELEMENT, lh);
      FlatMatrix<double> beta (6, DIM_ELEMENT * DIM_ELEMENT, lh);

      // Get standard shape functions
      shape.Rows (fel.GetNDof ()).Cols (DIM_ELEMENT) = 0.0;
      Cast (fel).CalcShape (mip.IP (), shape);

      // Piola transform shape functions
      piolashape = shape * Trans (mip.GetJacobian ());
      piolashape *= (1.0 / mip.GetJacobiDet ());

      // Linear combination of basis functions to ensure C0 continuity at Lagrange nodes
      const ElementTransformation &eltrans = mip.GetTransformation ();
      SetBeta (eltrans, beta, lh);
      shape.Rows (fel.GetNDof ()).Cols (DIM_ELEMENT) = 0.0;

      for (int j = 0; j < 6; j++)
      {
        shape.Row (j) = beta (j, 0) * piolashape.Row (j);
        shape.Row (j) += beta (j, 2) * piolashape.Row (j + 6);
        shape.Row (j + 6) = beta (j, 1) * piolashape.Row (j);
        shape.Row (j + 6) += beta (j, 3) * piolashape.Row (j + 6);
      }

      mat = Trans (shape);
    }

    // Helper function
    template <typename MAT>
    static void SetBeta (const ElementTransformation &eltrans, MAT &beta, LocalHeap &lh)
    {
      IntegrationPoint ip;
      FlatMatrix<double> F (DIM_ELEMENT, DIM_ELEMENT, lh);
      FlatMatrix<double> Finv (DIM_ELEMENT, DIM_ELEMENT, lh);
      double J;

      for (int i = 0; i < 6; i++)
      {
        SetRefNodeIP (ip, i);
        eltrans.CalcJacobian (ip, F);
        CalcInverse (F, Finv);
        J = Det (F);
        beta.Row (i) = Finv.AsVector ();
        beta.Row (i) *= J;
      }
    }

    // Helper function
    static void SetRefNodeIP (IntegrationPoint &ip, int &i)
    {
      switch (i)
      {
      case 0:
        ip = IntegrationPoint (1.0, 0.0);
        return;
      case 1:
        ip = IntegrationPoint (0.0, 1.0);
        return;
      case 2:
        ip = IntegrationPoint (0.0, 0.0);
        return;
      case 3:
        ip = IntegrationPoint (0.5, 0.0);
        return;
      case 4:
        ip = IntegrationPoint (0.0, 0.5);
        return;
      case 5:
        ip = IntegrationPoint (0.5, 0.5);
        return;
      default:
        throw Exception ("Unknown index");
      }
    }
  };

  // Gradient DiffOp: implements the chain rule for mapping gradients
  // from the reference element to the pysical element
  class MyVecDiffOpGradient : public DiffOp<MyVecDiffOpGradient>
  {
  public:
    static constexpr int DIM = 1;          // dimension of the input
    static constexpr int DIM_SPACE = 2;    // dimension of the space
    static constexpr int DIM_ELEMENT = 2;  // dimension of the element
    static constexpr int DIM_DMAT = 2 * 2; // dimension of the output
    static constexpr int DIFFORDER = 1;    // order of differentiation

    // so that you can call grad(u)
    static string Name () { return "grad"; }
    static IVec<2> GetDimensions () { return { DIM_SPACE, DIM_SPACE }; }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement &fel, const MIP &mip, MAT &mat, LocalHeap &lh)
    {
      HeapReset hr (lh);
      FlatMatrix<double> shape (fel.GetNDof (), DIM_ELEMENT, lh);
      FlatMatrix<double> piolashape (fel.GetNDof (), DIM_ELEMENT, lh);
      FlatMatrix<double> dshape (fel.GetNDof (), DIM_DMAT, lh);
      FlatMatrix<double> tmpdshape (fel.GetNDof (), DIM_DMAT, lh);
      FlatMatrix<double> pioladshape (fel.GetNDof (), DIM_DMAT, lh);
      FlatVector<double> alpha (DIM_SPACE, lh);
      FlatMatrix<double> beta (6, DIM_ELEMENT * DIM_ELEMENT, lh);

      // Get standard shape functions and piola map them
      shape.Rows (fel.GetNDof ()).Cols (DIM_ELEMENT) = 0.0;
      MyVecDiffOpId ::Cast (fel).CalcShape (mip.IP (), shape);
      piolashape = shape * Trans (mip.GetJacobian ());
      piolashape *= (1.0 / mip.GetJacobiDet ());

      // Get standard shape function derivatives and apply mapping
      pioladshape.Rows (fel.GetNDof ()).Cols (DIM_DMAT) = 0.0;
      MyVecDiffOpId ::Cast (fel).CalcDShape (mip.IP (), pioladshape);
      tmpdshape.Cols (0, 2) = pioladshape.Cols (0, 2) * mip.GetJacobianInverse ();
      tmpdshape.Cols (2, 4) = pioladshape.Cols (2, 4) * mip.GetJacobianInverse ();

      // Here, Piola map derivatives of each component rather than compontent derivatives
      dshape.Cols (0, 4) = tmpdshape.Cols (0, 4);
      dshape.Col (1) = tmpdshape.Col (2);
      dshape.Col (2) = tmpdshape.Col (1);

      pioladshape.Cols (0, 2) = dshape.Cols (0, 2) * Trans (mip.GetJacobian ());
      pioladshape.Cols (2, 4) = dshape.Cols (2, 4) * Trans (mip.GetJacobian ());
      pioladshape *= (1.0 / mip.GetJacobiDet ());

      // Swap back to component derivatives rather than derivatives of each component
      tmpdshape.Col (0) = pioladshape.Col (1);
      pioladshape.Col (1) = pioladshape.Col (2);
      pioladshape.Col (2) = tmpdshape.Col (0);

      // Isoparametric correction from chain rule
      const ElementTransformation &eltrans = mip.GetTransformation ();
      if (eltrans.IsCurvedElement ())
      {
        Vec<DIM_SPACE, Mat<DIM_SPACE, DIM_SPACE>> hesse = mip.CalcHesse ();
        Mat<DIM_SPACE, DIM_SPACE> tmp0 = hesse (0) * mip.GetJacobianInverse ();
        Mat<DIM_SPACE, DIM_SPACE> tmp1 = hesse (1) * mip.GetJacobianInverse ();
        alpha = 0;
        for (size_t j : Range (DIM_SPACE))
        {
          for (size_t r : Range (DIM_SPACE))
          {
            alpha (j) += mip.GetJacobianInverse () (r, 0) * tmp0 (r, j);
            alpha (j) += mip.GetJacobianInverse () (r, 1) * tmp1 (r, j);
          }
        }
        // fyi: alpha computed above identical to when computed with loops
        tmpdshape.Cols (0, 2) = shape * tmp0;
        tmpdshape.Cols (2, 4) = shape * tmp1;
        tmpdshape *= 1 / mip.GetJacobiDet ();
        pioladshape.Cols (0, 4) += tmpdshape.Cols (0, 4);

        pioladshape.Col (0) -= alpha(0) * piolashape.Col (0);
        pioladshape.Col (1) -= alpha(1) * piolashape.Col (0);
        pioladshape.Col (2) -= alpha(0) * piolashape.Col (1);
        pioladshape.Col (3) -= alpha(1) * piolashape.Col (1);
      }

      // Linear combination of basis functions to ensure C0 continuity at Lagrange nodes
      MyVecDiffOpId ::SetBeta (eltrans, beta, lh);
      dshape.Rows (fel.GetNDof ()).Cols (DIM_DMAT) = 0.0;

      size_t i0, i1;
      for (int j = 0; j < 6; j++)
      {
        for (size_t i = 0; i < DIM_SPACE; i++)
        {
          i0 = 2 * i;
          i1 = 2 * (i + 1); 
          dshape.Cols (i0, i1).Row (j) = beta (j, 0) * pioladshape.Cols (i0, i1).Row (j);
          dshape.Cols (i0, i1).Row (j) += beta (j, 2) * pioladshape.Cols (i0, i1).Row (j + 6);
          dshape.Cols (i0, i1).Row (j + 6) = beta (j, 1) * pioladshape.Cols (i0, i1).Row (j);
          dshape.Cols (i0, i1).Row (j + 6) += beta (j, 3) * pioladshape.Cols (i0, i1).Row (j + 6);
        }
      }
      mat = Trans (dshape);
    }
  };
} // namespace ngfem

#endif // FILE_PIOLAP2DIFFOP_HPP
