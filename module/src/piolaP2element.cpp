/*
Shape functions of vector valuer quadratic elements on triangles
*/

#include <fem.hpp>
#include "piolaP2element.hpp"

namespace ngfem
{

  // Shape functions -> 6x2 matrix
  void MyQuadraticVecTrig ::CalcShape (const IntegrationPoint &ip, BareSliceMatrix<> shape) const
  {
    // barycentric ccordinates
    double lam[3] = { ip (0), ip (1), 1 - ip (0) - ip (1) };
    const EDGE *edges = ElementTopology::GetEdges (ET_TRIG);

    // Same basis functions in each component
    for (int j = 0; j < 2; j++)
    {
      // vertex basis functions
      for (int i = 0; i < 3; i++)
        shape (i + j * 6, j) = lam[i] * (2 * lam[i] - 1);

      // edge basis functions
      for (int i = 0; i < 3; i++)
        shape (3 + i + j * 6, j) = 4 * lam[edges[i][0]] * lam[edges[i][1]];
    }
  }

  // Differentiated shape functions -> 6x4 matrix
  void MyQuadraticVecTrig ::CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const
  {
    AutoDiff<2> x (ip (0), 0); // value of x, grad is 0-th unit vec (1, 0)
    AutoDiff<2> y (ip (1), 1); // value of y, grad is 1-th unit vec (0, 1)
    AutoDiff<2> lam[3] = { x, y, 1 - x - y };
    const EDGE *edges = ElementTopology::GetEdges (ET_TRIG);

    // Same basis functions in each component
    for (int j = 0; j < 2; j++)
    {
      // vertex basis functions:
      for (int i = 0; i < 3; i++)
      {
        AutoDiff<2> shape = lam[i] * (2 * lam[i] - 1);
        dshape (i + 6 * j, 0 + 2 * j) = shape.DValue (0); // x-derivative
        dshape (i + 6 * j, 1 + 2 * j) = shape.DValue (1); // y-derivative
      }
      // edge basis functions:
      for (int i = 0; i < 3; i++)
      {
        AutoDiff<2> shape = 4 * lam[edges[i][0]] * lam[edges[i][1]];
        dshape (3 + i + 6 * j, 0 + 2 * j) = shape.DValue (0); // x-derivative
        dshape (3 + i + 6 * j, 1 + 2 * j) = shape.DValue (1); // y-derivative
      }
    }
  }

  // Shape functions -> 6x2 matrix
  void MyQuadraticVecSegm ::CalcShape (const IntegrationPoint &ip, BareSliceMatrix<> shape) const
  {
    // barycentric ccodinates on line
    double lam[2] = { ip (0), 1 - ip (0) };

    // Same basis functions in each component
    for (int j = 0; j < 2; j++)
    {
      shape (0 + 3 * j, j) = lam[0] * (2 * lam[0] - 1);
      shape (1 + 3 * j, j) = lam[1] * (2 * lam[1] - 1);
      shape (2 + 3 * j, j) = 4 * lam[0] * lam[1];
    }
  }

  // Differentiated shape functions -> 6x2 matrix
  void MyQuadraticVecSegm ::CalcDShape (const IntegrationPoint &ip, BareSliceMatrix<> dshape) const
  {
    AutoDiff<1> x (ip (0));
    AutoDiff<1> lam[2] = { x, 1 - x };
    AutoDiff<1> shape;
    for (int j = 0; j < 2; j++)
    {
      shape = lam[0] * (2 * lam[0] - 1);
      dshape (0 + 3 * j, j) = shape.DValue (0);

      shape = lam[1] * (2 * lam[1] - 1);
      dshape (1 + 3 * j, j) = shape.DValue (0);

      shape = 4 * lam[0] * lam[1];
      dshape (2 + 3 * j, j) = shape.DValue (0);
    }
  }

}
