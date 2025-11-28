#include <comp.hpp>
#include <python_comp.hpp>
#include "piolaP2fespace.hpp"

PYBIND11_MODULE (piola_lagrange_space, m)
{
  ngcomp::ExportFESpace<ngcomp::P2VecPiola> (m, "P2VecPiola");
}
