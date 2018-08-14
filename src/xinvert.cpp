#include "xinvert.h"
#include <math.h>
#include "xtensor/xmath.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xindex_view.hpp"

#define XEPSILON (1e-15)

template <class E, class F>
x_mat mprod(const E& A, const F& B)
{
  return xt::sum(xt::view(A, xt::all(), xt::newaxis()) * xt::transpose(B), {2});
}

template <class E>
reals mtrace(const E& A)
{
  return xt::sum(xt::view(xt::flatten(A), xt::range(xt::placeholders::_, xt::placeholders::_, A.shape()[0] + 1)))();
}

//Written like this in hopes of eventual vectorisation.
static x_lst Bell5(reals A1, reals A2, reals A3, reals A4, reals A5)
{
  reals C1, C2, C3, C4, C5;
  C1 = A1;
  C2 = C1 * A1 + A2;
  C3 = C2 * A1 + 2 * C1 * A2 + A3;
  C4 = C3 * A1 + 3 * C2 * A2 + 3 * C1 * A3 + A4;
  C5 = C4 * A1 + 4 * C3 * A2 + 6 * C2 * A3 + 4 * C1 * A4 + A5;
  x_lst ret = {C4 / 24, C3 / -6, C2 / 2, C1 / -1, 1, C5 / 120};
  return ret;
}

x_mat xinvert5(const x_mat& A)
{
  x_mat pw[6];
  pw[0] = xt::eye(5, 0);
  reals b_arg[5];
  for(int i = 1; i < 6; i++){
    pw[i] = mprod(pw[i-1], A);
    b_arg[i-1] = mtrace(pw[i]);
  }
  b_arg[1] *= -1;
  b_arg[2] *= 2;
  b_arg[3] *= -6;
  b_arg[4] *= 24;
  auto coef = Bell5(b_arg[0], b_arg[1], b_arg[2], b_arg[3], b_arg[4]);
  reals weight = (abs(coef[5]) < XEPSILON) ? 0 : 1 / coef[5];
  // I am NOT happay about this!!!!
  // This should be done by restricting views and adding two new axes to the scalar list
  // I wadded through many many compile and runtime errors...
  return weight * (coef[0] * pw[0] + coef[1] * pw[1] + coef[2] * pw[2] + coef[3] * pw[3] + coef[4] * pw[4]);
}

void clearepsilons(x_mat& A)
{
    xt::filtration(A, xt::fabs(A) < XEPSILON) = 0;
}

// Yeah, there's no fold support and expression building via accumilate is grr
/*
  x_mat Id = xt::eye(5, 0);
  xt::xarray<reals> fold_args = xt::stack(xt::xtuple(Id, A, A, A, A));
  auto powers = xt::accumulate([&A](E X, E Y) {return mprod(X, A);}, fold_args, 1);
*/
