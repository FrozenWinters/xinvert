#ifndef XINVERT_H
#define XINVERT_H

#include "xtensor/xtensor.hpp"

typedef long double reals;
typedef xt::xtensor<reals, 1> x_lst;
typedef xt::xtensor<reals, 2> x_mat;

x_mat xinvert5(const x_mat& A);

template <class E, class F>
x_mat mprod(const E& A, const F& B);

template <class E>
reals mtrace(const E& A);

void clearepsilons(x_mat& A);

#endif
