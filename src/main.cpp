#include "xinvert.h"
#include "xtensor/xio.hpp"
#include "xtensor/xrandom.hpp"
#include <iostream>
#include <time.h>

using namespace std;

int main()
{
  /*
  some tests:
  x_mat A = xt::eye(5,1); // A singular enample
  x_mat A = xt::eye(5,1) + xt::eye(5,0);
  x_mat A = xt::eye(5,1) * 15024 + xt::eye(5,0) * -234800 + xt::eye(5,-1) * 314159 + xt::eye(5, 4) * 2;
  //Numbers of different magnitudes
  */

  xt::random::seed(time(NULL));
  x_mat A = xt::random::rand<reals>({5, 5}, -100, 100);
  x_mat B = xinvert5(A);
  x_mat C = mprod(A, B);
  clearepsilons(C);
  cout << "A = " << endl << A << endl;
  cout << "A^{-1} = " << endl << B << endl;
  cout << "A*A^{-1} = " << endl << C << endl;
}
