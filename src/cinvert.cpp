/*
  First draft of the invert done without libraries.
*/

#include <string.h>
#include <iostream>
#include <math.h>

using namespace std;

#define at(i, j) (5 * i + j)

static void print5(long double *A){
  cout << "{ ";
  for(int i = 0; i < 5; i++){
    cout << "{ ";
    for(int j = 0; j < 5; j++){
      cout << (fabs(A[at(j,i)]) < 0.0001 ? 0 : A[at(j,i)]);
      if(j != 4)
        cout << ", ";
    }
    cout << " }";
    if(i != 4)
      cout << "," << endl << "  ";
  }
  cout << " }" << endl;;
}

static void list5(long double *x){
  for(int i = 0; i < 5; i++){
    cout << i << ": " << x[i] << endl;
  }
}

/*
  Bell_args: an input array of five long doubles
  char_poly_coef: an output array of five long doubles to be filled with the values
  of Bell polynomials, TO A FACTORIAL PRE-FACTOR, B_1 to B_5 given the bell_args
  With the correct bell_args, char_poly_coef will give the coefficients of
  the charistaristic polynomial of the matrix we're working with
*/

static void Bell5(long double *Bell_args, long double *char_poly_coef){
  #define ba Bell_args
  #define cp char_poly_coef

  cout << "Args:" << endl;
  list5(ba);

  cp[0] = ba[0];
  cp[1] = cp[0] * ba[0] + ba[1];
  cp[2] = cp[1] * ba[0] + 2 * cp[0] * ba[1] + ba[2];
  cp[3] = cp[2] * ba[0] + 3 * cp[1] * ba[1] + 3 * cp[0] * ba[2] + ba[3];
  cp[4] = cp[3] * ba[0] + 4 * cp[2] * ba[1] + 6 * cp[1] * ba[2]
    + 4 * cp[0] * ba[3] + ba[4];

  cout << "Initial vals:" << endl;
  list5(cp);

  cp[0] /= -1;
  cp[1] /= 2;
  cp[2] /= -6;
  cp[3] /= 24;
  cp[4] /= 120;

  cout << "Coefficients of char poly:" << endl;
  list5(cp);
  cout << endl;
}

static void mprod5(long double *A, long double *B, long double *prod){
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++){
      prod[at(i,j)] = 0;
      for(int k=0; k < 5; k++)
        prod[at(i,j)] += A[at(i,k)] * B[at(k,j)];
    }
}

static void msum5(long double *A, long double *B, long double *sum){
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++)
      sum[at(i,j)] = A[at(i,j)] + B[at(i,j)];
}

static void mscale5(long double lambda, long double *A, long double *res){
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++)
      res[at(i,j)] = lambda * A[at(i,j)];
}

static void mdivscale5(long double lambda, long double *A, long double *res){
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++)
      res[at(i,j)] = A[at(i,j)] / lambda;
}

static long double trace5(long double *A){
  long double res = 0;
  for(int i=0; i < 5; i++)
    res += A[at(i,i)];
  return res;
}

int invert5(long double *A, long double *res){
  long double powers[125];
  long double x[5];
  long double s[5];

  long double Id[25];
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++)
      Id[at(i,j)] = (i == j);


  memcpy(powers, A, 25 * sizeof(long double));
  for(int i = 0; i < 5; i++){
    if(i)
      mprod5(powers + 25 * (i - 1), A, powers + 25 * i);
    x[i] = trace5(powers + 25 * i);
    print5(powers + 25 * i);
  }

  x[1] *= -1;
  x[2] *= 2;
  x[3] *= -6;
  x[4] *= 24;

  Bell5(x, s);

  if(fabs(s[4]) < 0.0001)
    return 0;

  mscale5(s[3], Id, Id);
  memcpy(res, Id, 25 * sizeof(long double));

  print5(res);

  for(int i = 0; i < 3; i++){
    //cout << "+" << endl;
    mscale5(s[2 - i], powers + 25 * i, powers + 25 * i);
    //print5(powers + 25 * i);
    msum5(res, powers + 25 * i, res);
    //cout << "=" << endl;
    //print5(res);
  }

  //cout << "+" << endl;

  print5(powers + 75);

  //cout << "=" << endl;

  msum5(res, powers + 75, res);

  print5(res);

  cout << s[4] << endl;

  mdivscale5(s[4], res, res);
  //cout << "\\frac{1}{ " << s[4] << " }" << endl;

  return 1;
}

/*int main(){
  long double Mat[25];
  for(int i=0; i < 5; i++)
    for(int j=0; j < 5; j++)
      Mat[at(i,j)] = (i == j) + (i == (j-1));

  //print5(Mat);
  if(invert5(Mat, Mat)){
    print5(Mat);
  } else{
    cout << "Singular!" << endl;
  }
}*/
