#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
using namespace std;

struct node {
   double x, y, z;
   int number;
   double f;
};
struct EL {
   int node_n[8]{};
   int number;
   double hx = 0, hy = 0, hz = 0;
   double lambda, gamma;
};


extern int nodes_c, el_c, face_c, testNumber, maxiter;
extern double lam, gam, alpha, beta_beta, norma_pr, eps;
extern vector <node> nodes;
extern vector <EL> el;
extern vector <int> faces, ig, jg;
extern vector<double> u, dif, di, gg, r, b, q, x, y, z, val, Az, Ar, Mr, b_loc;
extern vector<vector<double>> A_loc, M_loc, D;

#endif // STRUCTURES_H