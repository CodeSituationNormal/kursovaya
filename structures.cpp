#include "common_includes.h"

int nodes_c, el_c, face_c, testNumber, maxiter = 10000;
double lam = 1, gam = 1, alpha, beta_beta, norma_pr, eps = 1e-16;
vector <node> nodes;
vector <EL> el;
vector <int> faces, ig, jg;
vector<double> u, dif, di, gg, r, b, q, x, y, z, val, Az, Ar, Mr, M, b_loc;
vector<vector<double>> A_loc, M_loc, D, alphaM;