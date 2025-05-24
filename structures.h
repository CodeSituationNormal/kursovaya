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

// extern int el_n = 0, nodes_n = 0, face_n = 0, maxiter = 100000;
// extern double eps = 1e-14, alpha, beta_beta, norma_pr;

// extern double* b_loc = nullptr;
// extern double** D = nullptr;
// extern double** alphaM = nullptr;
// extern double** A_loc = nullptr;
// extern double** M_loc = nullptr;

// extern vector<int> face, ig, jg;
// extern vector<double> val, gg, di , b, q, r, z, Az, Ar, Mr, M;

// extern node* nodes = nullptr; 
// extern EL* el = nullptr;