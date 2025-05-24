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