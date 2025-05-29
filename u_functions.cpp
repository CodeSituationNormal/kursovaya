#include "common_includes.h"

double u_a(int i) {
   u[i] = sin(nodes[i].x);
   return u[i];
}
void dif_u() {
   dif.resize(nodes_c);
   u.resize(nodes_c);
   ofstream difFile("../dif.txt");
   for (int i = 0; i < nodes_c; i++) {
      dif[i] = u_a(i) - q[i];
      difFile << scientific << setprecision(10) << dif[i] << endl;
   }
   difFile.close();
}