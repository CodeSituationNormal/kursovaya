#include "common_includes.h"

double u_a(int i) {
   u[i] = nodes[i].x; // modify manually if needed
   return u[i];
}

double u_c(double x, double y, double z) {
   return x;
}

// vector<double> genVal() {
//    for (int i = 0; i < face_c; i++) {

//    }
// }

double f_auto(double x, double y, double z) {
   return x; // modify manually if needed
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

void print_u() {
   ofstream uFile("../u.txt");
   cout << endl << "u ";
   for (int i = 0; i < nodes_c; i++) {
      cout << u[i] << " ";
      uFile << scientific << setprecision(10) << u[i] << endl;
   }
   cout << endl;
   uFile.close();
}