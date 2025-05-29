#include "common_includes.h"

int nx, ny, nz, nt, bc_left, bc_right, n;
double x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max, kx, ky, kz, kt, hx, hy, hz, ht;
vector<int> bc1;
vector<double> t; 

ofstream nodes_out, elements_out, f_out, faces_out;

void buildGrid() {
   nodes_out.open("../nodes_out.txt");
   elements_out.open("../elements_out.txt");
   f_out.open("../f.txt");
   faces_out.open("../faces.txt");

   vector<double> x_coords(nx);
   vector<double> y_coords(ny);
   vector<double> z_coords(nz);

   double sum;
   if (kx == 1.0) {
      double hx = (x_max - x_min) / (nx - 1);
      for (int i = 0; i < nx; ++i)
         x_coords[i] = x_min + i * hx;
   }
   else {
      sum = (1.0 - pow(kx, nx - 1)) / (1.0 - kx);
      double hx = (x_max - x_min) / sum;
      x_coords[0] = x_min;
      for (int i = 1; i < nx; ++i)
         x_coords[i] = x_coords[i - 1] + hx * pow(kx, i - 1);
   }

   if (ky == 1.0) {
      double hy = (y_max - y_min) / (ny - 1);
      for (int j = 0; j < ny; ++j)
         y_coords[j] = y_min + j * hy;
   }
   else {
      sum = (1.0 - pow(ky, ny - 1)) / (1.0 - ky);
      double hy = (y_max - y_min) / sum;
      y_coords[0] = y_min;
      for (int j = 1; j < ny; ++j)
         y_coords[j] = y_coords[j - 1] + hy * pow(ky, j - 1);
   }

   if (kz == 1.0) {
      double hz = (z_max - z_min) / (nz - 1);
      for (int k = 0; k < nz; ++k)
         z_coords[k] = z_min + k * hz;
   }
   else {
      sum = (1.0 - pow(kz, nz - 1)) / (1.0 - kz);
      double hz = (z_max - z_min) / sum;
      z_coords[0] = z_min;
      for (int k = 1; k < nz; ++k)
         z_coords[k] = z_coords[k - 1] + hz * pow(kz, k - 1);
   }

   nodes.clear(); 
   for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
         for (int i = 0; i < nx; ++i) {
            node n;
            n.x = x_coords[i];
            n.y = y_coords[j];
            n.z = z_coords[k];
            n.f = f_auto(n.x, n.y, n.z); 
            n.number = nodes.size();
            nodes.push_back(n);
         }
      }
   }

   double h_temp = 0;
   for (int j = 1; j < nt; ++j) {
      h_temp = ht * pow(kt, j - 1);
      t[j] = t[j - 1] + h_temp;
   }
   int n_xy = nx * ny;
   cout << "Elements grid:" << endl;
   
   for (int i = 0; i < nx - 1; i++) {
      for (int j = 0; j < ny - 1; j++) {
         for (int k = 0; k < nz - 1; k++) {
            EL el_cube{};
            el_cube.node_n[0] = k * n_xy + j * ny + i;
            el_cube.node_n[1] = k * n_xy + j * ny + i + 1;
            el_cube.node_n[2] = k * n_xy + (j + 1) * nx + i;
            el_cube.node_n[3] = k * n_xy + (j + 1) * nx + i + 1;
            el_cube.node_n[4] = (k + 1) * n_xy + j * ny + i;
            el_cube.node_n[5] = (k + 1) * n_xy + j * ny + i + 1;
            el_cube.node_n[6] = (k + 1) * n_xy + (j + 1) * nx + i;
            el_cube.node_n[7] = (k + 1) * n_xy + (j + 1) * nx + i + 1;
            el.push_back(el_cube);
            cout << el_cube.node_n[0] << " " << el_cube.node_n[1] << " " << el_cube.node_n[2] << " " << el_cube.node_n[3] << " " << el_cube.node_n[4] << " " << el_cube.node_n[5] << " " << el_cube.node_n[6] << " " << el_cube.node_n[7] << endl;
            elements_out << el_cube.node_n[0] << " " << el_cube.node_n[1] << " " << el_cube.node_n[2] << " " << el_cube.node_n[3] << " " << el_cube.node_n[4] << " " << el_cube.node_n[5] << " " << el_cube.node_n[6] << " " << el_cube.node_n[7] << endl;
         }
      }
   }
   cout << endl;
   cout << "Coordinates: " << endl;
   for (const auto& n : nodes) {
      cout << "Node #" << n.number << ": "
         << "x = " << n.x << ", "
         << "y = " << n.y << ", "
         << "z = " << n.z << endl;
      nodes_out << n.x << " " << n.y << " " << n.z << endl;
      f_out << scientific << setprecision(10) << n.f << endl;
   }
   for (const auto& node_n : nodes) if ((node_n.x == x_min) || (node_n.x == x_max) || (node_n.y == y_min) || (node_n.y == y_max) || (node_n.z == z_min) || (node_n.z == z_max))
      bc1.push_back(node_n.number);
   for (const auto& node_num : bc1) {
      cout << node_num << " ";
      faces_out << node_num << " " << u_c(nodes[node_num].x, nodes[node_num].y, nodes[node_num].z) << endl;
   }
   cout << endl;

   nodes_out.close();
   elements_out.close();
   f_out.close();
   faces_out.close();
}
void input() {
   ifstream inputGrid("../grid.txt");
   ifstream inputNodes("../nodes.txt");
   ifstream inputBc("../bc.txt");

   if (!inputGrid.is_open() || !inputNodes.is_open() || !inputBc.is_open()) {
      cerr << "Unable to open input files." << endl;
      exit(1);
   }

   inputGrid >> x_min >> x_max >>y_min>>y_max>>z_min>>z_max>> t_min >> t_max;
   inputNodes >> nx >>ny>>nz>> nt >> kx >>ky>>kz>> kt;
   // inputBc >> bc_left >> bc_right;
   t.resize(nt);

   buildGrid();

   inputBc.close();
   inputGrid.close();
   inputNodes.close();

   nodes.clear();
   el.clear();
   faces.clear();
}

int main() {
   input();
   return 0;
}