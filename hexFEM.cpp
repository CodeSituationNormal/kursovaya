#include "common_includes.h"

void portrait() {
   map<int, set<int>> list;
   for (int i = 0; i < el_n; ++i) {
      for (int j = 0; j < 4; ++j) {
         int node_index = el[i].node_n[j];
         list[node_index].insert(node_index);
         for (int k = 0; k < 4; ++k) {
            int neighbor_index = el[i].node_n[k];
            if (node_index != neighbor_index) {
               list[node_index].insert(neighbor_index);
            }
         }
      }
   }
   
  /* for (size_t i = 0; i < list.size(); ++i) {
      cout << "Node " << i << ": ";
      for (int node : list[i]) 
         cout << node << " ";
      cout << endl;
   }*/
   int k = 0;
   ig.push_back(0);
   for (int i = 0; i < nodes_n; i++) {
      k = 0;
      for (int node : list[i]) {
         if (node < i) {
            jg.push_back(node);
            k++;
         }
      }
      ig.push_back(ig.back() + k);
   }
  /* cout << "jg " << ": ";
   for (int node : jg) {
      cout << node << " ";
   }
   cout << endl;
   cout << "ig " << ": ";
   for (int node : ig) {
      cout << node << " ";
   }
   cout << endl;*/
}

void trans_Gauss(double** A, int i, int n) {
   int line = i;
   for (int j = i + 1; j < n; j++)
      if (fabs(A[j][i]) > fabs(A[line][i]))
         line = j;
   if (line != i) 
      for (int j = 0; j < 2 * n; j++)
         swap(A[i][j], A[line][j]);
}
void GaussJordan(double** A, double** alpha, int n) {
   double** matrix = new double* [n];
   for (int i = 0; i < n; i++) {
      matrix[i] = new double[2 * n];
      for (int j = 0; j < n; j++) 
         matrix[i][j] = A[i][j];
      for (int j = 0; j < n; j++) 
         matrix[i][j + n] = (j == i) ? 1.0 : 0.0;
   }
   for (int i = 0; i < n; i++) {
      trans_Gauss(matrix, i, n);
      double m_d = matrix[i][i];
      for (int j = 0; j < 2 * n; j++)
         matrix[i][j] /= m_d;
      for (int j = 0; j < n; j++) {
         if (j != i) {
            double m_j = matrix[j][i];
            for (int k = 0; k < 2 * n; k++)
               matrix[j][k] -= matrix[i][k] * m_j;
         }
      }
   }
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) 
         alpha[i][j] = matrix[i][j + n]; 
   }
   for (int i = 0; i < n; i++) delete[] matrix[i];
   delete[] matrix;
}

void local_el(int el_n) {
   double mes = 0;
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
         alphaM[i][j] = 0;
         A_loc[i][j] = 0;
         M_loc[i][j] = 0;

      }
      b_loc[i] = 0;
   }
   for (int j = 0; j < 4; j++)
      D[0][j] = 1;
   for (int j = 0; j < 4; j++)
      D[1][j] = nodes[el[el_n].node_n[j]].x;
   for (int j = 0; j < 4; j++)
      D[2][j] = nodes[el[el_n].node_n[j]].y;
   for (int j = 0; j < 4; j++)
      D[3][j] = nodes[el[el_n].node_n[j]].z;
  /* cout << "D[" << el_n << "]" << endl;
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++)
         cout << D[i][j] << " ";
      cout << endl;
   }*/
   mes = (D[1][1] - D[1][0]) * ((D[2][2] - D[2][0]) * (D[3][3] - D[3][0]) - (D[2][3] - D[2][0]) * (D[3][2] - D[3][0]))
      + (D[2][1] - D[2][0]) * ((D[3][2] - D[3][0]) * (D[1][3] - D[1][0]) - (D[3][3] - D[3][0]) * (D[1][2] - D[1][0]))
      + (D[3][1] - D[3][0]) * ((D[1][2] - D[1][0]) * (D[2][3] - D[2][0]) - (D[1][3] - D[1][0]) * (D[2][2] - D[2][0]));

   mes = fabs(mes) / 6;
   GaussJordan(D, alphaM, 4);

  /* cout << "alpha[" << el_n << "]" << endl;
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++)
         cout << alphaM[i][j] << " ";
      cout << endl;
   }*/
   double lambda = el[el_n].lambda;
   double gamma = el[el_n].gamma;
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
         if (i == j) M_loc[i][i] = mes / 10;
         else M_loc[i][j] = mes / 20;
      }
   }
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
         double sum = 0;
         for (int k = 1; k < 4; k++)
            sum += alphaM[i][k] * alphaM[j][k];
         A_loc[i][j] = (lambda * mes*sum + gamma * M_loc[i][j]);
         b_loc[i] += nodes[el[el_n].node_n[j]].f * M_loc[i][j];
      }
   }
   /*for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++)
         cout << A_loc[i][j] << " ";
      cout << endl;
   }*/
   /*for (int i = 0; i < 4; i++) 
      for (int j = 0; j < 4; j++) 
         b_loc[i] += nodes[el[el_n].node_n[j]].f * M_loc[i][j];*/
 
   //for (int i = 0; i < 4; i++) f_loc += alphaM[i][1] * nodes[el[el_n].node_n[i]].f;

   /*double f_loc_new = f_loc * mes / 4.0;
   fill(b_loc, b_loc + 4, f_loc_new);*/
}
int find_el_pos(int i, int j) {
   int s = ig[i], e = ig[i + 1];
   int cur = 0;
   bool find = false;
   for (int p = s; p < e && !find; p++) {
      if (jg[p] == j) {
         cur = p;
         find = true;
      }
   }
   return cur;
}
void global_A() {
   int i_gg, glob_i, glob_j;
   A_loc = new double* [4];
   M_loc = new double* [4];
   D = new double* [4];
   alphaM = new double* [4];
   b_loc = new double[4]{};
   di.resize(nodes_n, 0);
   gg.resize(ig[nodes_n], 0);
   b.resize(nodes_n, 0);

   for (int i = 0; i < 4; i++) {
      A_loc[i] = new double[4]{};
      M_loc[i] = new double[4]{};
      D[i] = new double[4];
      alphaM[i] = new double[4];
   }

   for (int k = 0; k < el_n; k++) {
      local_el(k);
      for (int i = 0; i < 4; i++) {
         glob_i = el[k].node_n[i];
         for (int j = 0; j < 4; j++) {
            glob_j = el[k].node_n[j];
            if (glob_i > glob_j) {
               i_gg = find_el_pos(glob_i, glob_j);
               gg[i_gg] += A_loc[i][j];
            }
         }
         di[glob_i] += A_loc[i][i];
         b[glob_i] += b_loc[i];
      }
   }
   for (int j = 0; j < face_n; j++) {
      glob_i = face[j];
      di[glob_i] = 1;
      b[glob_i] = val[j];
      for (int i = ig[glob_i]; i < ig[glob_i + 1]; i++) {
         b[jg[i]] -= gg[i] * val[j];
         gg[i] = 0;
      }
      for (int p = glob_i + 1; p < nodes_n; p++)
         for (int i = ig[p]; i < ig[p + 1]; i++) {
            if (jg[i] == glob_i) {
               b[p] -= gg[i] * val[j];
               gg[i] = 0;
            }
         }

   }
   
  /* for (int j = 0; j < face_n; j++) cout << face[j] << endl;
   cout << "di " << ": ";
   for (double node : di)
      cout << node << " ";
   cout << "gg " << ": ";
   for (double node : gg)
      cout << node << " ";
   cout << endl;
   cout << "b " << ": ";
   for (double node : b)
      cout << node << " ";
   cout << endl;*/
}

static void calc_Av(vector<double>& v, vector<double>& res) {
   for (int i = 0; i < nodes_n; i++) {
      res[i] = di[i] * v[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         res[i] += gg[k] * v[j];
         res[j] += gg[k] * v[i];
      }
   }
}
static void calc_r0() {
   for (int i = 0; i < nodes_n; i++) {
      r[i] = b[i] - di[i] * q[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         r[i] -= gg[k] * q[j];
         r[j] -= gg[k] * q[i];
      }
   }
}
static void calc_x() {
   for (int i = 0; i < nodes_n; i++) q[i] += alpha * z[i];
}
static void calc_r(vector<double>& x) {
   for (int i = 0; i < nodes_n; i++) r[i] -= alpha * x[i];
}
static void calc_z(vector<double>& x) {
   for (int i = 0; i < nodes_n; i++) z[i] = x[i] + beta_beta * z[i];
}
static double scMult(vector<double>& x, vector<double>& y) {
   double res = 0;
   for (int i = 0; i < nodes_n; i++) res += x[i] * y[i];
   return res;
}
static vector<double> vecMult(vector<double>& x, vector<double>& y, vector<double>& res) {
   for (int i = 0; i < nodes_n; i++) res[i] = x[i] * y[i];
   return res;
}
static void CGM() {
   double nev = 1;
   q.resize(nodes_n), r.resize(nodes_n), z.resize(nodes_n), Az.resize(nodes_n), Ar.resize(nodes_n), Mr.resize(nodes_n), M.resize(nodes_n);
   for (int i = 0; i < nodes_n; i++) q[i] = 0;  
 
   norma_pr = sqrt(scMult(b, b));
   cout << norma_pr << endl;
   calc_r0();
   double r_scMult = scMult(r, r);
   z = r;
   for (int k = 0; k < maxiter && nev > eps; k++) {
      double r_scMult_old = r_scMult;
      calc_Av(z, Az);
      double zAz_scMult = scMult(Az, z);
      alpha = r_scMult_old / zAz_scMult;
      calc_x();
      calc_r(Az);
      r_scMult = scMult(r, r);
      nev = sqrt(r_scMult) / norma_pr;
      beta_beta = r_scMult / r_scMult_old;
      calc_z(r);
    
   }
   ofstream qFile("q.txt");
  
   cout << "q ";
   for (int i = 0; i < nodes_n; i++) {
      cout << q[i] << " ";
      qFile << scientific << setprecision(10) << q[i] << endl;
   }
   qFile.close();
   r.clear(); z.clear(); Az.clear();
}
vector<double>u;
double u_a(int i) {
   u[i] = sin(nodes[i].x);
   return u[i];
}
vector<double>dif;
void dif_u() {
   dif.resize(nodes_n);
   u.resize(nodes_n);
   ofstream difFile("dif.txt");
   for (int i = 0; i < nodes_n; i++) {
      dif[i] = u_a(i) - q[i];
      difFile << scientific << setprecision(10) << dif[i] << endl;
   }
   difFile.close();
}
int main() {
   int testNumber;

   cout << "Enter the test number: ";
   cin >> testNumber;

 
   input_nodes(testNumber);
   input_el(testNumber);
   input_faces(testNumber);
   input_el_coef(testNumber);
   input_f(testNumber);

   portrait();
   global_A();
   CGM();
   dif_u();
   ofstream uFile("u.txt");
   for (int i = 0; i < nodes_n; i++) {

      uFile << scientific << setprecision(10) << u[i] << endl;
   }
   uFile.close();
   delete[] nodes;  
   delete[] el;
   return 0;
}