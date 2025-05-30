#include "common_includes.h"

void portrait() {
   map<int, set<int>> list;
   for (int i = 0; i < el_c; ++i) {
      for (int j = 0; j < 8; ++j) {
         int node_index = el[i].node_n[j];
         list[node_index].insert(node_index);
         for (int k = 0; k < 8; ++k) {
            int neighbor_index = el[i].node_n[k];
            if (node_index != neighbor_index) {
               list[node_index].insert(neighbor_index);
            }
         }
      }
   }
   
   // for (size_t i = 0; i < list.size(); ++i) {
   //    cout << "Node " << i << ": ";
   //    for (int node : list[i]) 
   //       cout << node << " ";
   //    cout << endl;
   // }
   int k = 0;
   ig.push_back(0);
   for (int i = 0; i < nodes_c; i++) {
      k = 0;
      for (int node : list[i]) {
         if (node < i) {
            jg.push_back(node);
            k++;
         }
      }
      ig.push_back(ig.back() + k);
   }
   // cout << "jg " << ": ";
   // for (int node : jg) {
   //    cout << node << " ";
   // }
   // cout << endl;
   // cout << "ig " << ": ";
   // for (int node : ig) {
   //    cout << node << " ";
   // }
   // cout << endl;
}

void trans_Gauss(vector<vector<double>> A, int i, int n) {
   int line = i;
   for (int j = i + 1; j < n; j++)
      if (fabs(A[j][i]) > fabs(A[line][i]))
         line = j;
   if (line != i) 
      for (int j = 0; j < 2 * n; j++)
         swap(A[i][j], A[line][j]);
}
void GaussJordan(vector<vector<double>> A, vector<vector<double>> alpha, int n) {
   vector<vector<double>> matrix;
   matrix.resize(n);
   for (int i = 0; i < n; i++) {
      matrix[i].resize(2 * n, 0.0);
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
   for (int i = 0; i < n; i++) matrix[i].clear();
   matrix.clear();
}

static void calc_h() {
	for (int i = 0; i < el_c; i++) {
		int node_min_index = el[i].node_n[0];
		int node_max_index = el[i].node_n[7];

		el[i].hx = nodes[node_max_index].x - nodes[node_min_index].x;
		el[i].hy = nodes[node_max_index].y - nodes[node_min_index].y;
		el[i].hz = nodes[node_max_index].z - nodes[node_min_index].z;
	}
}

static void local_el(int f_el_n) {
	double hx = el[f_el_n].hx;
	double hy = el[f_el_n].hy;
	double hz = el[f_el_n].hz;

	// double л = el[f_el_n].lambda;
	// double л = lam;
	// double г = el[f_el_n].gamma;
	// double г = gam;

	double dx2 = 1 / hx;
	double x2 = hx / 3;

	double dy2 = 1 / hy;
	double y2 = hy / 3;

	double dz2 = 1 / hz;
	double z2 = hz / 3;

	double** gx = nullptr;
	double** gy = nullptr;
	double** gz = nullptr;

   A_loc.resize(8, vector<double>(8));
   b_loc.clear();
   b_loc.resize(8, 0);
   M_loc.resize(8, vector<double>(8));
	gx = new double* [8] ();
	gy = new double* [8] ();
	gz = new double* [8] ();

	for (int i = 0; i < 8; i++) {
		gx[i] = new double[8] ();
		gy[i] = new double[8] ();
		gz[i] = new double[8] ();
	}
	double coefXd = 0, coefYd = 0, coefZd = 0;
	double coefX = 0, coefY = 0, coefZ = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if ((i % 2 != 0 && j % 2 != 0) || (i % 2 == 0 && j % 2 == 0)) coefXd = dx2, coefX = x2;
			else coefXd = -dx2, coefX = x2 / 2;
			if ((i < 4 && j < 4) || (i >= 4 && j >= 4)) coefZd = dz2, coefZ = z2;
			else coefZd = -dz2, coefZ = z2 / 2;
			if (((i == 1 || i == 0 || i == 4 || i == 5) && (j == 1 || j == 0 || j == 4 || j == 5)) || ((i != 1 && i != 0 && i != 4 && i != 5) && (j != 1 && j != 0 && j != 4 && j != 5))) coefYd = dy2, coefY = y2;
			else coefYd = -dy2, coefY = y2 / 2;

			gx[i][j] = coefXd * coefY * coefZ;
			gy[i][j] = coefX * coefYd * coefZ;
			gz[i][j] = coefX * coefY * coefZd;
			M_loc[i][j] = (hx * hy * hz) * coefX * coefY * coefZ; // x 64 less then in excel table???
			A_loc[i][j] = (lam * hx * hy * hz) * (gx[i][j] / (hx * hx) + gy[i][j] / (hy * hy) + gz[i][j] / (hz * hz)) + gam * M_loc[i][j];
			b_loc[i] += nodes[el[f_el_n].node_n[j]].f * M_loc[i][j];
      }
	}
	for (int i = 0; i < 8; i++)
	{
		delete[] gx[i];
		delete[] gy[i];
		delete[] gz[i];
	}
	delete[] gx;
	delete[] gy;
	delete[] gz;
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

   A_loc.resize(8, vector<double>(8, 0));
   M_loc.resize(8, vector<double>(8, 0));
   D.resize(8, vector<double>(8, 0));
   b_loc.resize(8, 0);

   di.resize(nodes_c, 0);
   gg.resize(ig[nodes_c], 0);
   b.resize(nodes_c, 0);


   for (int k = 0; k < el_c; k++) {
      local_el(k);
      for (int i = 0; i < 8; i++) {
         glob_i = el[k].node_n[i];
         for (int j = 0; j < 8; j++) {
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
   for (int j = 0; j < face_c; j++) {
      glob_i = faces[j];
      di[glob_i] = 1;
      b[glob_i] = val[j];
      for (int i = ig[glob_i]; i < ig[glob_i + 1]; i++) {
         b[jg[i]] -= gg[i] * val[j];
         gg[i] = 0;
      }
      for (int p = glob_i + 1; p < nodes_c; p++)
         for (int i = ig[p]; i < ig[p + 1]; i++) {
            if (jg[i] == glob_i) {
               b[p] -= gg[i] * val[j];
               gg[i] = 0;
            }
         }
   }
   // for (int j = 0; j < face_c; j++) cout << faces[j] << endl;
   // cout << "di " << ": ";
   // for (double node : di)
   //    cout << node << " ";
   // cout << "gg " << ": ";
   // for (double node : gg)
   //    cout << node << " ";
   // cout << endl;
   // cout << "b " << ": ";
   // for (double node : b)
   //    cout << node << " ";
   // cout << endl;
}

static void calc_Av(vector<double>& v, vector<double>& res) {
   for (int i = 0; i < nodes_c; i++) {
      res[i] = di[i] * v[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         res[i] += gg[k] * v[j];
         res[j] += gg[k] * v[i];
      }
   }
}
static void calc_r0() {
   for (int i = 0; i < nodes_c; i++) {
      r[i] = b[i] - di[i] * q[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         r[i] -= gg[k] * q[j];
         r[j] -= gg[k] * q[i];
      }
   }
}
static void calc_x() {
   for (int i = 0; i < nodes_c; i++) q[i] += alpha * z[i];
}
static void calc_r(vector<double>& x) {
   for (int i = 0; i < nodes_c; i++) r[i] -= alpha * x[i];
}
static void calc_z(vector<double>& x) {
   for (int i = 0; i < nodes_c; i++) z[i] = x[i] + beta_beta * z[i];
}
static double scMult(vector<double>& x, vector<double>& y) {
   double res = 0;
   for (int i = 0; i < nodes_c; i++) res += x[i] * y[i];
   return res;
}
static vector<double> vecMult(vector<double>& x, vector<double>& y, vector<double>& res) {
   for (int i = 0; i < nodes_c; i++) res[i] = x[i] * y[i];
   return res;
}
static void CGM() {
   double nev = 1;
   q.resize(nodes_c), r.resize(nodes_c), z.resize(nodes_c), Az.resize(nodes_c), Ar.resize(nodes_c), Mr.resize(nodes_c);
   for (int i = 0; i < nodes_c; i++) q[i] = 0;  
 
   norma_pr = sqrt(scMult(b, b));
   // cout << " NORMA " << norma_pr << endl;
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
   ofstream qFile("../q.txt");
  
   // cout << "q ";
   for (int i = 0; i < nodes_c; i++) {
      // cout << q[i] << " ";
      qFile << scientific << setprecision(10) << q[i] << endl;
   }
   qFile.close();
   r.clear(); z.clear(); Az.clear();
}

int main() {
   // eps = 1e-14;   
   // maxiter = 10000;

   int testNumber = 0;

   // cout << "Enter the test number: ";
   // cin >> testNumber;

   input_nodes(testNumber);
   input_el(testNumber);
   input_faces(testNumber);
   input_el_coef(testNumber);
   input_f(testNumber);

   portrait();
   calc_h();
   global_A();
   CGM();
   dif_u();
   print_u();

   return 0;
}