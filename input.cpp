#include "common_includes.h"

void input_nodes(int testNumber) {
   char filename[50];
   sprintf_s(filename, "nodes%d.txt", testNumber);
   FILE* nodes_f;
   if (fopen_s(&nodes_f, filename, "r") != 0) {
      perror("No such file");
      return;
   }
   fscanf_s(nodes_f, "%d", &nodes_n);
   nodes = new node[nodes_n]; 
   for (int i = 0; i < nodes_n; i++) {
      fscanf_s(nodes_f, "%lf %lf %lf", &nodes[i].x, &nodes[i].y, &nodes[i].z);
      nodes[i].number = i;
   }
   fclose(nodes_f);
}
void input_el(int testNumber) {
   char filename[50];
   sprintf_s(filename, "el%d.txt", testNumber);

   FILE* el_f;
   if (fopen_s(&el_f, filename, "r") != 0) {
      perror("No such file");
      return;
   }
   fscanf_s(el_f, "%d", &el_n);
   el = new EL[el_n]; 
   for (int i = 0; i < el_n; i++) {
      fscanf_s(el_f, "%d %d %d %d %d %d %d %d",
         &el[i].node_n[0],
         &el[i].node_n[1],
         &el[i].node_n[2],
         &el[i].node_n[3],
         &el[i].node_n[4],
         &el[i].node_n[5],
         &el[i].node_n[6],
         &el[i].node_n[7]
      );
      el[i].number = i;
   }
   fclose(el_f);
   //for (int i = 0; i < el_n; i++) {
   //   cout << "Element " << i << ": Nodes "
   //      << el[i].node_n[0] << ", "
   //      << el[i].node_n[1] << ", "
   //      << el[i].node_n[2] << ", "
   //      << el[i].node_n[3] << endl;
   //}
}
void input_faces(int testNumber) {
   char filename[50];
   sprintf_s(filename, "face%d.txt", testNumber);

   FILE* face_f;
   if (fopen_s(&face_f, filename, "r") != 0) {
      perror("No such file");
      return;
   }
   fscanf_s(face_f, "%d", &face_n);
   face.resize(face_n);
   val.resize(face_n);
   for (int i = 0; i < face_n; i++) {
      fscanf_s(face_f, "%d", &face[i]);
      fscanf_s(face_f, "%lf", &val[i]);
   }
   fclose(face_f);
}
void input_el_coef(int testNumber) {
   char filename[50];
   sprintf_s(filename, "coef%d.txt", testNumber); 

   FILE* el_coef_f;
   if (fopen_s(&el_coef_f, filename, "r") != 0) {
      perror("No such file");
      return;
   }
   for (int i = 0; i < el_n; i++) {
      fscanf_s(el_coef_f, "%lf %lf", &el[i].lambda, &el[i].gamma);
      el[i].number = i;
   }
   fclose(el_coef_f);
}
void input_f(int testNumber) {
   char filename[50];
   sprintf_s(filename, "f%d.txt", testNumber);

   FILE* f_f;
   if (fopen_s(&f_f, filename, "r") != 0) {
      perror("No such file");
      return;
   }
   for (int i = 0; i < nodes_n; i++) 
      fscanf_s(f_f, "%lf", &nodes[i].f);
   fclose(f_f);
} 
