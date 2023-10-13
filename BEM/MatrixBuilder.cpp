#include "MatrixBuilder.h"

void Build(
   vector<vector<double>>& matrixV, vector<vector<double>>& matrixK, vector<vector<double>>& matrixKt,
   vector<vector<double>>& matrixD, vector<Point>& points, vector<BEMElement>& elements)
{
   for (auto& e1 : elements)
      for (auto& e2 : elements)
      {
         Point a1 = points[e1.verts[0]];
         Point b1 = points[e1.verts[1]];
         Point a2 = points[e2.verts[0]];
         Point b2 = points[e2.verts[1]];

         double h1 = (b1 - a1).Norm();
         double h2 = (b2 - a2).Norm();

         int phiCount = phi.size();
         int psiCount = psi.size();

         vector<vector<QuadratureNode>> nodes(2);

         if (e1 == e2)
         {
            nodes[0] = SegmentGaussOrder23();
            nodes[1] = SegmentGaussOrder21();
         }

         else
         {
            nodes[0] = SegmentGaussOrder23();
            nodes[1] = SegmentGaussOrder23();
         }

         double v = 0.0;
         vector<double> dv(psiCount * psiCount);
         vector<double> k(phiCount * psiCount);

         Point normal1 = Normal(points[e1.verts[0]], points[e1.verts[1]]);
         Point normal2 = Normal(points[e2.verts[0]], points[e2.verts[1]]);

         for (int q1 = 0; q1 < nodes[0].size(); q1++)
         {
            Point p1 = Point(a1.r + nodes[0][q1].Node.r * (b1.r - a1.r), a1.z + nodes[0][q1].Node.r * (b1.z - a1.z));

            for (int q2 = 0; q2 < nodes[1].size(); q2++)
            {
               Point p2 = Point(a2.r + nodes[1][q2].Node.r * (b2.r - a2.r), a2.z + nodes[1][q2].Node.r * (b2.z - a2.z));

               double dr = p1.r - p2.r;
               double dz = p1.z - p2.z;

               Point R = Point(dr, dz);

               double a = p1.r * p1.r + p2.r * p2.r + dz * dz;
               double b = 2.0 * p1.r * p2.r;

               double m = 2.0 * b / (a + b);
               double m1 = 1 - m;

               vector<double> aK = { 1.3862936112, 0.09666344259, 0.03590092383, 0.03742563713, 0.01451196212 };
               vector<double> bK = { 0.5, 0.12498593597, 0.06880248576, 0.03328355346, 0.00441787012 };

               vector<double> aE = { 0.44325141463, 0.06260601220, 0.04757383546, 0.01736506451 };
               vector<double> bE = { 0.24998368310, 0.09200180037, 0.04069697526, 0.00526449639 };

               double K = 0, E = 1.0;

               for (int i = 0; i < 5; i++)
                  K += pow(m1, i) * (aK[i] + bK[i] * log(1 / m1));

               for (int i = 0; i < 4; i++)
                  E += pow(m1, i + 1) * (aE[i] + bE[i] * log(1 / m1));

               double G = K / (M_PI * sqrt(a + b));
               double C = p1.r * p2.r * h1 * h2 * nodes[0][q1].Weight * nodes[1][q2].Weight;

               double dGdr2 = 1.0 / 2.0 / p2.r * ((p1.r * p1.r - p2.r * p2.r + dz * dz) * E / (a - b) - K);
               double dGdr1 = 1.0 / 2.0 / p1.r * ((p2.r * p2.r - p1.r * p1.r + dz * dz) * E / (a - b) - K);
               double dG2dr = dz / 2.0 / p1.r / p2.r * ((p1.r * p1.r + p2.r * p2.r + dz * dz) * E / (a - b) - K);

               v += G * C;

               int ind = 0;

               for (int i = 0; i < phiCount; i++)
                  for (int j = 0; j < psiCount; j++)
                     k[ind++] += phi[i]() * psi[j](nodes[1][q2].Node.r) * C / (M_PI * sqrt(a + b)) *
                     (normal2.r * (E - K) / (2.0 * p2.r) + E * (normal2 * R) / (a - b));

               ind = 0;

               for (int i = 0; i < psiCount; i++)
                  for (int j = 0; j < psiCount; j++)
                     dv[ind++] += -C / (M_PI * sqrt(a + b)) * psi[i](nodes[0][q1].Node.r) *
                     (gradPsi[j]() / h2 * (-normal1.r * dG2dr + normal1.z * dGdr2));
            }
         }

         matrixV[e1.verts[0]][e2.verts[0]] += v;

         matrixD[e1.verts[0]][e2.verts[0]] += dv[0];
         matrixD[e1.verts[0]][e2.verts[1]] += dv[1];
         matrixD[e1.verts[1]][e2.verts[0]] += dv[2];
         matrixD[e1.verts[1]][e2.verts[1]] += dv[3];


         matrixK[e1.verts[0]][e2.verts[0]] += k[0];
         matrixK[e1.verts[0]][e2.verts[1]] += k[1];

         double h = Distance(a1, b1);

         if (e1 == e2)
         {
            for (int i = 0; i < phiCount; i++)
               for (int j = 0; j < psiCount; j++)
                  matrixK[e1.verts[i]][e2.verts[j]] += 0.5 * h * Gauss18([&](double ksi) { return phi[i]() * psi[j](ksi) * (a1.r + ksi * (b1.r - a1.r)); });
         }
      }
}

double* MatrixToVector(int n, int m, vector<vector<double>> vec)
{
   double* matrix = new double[n * m];

   for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
         matrix[i + j * n] = vec[i][j];

   return matrix;
}

vector<vector<double>> VectorToMatrix(int n, int m, double* matrix)
{
   vector<vector<double>> vec(n, vector<double>(m));

   for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
         vec[i][j] = matrix[i + j * n];

   return vec;
}


void ShurComplement(
   vector<vector<double>>& V, vector<vector<double>>& K,
   vector<vector<double>>& D, vector<vector<double>>& S,
   int phiCount, int psiCount)
{
   double* matrixV = MatrixToVector(phiCount, phiCount, V);
   double* matrixK = MatrixToVector(phiCount, psiCount, K);
   double* matrixD = MatrixToVector(psiCount, psiCount, D);

   LAPACKE_dpotrf(CBLAS_LAYOUT::CblasColMajor, 'L', phiCount, matrixV, phiCount); // V = LLT

   vector<vector<double>> temp(phiCount, vector<double>(phiCount));

   cblas_dtrsm(CBLAS_LAYOUT::CblasColMajor, CBLAS_SIDE::CblasLeft, CBLAS_UPLO::CblasLower,
      CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_DIAG::CblasNonUnit, phiCount, psiCount, 1, matrixV,
      phiCount, matrixK, phiCount); // K = L^(-1) * K
   cblas_dsyrk(CBLAS_LAYOUT::CblasColMajor, CBLAS_UPLO::CblasLower, CBLAS_TRANSPOSE::CblasTrans,
      psiCount, phiCount, 1, matrixK, phiCount, 1, matrixD, psiCount); // D = D + K^T * K = D + K^T * V^(-1) * K;

   temp = VectorToMatrix(phiCount, phiCount, matrixD);

   for (int i = 0; i < phiCount; i++)
      for (int j = 0; j < phiCount; j++)
         matrixD[i * phiCount + j] = j <= i ? temp[i][j] : temp[j][i];

   S = VectorToMatrix(psiCount, psiCount, matrixD);
}