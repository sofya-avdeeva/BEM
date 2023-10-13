#pragma once

void MatrixKtAndM(vector<vector<double>> K, vector<vector<double>>& Kt, vector<vector<double>>& M, vector<Point> points, vector<BEMElement> elements)
{
	int elementsCount = elements.size();
	int nodeCount = points.size();

	for (int i = 0; i < elementsCount; i++)
		for (int j = 0; j < nodeCount; j++)
			Kt[i][j] = K[j][i];

	for (auto& el : elements)
	{
		double h = Distance(points[el.verts[0]], points[el.verts[1]]);

		for (int i = 0; i < psi.size(); i++)
			for (int j = 0; j < phi.size(); j++)
				M[el.verts[i]][el.verts[j]] += h * Gauss18([&](double ksi) {return phi[j]() * psi[i](ksi) * (points[el.verts[0]].r + ksi * (points[el.verts[1]].r - points[el.verts[0]].r)); });
	
	}
}

void DirichletTest(vector<vector<double>> K, vector<vector<double>> V, int nodeCount, int elementsCount)
{
	vector<double> q(nodeCount);
	vector<double> p(nodeCount);

	q[0] = 0.1;
	q[1] = 0.1;
	q[2] = 0.1;
	q[3] = 0.2;
	q[4] = 0.3;
	q[5] = 0.3;
	q[6] = 0.3;
	q[7] = 0.2;

	p[0] = -1;
	p[1] = -1;
	p[2] = 0;
	p[3] = 0;
	p[4] = 1;
	p[5] = 1;
	p[6] = 0;
	p[7] = 0;

	vector<double> Vp(nodeCount);
	vector<double> Kq(nodeCount);

	MultMatrixVector(K, q, Kq);
	MultMatrixVector(V, p, Vp);

	cout << "Kq: " << endl;

	for (int i = 0; i < Kq.size(); i++)
		cout << Kq[i] << endl;

	cout << endl;

	cout << "Vp: " << endl;

	for (int i = 0; i < Vp.size(); i++)
		cout << Vp[i] << endl;

	Gauss(V, p, Kq);
	//Gauss(K, q, Vp);
	for (int i = 0; i < p.size(); i++)
		cout << p[i] << endl;
}

void NeumannTest(vector<vector<double>> D, vector<vector<double>> ddt, vector<vector<double>> Kt, vector<vector<double>> M, int nodeCount, int elementsCount)
{
	vector<double> q(nodeCount);

	q[0] =0.1;
	q[1] =0.1;
	q[2] =0.1;
	q[3] =0.2;
	q[4] =0.3;
	q[5] =0.3;
	q[6] =0.3;
	q[7] =0.2;

	vector<double> Dq(nodeCount);

	//for (int i = 0; i < D.size(); i++)
	//	for (int j = 0; j < D[i].size(); j++)
	//		D[i][j] = D[i][j] + ddt[i][j];

	MultMatrixVector(D, q, Dq);

	cout << "Dq: " << endl;

	for (int i = 0; i < Dq.size(); i++)
		cout << Dq[i] << endl;

	cout << endl;

	vector<double> p(elementsCount);

	p[0] = -1;
	p[1] = -1;
	p[2] = 0;
	p[3] = 0;
	p[4] = 1;
	p[5] = 1;
	p[6] = 0;
	p[7] = 0;

	vector<double> Ktp(nodeCount);

	//for (int i = 0; i < M.size(); i++)
	//	for (int j = 0; j < M[i].size(); j++)
	//		Kt[i][j] = M[i][j] - Kt[i][j];

	MultMatrixVector(Kt, p, Ktp);

	cout << "(M-Kt)p: " << endl;

	for (int i = 0; i < Ktp.size(); i++)
		cout << Ktp[i] << endl;

	cout << endl << endl;

	Gauss(Kt, p, Dq);
	//Gauss(K, q, Vp);
	for (int i = 0; i < p.size(); i++)
		cout << p[i] << endl;
}