#include "Boundary.h"
#include "MatrixBuilder.h"
#include "Tests.h"

void BlockMatrix(
	vector<vector<double>>& V, vector<vector<double>>& K, vector<vector<double>>& Kt, vector<vector<double>>& D,
	vector<Point>& points, vector<Edge> bound1, vector<Edge> bound2, int elementsCount, int nodeCount)
{
	vector<vector<double>> S(nodeCount + elementsCount, vector<double>(nodeCount + elementsCount));

	//Сборка блочной матрицы из матриц V, K и D
	for (int i = 0; i < elementsCount; i++)
		for (int j = 0; j < elementsCount; j++)
			S[i][j] += V[i][j];

	for (int i = 0; i < elementsCount; i++)
		for (int j = 0; j < nodeCount; j++)
			S[i][j + elementsCount] += -K[i][j];

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < elementsCount; j++)
			S[i + elementsCount][j] += Kt[i][j];

	for (int i = 0; i < nodeCount; i++)
		for (int j = 0; j < elementsCount; j++)
			S[i + elementsCount][j + elementsCount] += D[i][j];

	vector<double> b(nodeCount + elementsCount);

	//Учет краевых условий
	Boundary2(b, bound2, points, true);
	Boundary1(S, b, bound1, points, true);

	vector<double> pq(nodeCount + elementsCount);
	Gauss(S, pq, b);

	for (int i = 0; i < pq.size(); i++)
		cout << pq[i] << endl;
}

void Shur(vector<vector<double>>& V, vector<vector<double>>& K, vector<vector<double>>& D,
	vector<Point>& points, vector<Edge> bound1, vector<Edge> bound2, int elementsCount, int nodeCount)
{
	vector<vector<double>> S(nodeCount, vector<double>(nodeCount));

	ShurComplement(V, K, D, S, elementsCount, nodeCount);
	vector<double> b(nodeCount);

	Boundary2(b, bound2, points, false);
	Boundary1(S, b, bound1, points, false);

	vector<double> q(nodeCount);
	vector<double> res(nodeCount);
	Gauss(S, q, b);

	MultMatrixVector(S, q, res);

	cout << "-------------------- q --------------------" << endl;

	for (int i = 0; i < q.size(); i++)
		cout << q[i] << endl;

	vector<double> Kq(elementsCount);
	MultMatrixVector(K, q, Kq);

	vector<double> p(elementsCount);
	Gauss(V, p, Kq);

	cout << "-------------------- p --------------------" << endl;

	for (int i = 0; i < p.size(); i++)
		cout << p[i] << endl;
}

int main()
{
	// Чтение входных данных сетки
	Interval intervalX, intervalY;
	InputGrid("xy.txt", intervalX, intervalY);

	// Чтение граничных условий
	vector<Edge> bound1, bound2;
	InputBound(bound1, "boundary1.txt");
	InputBound(bound2, "boundary2.txt");

	// Построение сетки по осям x и y
	vector<double> x, y;
	BuildRectMesh(intervalX, x);
	BuildRectMesh(intervalY, y);

	// Создание координат точки
	vector<Point> points;
	CreatePoints(points, x, y);
	int nodeCount = points.size();

	vector<BEMElement> elements;

	// Создание граничных элементов
	CreateRectBEMElements(elements, nodeCount);

	int elementsCount = elements.size();

	vector<vector<double>> V(elementsCount, vector<double>(elementsCount));
	vector<vector<double>> K(elementsCount, vector<double>(nodeCount));
	vector<vector<double>> Kt(nodeCount, vector<double>(elementsCount));
	vector<vector<double>> D(nodeCount, vector<double>(nodeCount));

	vector<vector<double>> M(nodeCount, vector<double>(elementsCount));
	vector<vector<double>> d(nodeCount, vector<double>(1));
	vector<vector<double>> dt(1, vector<double>(nodeCount));
	vector<vector<double>> ddt(nodeCount, vector<double>(nodeCount));

	// Построение отдельных матриц V, K и D
	Build(V, K, Kt, D, points, elements);

	// Создание блочной матрицы
	//BlockMatrix(V, K, Kt, D, points, bound1, bound2, elementsCount, nodeCount);

	// Дополнение Шура
	Shur(V, K, D, points, bound1, bound2, elementsCount, nodeCount);

	// Тесты Дирихле и Неймана
	//DirichletTest(K, V, nodeCount, elementsCount);

	//MatrixKtAndM(K, Kt, M, points, elements);

	//NeumannTest(D, ddt, Kt, M, nodeCount, elementsCount);

	return 0;
}