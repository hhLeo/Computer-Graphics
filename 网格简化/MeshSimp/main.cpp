#include "Vec3f.h"
#include "SimpleObject.h"
#include "StdAfx.h"
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
using namespace std;
using namespace SimpleOBJ;

CSimpleObject obj;

vector<Pnt> pnts;
vector<Triangle> triangles;
vector<Matrix> pntMatrix, triangleMatrix;
FindSet pntLink;
Heap heap;
vector<vector<int> > pnt2triangle;
map<pair<int, int>, int> edgeHash;

int nPnt, nTriangle, nEdge;
int workTimes = 44000;
double ratio;

//Pnt
Pnt::Pnt() {
	x = y = z = 0;
}

Pnt::Pnt(double x, double y, double z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

void Pnt::unitization(double len) {
	double tmp = sqrt(x * x + y * y + z * z);
	if (tmp < epsilon)	return;
	x = x / tmp * len;
	y = y / tmp * len;
	z = z / tmp * len;
}

double Pnt::dot(const Pnt& B) const {
	return x * B.x + y * B.y + z * B.z;
}

double Pnt::dot(const Pnt& A, const Pnt& B) {
	return A.dot(B);
}

Pnt Pnt::cross(const Pnt& B) const {
	Pnt p;
	p.x = y * B.z - z * B.y;
	p.y = z * B.x - x * B.z;
	p.z = x * B.y - y * B.x;
	return p;
}

Pnt Pnt::cross(const Pnt& A, const Pnt& B) {
	return A.cross(B);
}

void Pnt::plus(const Pnt& B) {
	x += B.x;
	y += B.y;
	z += B.z;
}

Pnt Pnt::plus(const Pnt& A, const Pnt& B) {
	Pnt p = A;
	p.plus(B);
	return p;
}

void Pnt::minus(const Pnt& B) {
	x -= B.x;
	y -= B.y;
	z -= B.z;
}

Pnt Pnt::minus(const Pnt& A, const Pnt& B) {
	Pnt p = A;
	p.minus(B);
	return p;
}

void Pnt::multiply(double u) {
	x *= u;
	y *= u;
	z *= u;
}

Pnt Pnt::multiply(const Pnt &A, double u) {
	Pnt p = A;
	p.multiply(u);
	return p;
}

//Edge
Edge::Edge(int u, int v) {
	if (u < v) {
		this->u = u;
		this->v = v;
	}
	else {
		this->u = v;
		this->v = u;
	}
}

//Heap
Edge& Heap::top() {
	Edge e;
	if (element.size() == 0)	return e;
	return element[0];
}

void Heap::modify(int u) {
	int n = element.size();
	if (n < u + 1)	return;
	while (u > 0) {
		int v = (u - 1) >> 1;
		if (element[v].deltaV > element[u].deltaV) {
			swap(element[v], element[u]);
			edgeMap[element[u].match] = u;
			edgeMap[element[v].match] = v;
		}
		else {
			break;
		}
		
		u = v;
	}

	while (u * 2 + 1 < n) {
		int v = u * 2 + 1, v2 = u * 2 + 2;
		if (v2 < n && element[v].deltaV > element[v2].deltaV) {
			v = v2;
		}
		if (element[v].deltaV < element[u].deltaV) {
			swap(element[v], element[u]);
			edgeMap[element[u].match] = u;
			edgeMap[element[v].match] = v;
		}
		else {
			break;
		}

		u = v;
	}
}

void Heap::del(int u) {
	if (u < 0)	return;
	if (element.size() < u + 1)	return;
	edgeMap[element[u].match] = -1;
	if (element.size() == u + 1) {
		element.pop_back();
	}
	else {
		swap(element[u], element[element.size() - 1]);
		element.pop_back();
		edgeMap[element[u].match] = u;
		modify(u);
	}
}

void Heap::pop() {
	if (element.size()) {
		del(0);
	}
}

void Heap::push(const Edge& e) {
	element.push_back(e);
	if (edgeMap.size() < e.match + 1) {
		edgeMap.resize(e.match + 1);
	}
	edgeMap[e.match] = element.size() - 1;
	modify(element.size() - 1);
}

//FindSet
FindSet::FindSet(int n) {
	f.resize(n);
}

int FindSet::find(int u) {
	if (f[u] == u) {
		return u;
	}
	return f[u] = find(f[u]);
}

void FindSet::add(int u) {
	if (u == -1) {
		f.push_back(f.size());
	}
	else if (f.size() < u + 1) {
		f.resize(u + 1);
		f[u] = u;
	}
	else {
		f[u] = u;
	}
}

//Triangle
Triangle::Triangle(int p0, int p1, int p2) {
	pnt[0] = p0;
	pnt[1] = p1;
	pnt[2] = p2;
}

bool Triangle::valid() {
	bool tmp = true;
	if (pnt[0] == pnt[1] || pnt[1] == pnt[2] || pnt[2] == pnt[0])
		tmp = false;
	return tmp;
}

void Triangle::standard() {
	double tmp = sqrt(a * a + b * b + c * c);
	a /= tmp;
	b /= tmp;
	c /= tmp;
	d /= tmp;
}

//Matrix
Matrix::Matrix() {
	clear();
}

void Matrix::add(const Matrix& m) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			element[i][j] += m.element[i][j];
		}
	}
}

void Matrix::clear() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			element[i][j] = 0;
		}
	}
}

double Matrix::det() {
	double res = 0;
	for (int i = 0; i < 4; ++i) {
		int t[3];
		for (int j = 0, k = 0; j < 4; ++j) {
			if (i != j) {
				t[k++] = j;
			}
		}

		double tmp = element[1][t[0]] * element[2][t[1]] * element[3][t[2]] + element[1][t[1]] * element[2][t[2]] * element[3][t[0]] + element[1][t[2]] * element[2][t[0]] * element[3][t[1]] - element[1][t[0]] * element[2][t[2]] * element[3][t[1]] - element[1][t[1]] * element[2][t[0]] * element[3][t[2]] - element[1][t[2]] * element[2][t[1]] * element[3][t[0]];
		if (i & 1) {
			res -= tmp;
		}
		else {
			res += tmp;
		}
	}
	return res;
}

//Functions
void getTriangle(Triangle &p) {
	for (int i = 0; i < 3; ++i) {
		p.pnt[i] = pntLink.find(p.pnt[i]);
	}
	if (!p.valid()) {
		p.a = p.b = p.c = p.d = 0;
		return;
	}
	Pnt d = Pnt::minus(pnts[p.pnt[1]], pnts[p.pnt[0]]).cross(Pnt::minus(pnts[p.pnt[2]], pnts[p.pnt[0]]));

	p.a = d.x;
	p.b = d.y;
	p.c = d.z;
	p.d = -d.dot(pnts[p.pnt[0]]);

	p.standard();
}

pair<int, int> makeEdgePair(int x, int y) {
	if (x < y) {
		return make_pair(x,y);
	}
	return make_pair(y,x);
}

void getMatrix(const Triangle& p, Matrix& m) {
	const double *t = & p.a;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			m.element[i][j] = t[i] * t[j];
		}
	}
}

double getQ(const Pnt& p, const Matrix& m) {
	double t[4] = {p.x, p.y, p.z, 1}, res = 0;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			res += t[i] * t[j] * m.element[i][j];
		}
	}
	return res;
}

Pnt getNewPnt(const Pnt& u, const Pnt& v, const Matrix& m) {
	Pnt res;
	Matrix m1, m2;
	m1.element[0][0] = m.element[0][0] + m.element[0][0];
	m1.element[0][1] = m.element[0][1] + m.element[1][0];
	m1.element[0][2] = m.element[0][2] + m.element[2][0];
	m1.element[0][3] = m.element[0][3] + m.element[3][0];
	m1.element[1][0] = m.element[0][1] + m.element[1][0];
	m1.element[1][1] = m.element[1][1] + m.element[1][1];
	m1.element[1][2] = m.element[1][2] + m.element[2][1];
	m1.element[1][3] = m.element[1][3] + m.element[3][1];
	m1.element[2][0] = m.element[0][2] + m.element[2][0];
	m1.element[2][1] = m.element[1][2] + m.element[2][1];
	m1.element[2][2] = m.element[2][2] + m.element[2][2];
	m1.element[2][3] = m.element[2][3] + m.element[3][2];
	m1.element[3][3] = 1;

	double det = m1.det();
	if (fabs(det) < epsilon) {
		return Pnt::multiply(Pnt::plus(u, v), 0.5);
	}

	m2 = m1;
	m2.element[0][0] = 0;
	m2.element[1][0] = 0;
	m2.element[2][0] = 0;
	m2.element[3][0] = 1;
	res.x = m2.det() / det;

	m2 = m1;
	m2.element[0][1] = 0;
	m2.element[1][1] = 0;
	m2.element[2][1] = 0;
	m2.element[3][1] = 1;
	res.y = m2.det() / det;

	m2 = m1;
	m2.element[0][2] = 0;
	m2.element[1][2] = 0;
	m2.element[2][2] = 0;
	m2.element[3][2] = 1;
	res.z = m2.det() / det;

	if (getQ(res, m) >  getQ(Pnt::multiply(Pnt::plus(u, v), 0.5), m)) {
		return Pnt::multiply(Pnt::plus(u, v), 0.5);
	}

	return res;
}

void Read() {
	nPnt = obj.m_nVertices;
	nTriangle = obj.m_nTriangles;
	
	pnt2triangle.resize(nPnt);
	pntMatrix.resize(nPnt);
	triangleMatrix.resize(nTriangle);

	for (int i = 0; i < nPnt; ++i) {
		pnts.push_back(Pnt(obj.m_pVertexList[i].x, obj.m_pVertexList[i].y, obj.m_pVertexList[i].z));
		pntLink.add(); //并查集添加
	}

	for (int i = 0; i < nTriangle; ++i) {
		triangles.push_back(Triangle(obj.m_pTriangleList[i][0], obj.m_pTriangleList[i][1], obj.m_pTriangleList[i][2]));
		getTriangle(triangles[i]);
		getMatrix(triangles[i], triangleMatrix[i]);
		for (int j = 0; j < obj.m_pTriangleList[i]._len; ++j) {
			pnt2triangle[obj.m_pTriangleList[i][j]].push_back(i);
			pntMatrix[obj.m_pTriangleList[i][j]].add(triangleMatrix[i]);
		}
	}

	for (int i = 0; i < nTriangle; ++i) {
		for (int j = 0; j < 3; ++j) {
			int k = (j + 1) % 3;
			pair<int, int> e = makeEdgePair(triangles[i].pnt[j], triangles[i].pnt[k]);
			if (edgeHash.count(e) == 0) {
				Matrix m = pntMatrix[e.first];
				m.add(pntMatrix[e.second]);

				edgeHash[e] = nEdge;
				Edge edge(e.first, e.second);
				edge.newPnt = getNewPnt(pnts[e.first], pnts[e.second], m);
				edge.match = nEdge;
				edge.deltaV = getQ(edge.newPnt, m);

				heap.push(edge);

				nEdge++;
			}
		}
	}
	workTimes = int((1.0 - ratio) * (double)nPnt);
}

void Simplify() {
	while (workTimes -- ) {
		Edge edge = heap.top();
		heap.pop();

		//更新信息
		pnts[edge.u] = edge.newPnt; //点u更新成新节点
		pntLink.f[edge.v] = edge.u; //点v就是点u
		
		//更新这两个点相关的所有三角面
		map<int, bool> h, c;
		vector<int> newTriangles;
		for (int tmp = 0, u = edge.u; tmp < 2; ++tmp, u = edge.v) {
			for (int i = 0; i < pnt2triangle[u].size(); ++i) {
				int p = pnt2triangle[u][i];
				if (h.count(p) > 0)	continue;
				h[p] = true;
			
				for (int j = 0; j < 3; ++j) {
					int k = (j + 1) % 3;
					pair<int,int> e = makeEdgePair(triangles[p].pnt[j], triangles[p].pnt[k]);
					if (tmp == 1 && (e.first == u || e.second == u) && edgeHash.count(e) > 0) {
						heap.del(heap.edgeMap[edgeHash[e]]);
						edgeHash.erase(e);
					}
				}
			
				getTriangle(triangles[p]);
				getMatrix(triangles[p], triangleMatrix[p]);

				if (triangles[p].valid()) newTriangles.push_back(p);
			}
		}

		pnt2triangle[edge.u] = newTriangles;
		pnt2triangle[edge.v].clear();
		for (int i = 0; i < newTriangles.size(); ++i)
			for (int j = 0; j < 3; ++j)
				c[triangles[newTriangles[i]].pnt[j]] = true;

		//去除无用的面
		for (map<int, bool>::iterator i = c.begin(); i != c.end(); i++) {
			int x = i->first;
			vector<int> newp;
			for (int j = pnt2triangle[x].size() - 1; j >= 0; --j)
				if (triangles[pnt2triangle[x][j]].valid())	newp.push_back(pnt2triangle[x][j]);
			pnt2triangle[x] = newp;
		}
		//重新计算相关点的矩阵
		for (map<int, bool>::iterator i = c.begin(); i != c.end(); i++) {
			int x = i->first;
			pntMatrix[x].clear();
			for (int j = 0; j < pnt2triangle[x].size(); ++j)
				pntMatrix[x].add(triangleMatrix[pnt2triangle[x][j]]);
		}
		//重新计算相关边的值
		map<pair<int,int>, bool> tmpHash;
		for (map<int, bool>::iterator i = c.begin(); i != c.end(); i++) {
			int x = i->first;
			for (int i = 0; i < pnt2triangle[x].size(); ++i)
				for (int j = 0; j < 3; ++j) {
					int k = (j + 1) % 3;
					pair<int, int> e = makeEdgePair(triangles[pnt2triangle[x][i]].pnt[j], triangles[pnt2triangle[x][i]].pnt[k]);

					if (tmpHash[e] == true) continue;
					tmpHash[e] = true;

					if (edgeHash.count(e) == 0) {
						//未插入的边
						Matrix m = pntMatrix[e.first];
						m.add(pntMatrix[e.second]);

						edgeHash[e] = nEdge;
						Edge edge(e.first, e.second);
						edge.newPnt = getNewPnt(pnts[e.first], pnts[e.second], m);
						edge.match = nEdge;
						edge.deltaV = getQ(edge.newPnt, m);

						heap.push(edge);

						++nEdge;
					}
					else {
						//已插入的边
						Matrix m = pntMatrix[e.first];
						m.add(pntMatrix[e.second]);

						int y = heap.edgeMap[edgeHash[e]];
						heap.element[y].newPnt = getNewPnt(pnts[e.first], pnts[e.second], m);
						heap.element[y].deltaV = getQ(heap.element[y].newPnt, m);

						heap.modify(y);
					}
				}
		}

	}
}

void Write(const char *file) {
	vector<int> rank(nPnt, 0);
	CSimpleObject obj;

	int V = 0, T = 0;
	for (int i = 0; i < nPnt; ++i)
		if (pntLink.find(i) == i)
			rank[i] = V++;
	for (int i = 0; i < nTriangle; ++i)
		if (triangles[i].valid()) ++T;

	obj.m_nTriangles = T;
	obj.m_nVertices = V;
	obj.m_pTriangleList = new Array<int, 3>[T];
	obj.m_pVertexList = new Vec3f[V];

	for (int i = 0, j = 0; i < nPnt; ++i)
		if (pntLink.find(i) == i) {
			obj.m_pVertexList[j] = Vec3f(pnts[i].x, pnts[i].y, pnts[i].z);
			j++;
		}
	for (int i = 0, j = 0; i < nTriangle; ++i)
		if (triangles[i].valid()) {
			for (int k = 0; k < 3; ++k)
				obj.m_pTriangleList[j][k] = rank[triangles[i].pnt[k]];
			j++;
		}

	obj.SaveToObj(file);
}

int main(int argc, char const *argv[]) {
	double t1 = clock();
	ratio = atof(argv[3]);//char* -> double
	obj.LoadFromObj(argv[1]);
	Read();
	Simplify();
	Write(argv[2]);
	double t2 = clock();
	printf("time: %lf\n", (t2 - t1) / 1000.0);
	return 0;
}
