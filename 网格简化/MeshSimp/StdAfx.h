#pragma once
#include <vector>
const double epsilon = 1E-4;

//Pnt
class Pnt
{
public:
	double x, y, z;

	Pnt();
	Pnt(double x, double y, double z);
	~Pnt(){}

	void unitization(double = 1);

	double dot(const Pnt&) const;
	static double dot(const Pnt&, const Pnt&);

	Pnt cross(const Pnt&) const;
	static Pnt cross(const Pnt&, const Pnt&);

	void plus(const Pnt&);
	static Pnt plus(const Pnt&, const Pnt&);

	void minus(const Pnt&);
	static Pnt minus(const Pnt&, const Pnt&);

	void multiply(double);
	static Pnt multiply(const Pnt&, double);
};

//Edge
class Edge
{
public:
	int u, v;
	Pnt newPnt;
	double deltaV;
	int match;

	Edge(){}
	Edge(int u, int v);
	~Edge(){}
};

//Heap
class Heap
{
public:
	std::vector<int> edgeMap;
	std::vector<Edge> element;
	int n;

	Heap(){}
	~Heap(){}

	Edge& top();
	void modify(int u);
	void del(int u);
	void pop();
	void push(const Edge& e);
};

//FindSet
class FindSet
{
public:
	std::vector<int> f;

	FindSet(int n = 0);
	~FindSet(){}

	int find(int u);
	void add(int u = -1);
};

//Triangle
class Triangle
{
public:
	double a, b, c, d;
	int pnt[3];

	Triangle(){}
	Triangle(int p0, int p1, int p2);
	~Triangle(){}
	
	bool valid();
	void standard();
};

//Matrix
class Matrix
{
public:
	double element[4][4];

	Matrix();
	~Matrix(){}

	void add(const Matrix& m);
	void clear();
	double det();
};

