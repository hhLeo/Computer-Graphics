#pragma once

#define EPS 1e-6
#define PI 3.14159265359
#define ran() (double(rand() % 32768) / 32768)

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <sstream>
#include <ctime>

using namespace std;

//Vector3
class Vector3 {
public:
	double x , y , z;
	
	Vector3(double X = 0 , double Y = 0 , double Z = 0) : x(X) , y(Y) , z(Z) {}
	~Vector3() {}

	friend Vector3 operator + (const Vector3& , const Vector3&);
	friend Vector3& operator += (Vector3& , const Vector3&);
	friend Vector3 operator - (const Vector3& , const Vector3&);
	friend Vector3& operator -= (Vector3& , const Vector3&);
	friend Vector3 operator - (const Vector3&);
	friend Vector3 operator * (const Vector3& , const double&);
	friend Vector3 operator * (const Vector3& , const Vector3&); //cross product叉乘
	friend Vector3& operator *= (Vector3& , const double&);
	friend Vector3& operator *= (Vector3& , const Vector3&);
	friend Vector3 operator / (const Vector3& , const double&);
	friend Vector3& operator /= (Vector3& , const double&);
	double Dot(const Vector3&); //dot product点乘
	double Module();//模长
	double Module2();//模长的平方
	Vector3 GetUnitVector();//获得单位向量
	bool IsZeroVector();//判断是否为零向量（EPS）
	void Input(std::stringstream& fin);

	Vector3 Reflect(Vector3 N);//计算该光线遇到面N的反射光线
	Vector3 Refract(Vector3 N , double n);//计算该光线遇到面N的折射光线（n表示折射率之比）
};

//Color
class Color {
public:
	double r , g , b;

	Color(double R = 0 , double G = 0 , double B = 0) : r(R) , g(G) , b(B) {}
	~Color() {}

	friend Color operator + (const Color& , const Color&);
	friend Color& operator += (Color& , const Color&);
	friend Color operator - (const Color& , const Color&);
	friend Color& operator -= (Color& , const Color&);
	friend Color operator * (const Color& , const Color&);
	friend Color operator * (const Color& , const double&);
	friend Color& operator *= (Color& , const double&);
	friend Color operator / (const Color& , const double&);
	friend Color& operator /= (Color& , const double&);
	void Valid(); //使色光亮度分量<=1
	void Input(std::stringstream&);
};

//Bmp
typedef unsigned char byte;
typedef unsigned short word;
typedef unsigned int dword;

struct BITMAPFILEHEADER {
	dword bfSize;
	word bfReserved1;
	word bfReserved2;
	dword bfOffBits;
};

struct BITMAPINFOHEADER {
	dword biSize;
	long biWidth;
	long biHeight;
	word biPlanes;
	word biBitCount;
	dword biCompression;
	dword biSizeImage;
	long biXPelsPerMeter;
	long biYPelsPerMeter;
	dword biClrUsed;
	dword biClrImportant;
};

struct RGBQUAD {
	byte rgbBlue;
	byte rgbGreen;
	byte rgbRed;
	byte rgbReserved;
};

struct IMAGEDATA {
	byte red;
	byte green;
	byte blue;
	Color GetColor() {
		return Color(red , green , blue) / 256;
	}
};

class Bmp {
	BITMAPFILEHEADER strHead;
	BITMAPINFOHEADER strInfo;
	bool ima_created;
	IMAGEDATA** ima;

	void Release();
	
public:
	Bmp(int H = 0 , int W = 0) {
		//H*W为规定的图片大小
		Initialize(H , W);
	}
	~Bmp() {
		Release();
	}

	int GetH() {
		return strInfo.biHeight;
	}
	int GetW() {
		return strInfo.biWidth;
	}
	Color GetColor(int i , int j) {
		return Color(ima[i][j].red , ima[i][j].green , ima[i][j].blue) / 256;
	}
	void SetColor(int i , int j , Color);

	void Initialize(int H , int W);//初始化
	void Input(std::string file);//从bmp读入数据
	void Output(std::string file);//输出数据到bmp
	Color GetSmoothColor(double u , double v);//得到加权平均的像素（u，v）处的色光
};

//Material
//Intersect
//Primitive: Sphere & Plane
class Material {
public:
	Material();
	~Material() {}

	Color color , absor;//色彩，透明物品吸收的色光
	double refl , refr;//反射、折射[0,1]
	double rindex;//折射率
	double diff , spec;//漫反射、镜面漫反射
	double drefl;//蒙特卡洛算法模拟镜面漫反射时，反射光线的偏差半径（若等于0，则不使用该算法）
	Bmp* texture;//不为NULL时，纹理相关的bmp文件；等于NULL时无纹理

	void Input(std::string , std::stringstream&);
};

struct Intersect {
	Vector3 N , C;//碰撞处的法平面、坐标
	bool front;//判断是否从物体的正面碰撞
	double distant;//碰撞前光线走过的距离
};

class Primitive {
protected:
	int sample;//每个物品都有的一个随机数，用于渲染时hash判断是否可能有锯齿
	Material* material;//材质
	Primitive* next;//物品链表中的下一个物品

public:
	Primitive();
	Primitive(const Primitive&);
	~Primitive();
	
	virtual void Input(std::string , std::stringstream&);
	virtual bool Collide(Vector3 ray_O , Vector3 ray_V) = 0;//判断是否碰撞
	Intersect intersect;//相交的信息
	Primitive* GetNext() {
		return next;
	}
	void SetNext(Primitive* primitive) {
		next = primitive;
	}
	Material* GetMaterial() {
		return material;
	}
	virtual Color GetTexture() = 0;//根据intersect计算纹理色彩
	int GetSample() {
		return sample;
	}
	virtual Primitive* PrimitiveCopy() = 0;//get一个复制物品的指针
};

class Sphere : public Primitive {
	Vector3 O , De , Dc;//球心坐标和球的坐标轴（z轴和与之垂直的辅角为0的轴），用于计算纹理
	double R;

public:
	Sphere();
	~Sphere() {}

	void Input(std::string , std::stringstream&);
	bool Collide(Vector3 ray_O , Vector3 ray_V);
	Primitive* PrimitiveCopy();
	Color GetTexture();
};

class Plane : public Primitive {
	Vector3 N , Dx , Dy;//平面法向量、坐标轴
	double R;//平面与原点的距离

public:
	Plane() : Primitive() {}
	~Plane() {}

	void Input(std::string , std::stringstream&);
	bool Collide(Vector3 ray_O , Vector3 ray_V);
	Primitive* PrimitiveCopy();
	Color GetTexture();
};


//Light: PointLight & AreaLight
class Light {
protected:
	Color color;//色光
	Light* next;//下一个光源
	int sample;//判断锯齿的随机数

public:
	Light();
	~Light() {}
	
	double intersect_distant;//距离
	virtual void Input(std::string , std::stringstream&);
	virtual bool IsPointLight() = 0;//判断是否为点光源
	virtual Vector3 GetO() = 0;
	Light* GetNext() {
		return next;
	}
	void SetNext(Light* light) {
		next = light;
	}
	virtual bool Collide(Vector3 ray_O , Vector3 ray_V) = 0;//判断是否相交
	Color GetColor() {
		return color;
	}
	int GetSample() {
		return sample;
	}
	virtual double CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality) = 0;//蒙特卡洛，计算接受该光源光线的百分比
};

class PointLight : public Light {
	Vector3 O;
public:
	PointLight() : Light() {}
	~PointLight() {}
	
	void Input(std::string , std::stringstream&);
	bool IsPointLight() {return true;}
	Vector3 GetO() {return O;}
	bool Collide(Vector3 ray_O , Vector3 ray_V);
	double CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality);
};

class AreaLight : public Light {
	Vector3 O , Dx , Dy;
public:
	AreaLight() : Light() {}
	~AreaLight() {}
	
	bool IsPointLight() {return false;}
	Vector3 GetO() {return O;}
	void Input(std::string , std::stringstream&);
	bool Collide(Vector3 ray_O , Vector3 ray_V);
	double CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality);
};

//Camera
extern const double STD_LENS_WIDTH; //the width of lens in the scene
extern const double STD_LENS_HEIGHT;
extern const int STD_IMAGE_WIDTH;
extern const int STD_IMAGE_HEIGTH;
extern const int STD_SHADE_QUALITY; //caln shade :: how many times will be run (*16)
extern const int STD_DREFL_QUALITY; //caln drefl :: how many times will be run (*16)

class Camera {
	Vector3 O , N , Dx , Dy;//感光点的位置、朝向、镜头的长宽半轴
	double lens_W , lens_H;//镜头的长、宽与感光点到镜头距离的比值
	int W , H;//照片像素长与宽
	Color** data;//存储每个像素的色彩
	double shade_quality;//~*16=蒙特卡洛计算阴影时的计算次数
	double drefl_quality;//~*16=蒙特卡洛计算镜面漫反射时的计算次数

public:
	Camera();
	~Camera();
	
	Vector3 GetO() {
		return O;
	}
	int GetW() {
		return W;
	}
	int GetH() {
		return H;
	}
	void SetColor(int i , int j , Color color) {
		data[i][j] = color;
	}
	double GetDreflQuality() {
		return drefl_quality;
	}
	double GetShadeQuality() {
		return shade_quality;
	}

	Vector3 Emit(double i , double j);//得到像素对应的射出光线
	void Initialize();
	void Input(std::string var , std::stringstream& fin);
	void Output(Bmp*);
};

//Scene
class Scene {
	Primitive* primitive_head;//物品的链表头
	Light* light_head;//光源的链表头
	Color background_color;
	Camera* camera;

public:
	Scene();
	~Scene();
	
	Color GetBackgroundColor() {
		return background_color;
	}
	Primitive* GetPrimitiveHead() {
		return primitive_head;
	}
	Light* GetLightHead() {
		return light_head;
	}
	Camera* GetCamera() {
		return camera;
	}

	Primitive* FindNearestPrimitive(Vector3 ray_O , Vector3 ray_V);
	Light* FindNearestLight(Vector3 ray_O , Vector3 ray_V);
	void CreateScene(std::string file);
};

//Raytracer
extern const double SPEC_POWER;
extern const int MAX_DREFL_DEP;
extern const int MAX_RAYTRACING_DEP;
extern const int HASH_FAC;
extern const int HASH_MOD;

class Raytracer {
	std::string input , output;
	Scene scene;
	Color CalnDiffusion(Primitive* pri , int* hash);//计算漫反射色光，其中pri是刚碰撞到的物品的指针，碰撞信息储存在intersect变量中，hash用于记录光线经过的物品（hash值不同的相邻象素要使用超级取样去锯齿）
	Color CalnReflection(Primitive* pri , Vector3 ray_V , int dep , int* hash);//计算反射后得到的色光，其中ray_V是入射光线的方向，dep是迭代层数
	Color CalnRefraction(Primitive* pri , Vector3 ray_V , int dep , int* hash);//计算折射后得到的色光
	Color RayTracing(Vector3 ray_O , Vector3 ray_V , int dep , int* hash);//光线追踪的主要迭代函数，衔接不同层次的光线

public:
	Raytracer() {}
	~Raytracer() {}
	
	void SetInput(std::string file) {
		input = file;
	}
	void SetOutput(std::string file) {
		output = file;
	}
	void Run();
};
