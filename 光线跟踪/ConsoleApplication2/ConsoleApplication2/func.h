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
	friend Vector3 operator * (const Vector3& , const Vector3&); //cross product���
	friend Vector3& operator *= (Vector3& , const double&);
	friend Vector3& operator *= (Vector3& , const Vector3&);
	friend Vector3 operator / (const Vector3& , const double&);
	friend Vector3& operator /= (Vector3& , const double&);
	double Dot(const Vector3&); //dot product���
	double Module();//ģ��
	double Module2();//ģ����ƽ��
	Vector3 GetUnitVector();//��õ�λ����
	bool IsZeroVector();//�ж��Ƿ�Ϊ��������EPS��
	void Input(std::stringstream& fin);

	Vector3 Reflect(Vector3 N);//����ù���������N�ķ������
	Vector3 Refract(Vector3 N , double n);//����ù���������N��������ߣ�n��ʾ������֮�ȣ�
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
	void Valid(); //ʹɫ�����ȷ���<=1
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
		//H*WΪ�涨��ͼƬ��С
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

	void Initialize(int H , int W);//��ʼ��
	void Input(std::string file);//��bmp��������
	void Output(std::string file);//������ݵ�bmp
	Color GetSmoothColor(double u , double v);//�õ���Ȩƽ�������أ�u��v������ɫ��
};

//Material
//Intersect
//Primitive: Sphere & Plane
class Material {
public:
	Material();
	~Material() {}

	Color color , absor;//ɫ�ʣ�͸����Ʒ���յ�ɫ��
	double refl , refr;//���䡢����[0,1]
	double rindex;//������
	double diff , spec;//�����䡢����������
	double drefl;//���ؿ����㷨ģ�⾵��������ʱ��������ߵ�ƫ��뾶��������0����ʹ�ø��㷨��
	Bmp* texture;//��ΪNULLʱ��������ص�bmp�ļ�������NULLʱ������

	void Input(std::string , std::stringstream&);
};

struct Intersect {
	Vector3 N , C;//��ײ���ķ�ƽ�桢����
	bool front;//�ж��Ƿ�������������ײ
	double distant;//��ײǰ�����߹��ľ���
};

class Primitive {
protected:
	int sample;//ÿ����Ʒ���е�һ���������������Ⱦʱhash�ж��Ƿ�����о��
	Material* material;//����
	Primitive* next;//��Ʒ�����е���һ����Ʒ

public:
	Primitive();
	Primitive(const Primitive&);
	~Primitive();
	
	virtual void Input(std::string , std::stringstream&);
	virtual bool Collide(Vector3 ray_O , Vector3 ray_V) = 0;//�ж��Ƿ���ײ
	Intersect intersect;//�ཻ����Ϣ
	Primitive* GetNext() {
		return next;
	}
	void SetNext(Primitive* primitive) {
		next = primitive;
	}
	Material* GetMaterial() {
		return material;
	}
	virtual Color GetTexture() = 0;//����intersect��������ɫ��
	int GetSample() {
		return sample;
	}
	virtual Primitive* PrimitiveCopy() = 0;//getһ��������Ʒ��ָ��
};

class Sphere : public Primitive {
	Vector3 O , De , Dc;//�����������������ᣨz�����֮��ֱ�ĸ���Ϊ0���ᣩ�����ڼ�������
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
	Vector3 N , Dx , Dy;//ƽ�淨������������
	double R;//ƽ����ԭ��ľ���

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
	Color color;//ɫ��
	Light* next;//��һ����Դ
	int sample;//�жϾ�ݵ������

public:
	Light();
	~Light() {}
	
	double intersect_distant;//����
	virtual void Input(std::string , std::stringstream&);
	virtual bool IsPointLight() = 0;//�ж��Ƿ�Ϊ���Դ
	virtual Vector3 GetO() = 0;
	Light* GetNext() {
		return next;
	}
	void SetNext(Light* light) {
		next = light;
	}
	virtual bool Collide(Vector3 ray_O , Vector3 ray_V) = 0;//�ж��Ƿ��ཻ
	Color GetColor() {
		return color;
	}
	int GetSample() {
		return sample;
	}
	virtual double CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality) = 0;//���ؿ��壬������ܸù�Դ���ߵİٷֱ�
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
	Vector3 O , N , Dx , Dy;//�й���λ�á����򡢾�ͷ�ĳ������
	double lens_W , lens_H;//��ͷ�ĳ�������й�㵽��ͷ����ı�ֵ
	int W , H;//��Ƭ���س����
	Color** data;//�洢ÿ�����ص�ɫ��
	double shade_quality;//~*16=���ؿ��������Ӱʱ�ļ������
	double drefl_quality;//~*16=���ؿ�����㾵��������ʱ�ļ������

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

	Vector3 Emit(double i , double j);//�õ����ض�Ӧ���������
	void Initialize();
	void Input(std::string var , std::stringstream& fin);
	void Output(Bmp*);
};

//Scene
class Scene {
	Primitive* primitive_head;//��Ʒ������ͷ
	Light* light_head;//��Դ������ͷ
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
	Color CalnDiffusion(Primitive* pri , int* hash);//����������ɫ�⣬����pri�Ǹ���ײ������Ʒ��ָ�룬��ײ��Ϣ������intersect�����У�hash���ڼ�¼���߾�������Ʒ��hashֵ��ͬ����������Ҫʹ�ó���ȡ��ȥ��ݣ�
	Color CalnReflection(Primitive* pri , Vector3 ray_V , int dep , int* hash);//���㷴���õ���ɫ�⣬����ray_V��������ߵķ���dep�ǵ�������
	Color CalnRefraction(Primitive* pri , Vector3 ray_V , int dep , int* hash);//���������õ���ɫ��
	Color RayTracing(Vector3 ray_O , Vector3 ray_V , int dep , int* hash);//����׷�ٵ���Ҫ�����������νӲ�ͬ��εĹ���

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
