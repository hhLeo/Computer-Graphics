#include "func.h"

void kitten(){
	//Plane
	//Texture

	//plane->material = new PhongMaterial
	//light

	//Scene *scene = new Scene
	
	//scene->addLight
	//scene->addMesh
	
	//Camera

	//RayTracer rayTracer

	//img = rayTrace(scene, camera, width, height);

	//imwrite("kitten.png", img);
}

//Vector3
Vector3 operator + (const Vector3& A , const Vector3& B) {
	return Vector3(A.x + B.x , A.y + B.y , A.z + B.z);
}

Vector3 operator - (const Vector3& A , const Vector3& B) {
	return Vector3(A.x - B.x , A.y - B.y , A.z - B.z);
}

Vector3 operator * (const Vector3& A , const double& k) {
	return Vector3(A.x * k , A.y * k , A.z * k);
}

Vector3 operator / (const Vector3& A , const double& k) {
	return Vector3(A.x / k , A.y / k , A.z / k);
}

Vector3 operator * (const Vector3& A , const Vector3& B) {
	return Vector3(A.y * B.z - A.z * B.y , A.z * B.x - A.x * B.z , A.x * B.y - A.y * B.x);
}

Vector3& operator += (Vector3& A , const Vector3& B) {
	A = A + B;
	return A;
}

Vector3& operator -= (Vector3& A , const Vector3& B) {
	A = A - B;
	return A;
}

Vector3& operator *= (Vector3& A , const double& k) {
	A = A * k;
	return A;
}

Vector3& operator += (Vector3& A , const double& k) {
	A = A / k;
	return A;
}

Vector3& operator *= (Vector3& A , const Vector3& B) {
	A = A * B;
	return A;
}

Vector3 operator - (const Vector3& A) {
	return Vector3(-A.x , -A.y , -A.z);
}

double Vector3::Dot(const Vector3& term) {
	return x * term.x + y * term.y + z * term.z;
}

double Vector3::Module2() {
	return x * x + y * y + z * z;
}

double Vector3::Module() {
	return sqrt(x * x + y * y + z * z);
}

Vector3 Vector3::GetUnitVector() {
	return *this / Module();
}

bool Vector3::IsZeroVector() {
	return fabs(x) < EPS && fabs(y) < EPS && fabs(z) < EPS;
}

void Vector3::Input(std::stringstream& fin) {
	fin >> x >> y >> z;
}

Vector3 Vector3::Reflect(Vector3 N) {
	return *this - N * (2 * Dot(N));
}

Vector3 Vector3::Refract(Vector3 N , double n) {
	Vector3 V = GetUnitVector();
	double cosI = -N.Dot(V) , cosT2 = 1 - (n * n) * (1 - cosI * cosI); 
	if (cosT2 > EPS) return V * n + N * (n * cosI - sqrt(cosT2));
	return V.Reflect(N);
}

//Color
Color operator + (const Color& A , const Color& B) {
	return Color(A.r + B.r , A.g + B.g , A.b + B.b);
}

Color operator - (const Color& A , const Color& B) {
	return Color(A.r - B.r , A.g - B.g , A.b - B.b);
}

Color operator * (const Color& A , const Color& B) {
	return Color(A.r * B.r , A.g * B.g , A.b * B.b);
}

Color operator * (const Color& A , const double& k) {
	return Color(A.r * k , A.g * k , A.b * k);
}

Color operator / (const Color& A , const double& k) {
	return Color(A.r / k , A.g / k , A.b / k);
}

Color& operator += (Color& A , const Color& B) {
	A = A + B;
	return A;
}

Color& operator -= (Color& A , const Color& B) {
	A = A - B;
	return A;
}

Color& operator *= (Color& A , const double& k) {
	A = A * k;
	return A;
}

Color& operator /= (Color& A , const double& k) {
	A = A / k;
	return A;
}

//使A的色光亮度分量不大于1（实际上有可能大于1但屏幕不能显示）
void Color::Valid() {
	if (r > 1) r = 1;
	if (g > 1) g = 1;
	if (b > 1) b = 1;
}

void Color::Input(std::stringstream& fin) {
	fin >> r >> g >> b;
}

//Bmp
void Bmp::Initialize(int H , int W) {
	strHead.bfReserved1 = 0;
	strHead.bfReserved2 = 0;
	strHead.bfOffBits = 54;

	strInfo.biSize = 40;
	strInfo.biPlanes = 1;
	strInfo.biHeight = H;
	strInfo.biWidth = W;
	strInfo.biBitCount = 24;
	strInfo.biCompression = 0;
	strInfo.biSizeImage = H * W * 3;
	strInfo.biXPelsPerMeter = 0;
	strInfo.biYPelsPerMeter = 0;
	strInfo.biClrUsed = 0;
	strInfo.biClrImportant = 0;

	strHead.bfSize = strInfo.biSizeImage + strInfo.biBitCount;
	
	ima = new IMAGEDATA*[H];
	for (int i = 0 ; i < H ; i++)
		ima[i] = new IMAGEDATA[W];
}

void Bmp::Release() {
	for (int i = 0 ; i < strInfo.biHeight ; i++)
		delete[] ima[i];

	delete[] ima;
}

void Bmp::Input(std::string file) {
	Release();

	FILE *fpi = fopen(file.c_str() , "rb");
	word bfType;
	fread(&bfType , 1 , sizeof(word) , fpi);
	fread(&strHead , 1 , sizeof(BITMAPFILEHEADER) , fpi);
	fread(&strInfo , 1 , sizeof(BITMAPINFOHEADER) , fpi);
	
	RGBQUAD Pla;
	for (int i = 0 ; i < (int) strInfo.biClrUsed ; i++) {
		fread((char *) & (Pla.rgbBlue) , 1 , sizeof(byte) , fpi);
		fread((char *) & (Pla.rgbGreen) , 1 , sizeof(byte) , fpi);
		fread((char *) & (Pla.rgbRed) , 1 , sizeof(byte) , fpi);
	}
	
	Initialize(strInfo.biHeight , strInfo.biWidth);
	for(int i = 0 ; i < strInfo.biHeight ; i++)
		for(int j = 0 ; j < strInfo.biWidth ; j++) {
			fread(&ima[i][j].blue , 1 , sizeof(byte) , fpi);
			fread(&ima[i][j].green , 1 , sizeof(byte) , fpi);
			fread(&ima[i][j].red , 1 , sizeof(byte) , fpi);
		}

	fclose(fpi);
}

void Bmp::Output(std::string file) {
	FILE *fpw = fopen(file.c_str() , "wb");

	word bfType = 0x4d42;
	fwrite(&bfType , 1 , sizeof(word) , fpw);
	fwrite(&strHead , 1 , sizeof(BITMAPFILEHEADER) , fpw);
	fwrite(&strInfo , 1 , sizeof(BITMAPINFOHEADER) , fpw);

	for (int i = 0 ; i < strInfo.biHeight ; i++)
		for (int j = 0 ; j < strInfo.biWidth ; j++) {
			fwrite(&ima[i][j].blue , 1 , sizeof(byte) , fpw);
			fwrite(&ima[i][j].green , 1 , sizeof(byte) , fpw);
			fwrite(&ima[i][j].red , 1 , sizeof(byte) , fpw);
		}
	
	fclose(fpw);
}

void Bmp::SetColor(int i , int j , Color col) {
	ima[i][j].red = (int) (col.r * 255);
	ima[i][j].green = (int) (col.g * 255);
	ima[i][j].blue = (int) (col.b * 255);
}

//得到在象素(u,v)的色光，取周围整数象素色光的加权平均值
Color Bmp::GetSmoothColor(double u , double v) {
	double U = (u - floor(u)) * strInfo.biHeight;
	double V = (v - floor(v)) * strInfo.biWidth;
	int U1 = (int) floor(U + EPS) , U2 = U1 + 1;
	int V1 = (int) floor(V + EPS) , V2 = V1 + 1;
	double rat_U = U2 - U;
	double rat_V = V2 - V;
	if (U1 < 0) U1 = strInfo.biHeight - 1;
	if (U2 == strInfo.biHeight) U2 = 0;
	if (V1 < 0) V1 = strInfo.biWidth - 1;
	if (V2 == strInfo.biWidth) V2 = 0;
	Color ret;
	ret += ima[U1][V1].GetColor() * rat_U * rat_V;
	ret += ima[U1][V2].GetColor() * rat_U * (1 - rat_V);
	ret += ima[U2][V1].GetColor() * (1 - rat_U) * rat_V;
	ret += ima[U2][V2].GetColor() * (1 - rat_U) * (1 - rat_V);
	return ret;
}

//Material
//Intersect
//Primitive: Sphere & Plane
Material::Material() {
	color = absor = Color();
	refl = refr = 0;
	diff = spec = 0;
	rindex = 0;
	drefl = 0;
	texture = NULL;
}

void Material::Input(std::string var , std::stringstream& fin) {
	if (var == "color=") color.Input(fin);
	if (var == "absor=") absor.Input(fin);
	if (var == "refl=") fin >> refl;
	if (var == "refr=") fin >> refr;
	if (var == "diff=") fin >> diff;
	if (var == "spec=") fin >> spec;
	if (var == "drefl=") fin >> drefl;
	if (var == "rindex=") fin >> rindex;
	if (var == "texture=") {
		std::string file; fin >> file;
		texture = new Bmp;
		texture->Input(file);
	}
}

Primitive::Primitive() {
	sample = rand();
	material = new Material;
	next = NULL;
}

Primitive::Primitive(const Primitive& primitive) {
	*this = primitive;
	material = new Material;
	*material = *primitive.material;
}

Primitive::~Primitive() {
	delete material;
}

void Primitive::Input(std::string var , std::stringstream& fin) {
	material->Input(var , fin);
}

Sphere::Sphere() : Primitive() {
	De = Vector3(0 , 0 , 1);
	Dc = Vector3(0 , 1 , 0);
}

void Sphere::Input(std::string var , std::stringstream& fin) {
	if (var == "O=") O.Input(fin);
	if (var == "R=") fin >> R;
	if (var == "De=") De.Input(fin);
	if (var == "Dc=") Dc.Input(fin);
	Primitive::Input(var , fin);
}

//计算光线(ray_O,ray_V)到物品碰撞情况，碰撞数据存在intersect中，返回是否有碰撞
bool Sphere::Collide(Vector3 ray_O , Vector3 ray_V) {
	ray_V = ray_V.GetUnitVector();
	Vector3 P = ray_O - O;
	double b = -P.Dot(ray_V);
	double det = b * b - P.Module2() + R * R;

	if (det > EPS) {
		det = sqrt(det);
		double x1 = b - det  , x2 = b + det;

		if (x2 < EPS) return false;
		if (x1 > EPS) {
			intersect.distant = x1;
			intersect.front = true;
		} else {
			intersect.distant = x2;
			intersect.front = false;
		} 
	} else 
		return false;

	intersect.C = ray_O + ray_V * intersect.distant;
	intersect.N = (intersect.C - O).GetUnitVector();
	if (intersect.front == false) intersect.N = -intersect.N;
	return true;
}

//根据碰撞信息intersect计算碰撞点的纹理色彩
Color Sphere::GetTexture() {
	Vector3 I = (intersect.C - O).GetUnitVector();
	double a = acos(-I.Dot(De));
	double b = acos(std::min(std::max(I.Dot(Dc) / sin(a) , -1.0) , 1.0));
	double u = a / PI , v = b / 2 / PI;
	if (I.Dot(Dc * De) < 0) v = 1 - v;
	return material->texture->GetSmoothColor(u , v);
}

//得到一个复制的物品的指针
Primitive* Sphere::PrimitiveCopy() {
	Sphere* ret = new Sphere(*this);
	return ret;
}

void Plane::Input(std::string var , std::stringstream& fin) {
	if (var == "N=") N.Input(fin);
	if (var == "R=") fin >> R;
	if (var == "Dx=") Dx.Input(fin);
	if (var == "Dy=") Dy.Input(fin);
	Primitive::Input(var , fin);
}

bool Plane::Collide(Vector3 ray_O , Vector3 ray_V) {
	ray_V = ray_V.GetUnitVector();
	N = N.GetUnitVector();
	double d = N.Dot(ray_V);
	if (fabs(d) < EPS) return false;
	double l = (N * R - ray_O).Dot(N) / d;
	if (l < EPS) return false;

	intersect.distant = l;
	intersect.front = (d < 0);
	intersect.C = ray_O + ray_V * intersect.distant;
	intersect.N = (intersect.front) ? N : -N;
	return true;
}

Color Plane::GetTexture() {
	double u = intersect.C.Dot(Dx) / Dx.Module2();
	double v = intersect.C.Dot(Dy) / Dy.Module2();
	return material->texture->GetSmoothColor(u , v);
}

Primitive* Plane::PrimitiveCopy() {
	Plane* ret = new Plane(*this);
	return ret;
}


//Light: PointLight & AreaLight
Light::Light() {
	sample = rand();
	next = NULL;
}

void Light::Input(std::string var , std::stringstream& fin) {
	if (var == "color=") 
		color.Input(fin);
}

void PointLight::Input(std::string var , std::stringstream& fin) {
	if (var == "O=")
		O.Input(fin);
	Light::Input(var , fin);
}

//计算光线(ray_O,ray_V)到物品碰撞情况，碰撞前走过距离存在intersect_distant中，返回是否有碰撞
bool PointLight::Collide(Vector3 ray_O , Vector3 ray_V) {
	return false;
}


//计算点C的接受到该光源光线的百分比，蒙特卡罗算法，shade_quality*16是计算次数
double PointLight::CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality) {
	Vector3 V = O - C;
	double distant = V.Module();

	for (Primitive* now = primitive_head ; now != NULL ; now = now->GetNext())
		if (now->Collide(C , V) && (now->intersect.distant < distant))
			return 0;

	return 1;
}

void AreaLight::Input(std::string var , std::stringstream& fin) {
	if (var == "O=") O.Input(fin);
	if (var == "Dx=") Dx.Input(fin);
	if (var == "Dy=") Dy.Input(fin);
	Light::Input(var , fin);
}

bool AreaLight::Collide(Vector3 ray_O , Vector3 ray_V) {
	ray_V = ray_V.GetUnitVector();
	Vector3 N = (Dx * Dy).GetUnitVector();
	double d = N.Dot(ray_V);
	if (fabs(d) < EPS) return false;
	double l = (N * O.Dot(N) - ray_O).Dot(N) / d;
	if (l < EPS) return false;

	Vector3 C = (ray_O + ray_V * l) - O;
	if (fabs(Dx.Dot(C)) > Dx.Dot(Dx)) return false;
	if (fabs(Dy.Dot(C)) > Dy.Dot(Dy)) return false;

	intersect_distant = l;
	return true;
}

double AreaLight::CalnShade(Vector3 C , Primitive* primitive_head , int shade_quality) {
	int shade = 0;
	
	for (int i = -2 ; i < 2 ; i++) {
		for (int j = -2 ; j < 2 ; j++) {
			for (int k = 0 ; k < shade_quality ; k++) {
				Vector3 V = O - C + Dx * ((ran() + i) / 2) + Dy * ((ran() + j) / 2);
				double distant = V.Module();

				for (Primitive* now = primitive_head ; now != NULL ; now = now->GetNext())
					if (now->Collide(C , V) && (now->intersect.distant < distant)) {
						shade++;
						break;
					}
			}
		}
	}
	
	return 1 - (double) shade / (16.0 * shade_quality);
}

//Camera
//标准镜头的长（或宽）与感光点到镜头距离的比值
const double STD_LENS_WIDTH = 0.88;
const double STD_LENS_HEIGHT = 0.88;
//标准照片的象素长（或宽）
const int STD_IMAGE_WIDTH = 420;
const int STD_IMAGE_HEIGHT = 420;
//~*16=标准情况下，蒙特卡罗计算阴影时的计算次数
const int STD_SHADE_QUALITY = 4;
//~*16=标准情况下，蒙特卡罗计算镜面漫反射时的计算次数
const int STD_DREFL_QUALITY = 4;

Camera::Camera() {
	O = Vector3(0 , 0 , 0);
	N = Vector3(0 , 1 , 0);
	lens_W = STD_LENS_WIDTH;
	lens_H = STD_LENS_HEIGHT;
	W = STD_IMAGE_WIDTH;
	H = STD_IMAGE_HEIGHT;
	shade_quality = STD_SHADE_QUALITY;
	drefl_quality = STD_DREFL_QUALITY;
	data = NULL;
}

Camera::~Camera() {
	if (data == NULL) {
		for (int i = 0 ; i < H ; i++)
			delete[] data[i];
		delete[] data;
	}
}

void Camera::Initialize() {
	N = N.GetUnitVector();
	Dx = N * Vector3(0 , 0 , 1);
	Dy = Dx * N;
	Dx = Dx * lens_W / 2;
	Dy = Dy * lens_H / 2;

	data = new Color*[H];
	for (int i = 0 ; i < H ; i++)
		data[i] = new Color[W];
}

//得到象素(i,j)对应的射出光线
Vector3 Camera::Emit(double i , double j) {
	return N + Dy * (2 * i / H - 1) + Dx * (2 * j / W - 1);
}

void Camera::Input(std::string var , std::stringstream& fin) {
	if (var == "O=") O.Input(fin);
	if (var == "N=") N.Input(fin);
	if (var == "lens_W=") fin >> lens_W;
	if (var == "lens_H=") fin >> lens_H;
	if (var == "image_W=") fin >> W;
	if (var == "image_H=") fin >> H;
	if (var == "shade_quality=") fin >> shade_quality;
	if (var == "drefl_quality=") fin >> drefl_quality;
}

void Camera::Output(Bmp* bmp) {
	bmp->Initialize(H , W);

	for (int i = 0 ; i < H ; i++)
		for (int j = 0 ; j < W ; j++)
			bmp->SetColor(i , j , data[i][j]);
}


//Scene
Scene::Scene() {
	primitive_head = NULL;
	light_head = NULL;
	background_color = Color();
	camera = new Camera;
}

Scene::~Scene() {
	while (primitive_head != NULL) {
		Primitive* next_head = primitive_head->GetNext();
		if (primitive_head->GetMaterial()->texture != NULL) {
			delete primitive_head->GetMaterial()->texture;
		}
		delete primitive_head;
		primitive_head = next_head;
	}

	while (light_head != NULL) {
		Light* next_head = light_head->GetNext();
		delete light_head;
		light_head = next_head;
	}

	delete camera;
}

//ReadIn
void Scene::CreateScene(std::string file) {
	srand(1995 - 05 - 12);
	std::ifstream fin(file.c_str());

	std::string obj;
	while (fin >> obj) {
		Primitive* new_primitive = NULL;
		Light* new_light = NULL;

		if (obj == "primitive") {
			std::string type; fin >> type;
			if (type == "sphere") new_primitive = new Sphere;
			if (type == "plane") new_primitive = new Plane;
			if (new_primitive != NULL) {
				new_primitive->SetNext(primitive_head);
				primitive_head = new_primitive;
			}
		}
		else if (obj == "light") {
			std::string type; fin >> type;
			if (type == "point") new_light = new PointLight;
			if (type == "area") new_light = new AreaLight;
			if (new_light != NULL) {
				new_light->SetNext(light_head);
				light_head = new_light;
			}
		}
		else if (obj != "background" && obj != "camera")
			continue;

		fin.ignore(1024 , '\n');
		
		std::string order;
		while (getline(fin , order , '\n')) {
			std::stringstream fin2(order);
			std::string var; fin2 >> var;
			if (var == "end") break;

			if (obj == "background" && var == "color=") background_color.Input(fin2);
			if (obj == "primitive" && new_primitive != NULL) new_primitive->Input(var , fin2);
			if (obj == "light" && new_light != NULL) new_light->Input(var , fin2);
			if (obj == "camera") camera->Input(var , fin2);
		}
	}

	camera->Initialize();
}

//计算光线(ray_O,ray_V)碰到的第一个物体，碰撞信息储存在物体的intersect变量中（=NULL时表示碰不到任何物体）
Primitive* Scene::FindNearestPrimitive(Vector3 ray_O , Vector3 ray_V) {
	Primitive* ret = NULL;

	for (Primitive* now = primitive_head ; now != NULL ; now = now->GetNext())
		if (now->Collide(ray_O , ray_V) && (ret == NULL || now->intersect.distant < ret->intersect.distant)) ret = now;

	return ret;
}

//计算光线(ray_O,ray_V)碰到的第一个光源，碰撞前光走过的距离储存在物体的intersect_distant变量中（=NULL时表示碰不到任何光源）
Light* Scene::FindNearestLight(Vector3 ray_O , Vector3 ray_V) {
	Light* ret = NULL;

	for (Light* now = light_head ; now != NULL ; now = now->GetNext())
		if (now->Collide(ray_O , ray_V) && (ret == NULL || now->intersect_distant < ret->intersect_distant)) ret = now;

	return ret;
}


//Raytracer
const double SPEC_POWER = 20;
const int MAX_DREFL_DEP = 2;
const int MAX_RAYTRACING_DEP = 10;
const int HASH_FAC = 7;
const int HASH_MOD = 10000007;

//计算漫反射色光，其中pri是刚碰撞到的物品的指针，碰撞信息储存在intersect变量中，hash用于记录光线经过的物品（hash值不同的相邻象素要使用超级取样去锯齿）
Color Raytracer::CalnDiffusion(Primitive* pri , int* hash) {
	Color color = pri->GetMaterial()->color;
	if (pri->GetMaterial()->texture != NULL) color = color * pri->GetTexture();
	
	Color ret = color * scene.GetBackgroundColor() * pri->GetMaterial()->diff;

	for (Light* light = scene.GetLightHead() ; light != NULL ; light = light->GetNext()) {
		double shade = light->CalnShade(pri->intersect.C , scene.GetPrimitiveHead() , scene.GetCamera()->GetShadeQuality());
		if (shade < EPS) continue;
		
		Vector3 R = (light->GetO() - pri->intersect.C).GetUnitVector();
		double dot = R.Dot(pri->intersect.N);
		if (dot > EPS) {
			if (hash != NULL && light->IsPointLight()) *hash = (*hash + light->GetSample()) & HASH_MOD;

			if (pri->GetMaterial()->diff > EPS) {
				double diff = pri->GetMaterial()->diff * dot * shade;
				ret += color * light->GetColor() * diff;
			}
			if (pri->GetMaterial()->spec > EPS) {
				double spec = pri->GetMaterial()->spec * pow(dot , SPEC_POWER) * shade;
				ret += color * light->GetColor() * spec;
			}
		}
	}

	return ret;
}

//计算反射后得到的色光，其中ray_V是入射光线的方向，dep是迭代层数
Color Raytracer::CalnReflection(Primitive* pri , Vector3 ray_V , int dep , int* hash) {
	ray_V = ray_V.Reflect(pri->intersect.N);

	if (pri->GetMaterial()->drefl < EPS || dep > MAX_DREFL_DEP)
		return RayTracing(pri->intersect.C , ray_V , dep + 1 , hash) * pri->GetMaterial()->color * pri->GetMaterial()->refl;

	Vector3 Dx = ray_V * Vector3(1 , 0 , 0);
	if (Dx.IsZeroVector()) Dx = Vector3(1 , 0 , 0);
	Vector3 Dy = ray_V * Dx;
	Dx = Dx.GetUnitVector() * pri->GetMaterial()->drefl;
	Dy = Dy.GetUnitVector() * pri->GetMaterial()->drefl;

	Color ret;
	for (int k = 0 ; k < 16 * scene.GetCamera()->GetDreflQuality() ; k++) {
		double x , y;
		do {
			x = ran() * 2 - 1;
			y = ran() * 2 - 1;
		} while (x * x + y * y > 1);
		x *= pri->GetMaterial()->drefl;
		y *= pri->GetMaterial()->drefl;

		ret += RayTracing(pri->intersect.C , ray_V + Dx * x + Dy * y , dep + MAX_DREFL_DEP , NULL);
	}

	ret = ret * pri->GetMaterial()->color * pri->GetMaterial()->refl / (16 * scene.GetCamera()->GetDreflQuality());
	return ret;
}

//计算折射后得到的色光
Color Raytracer::CalnRefraction(Primitive* pri , Vector3 ray_V , int dep , int* hash) {
	double n = pri->GetMaterial()->rindex;
	if (pri->intersect.front) n = 1 / n;
	
	ray_V = ray_V.Refract(pri->intersect.N , n);
	
	Color rcol = RayTracing(pri->intersect.C , ray_V , dep + 1 , hash);
	if (pri->intersect.front) return rcol * pri->GetMaterial()->refr;
	Color absor = pri->GetMaterial()->absor * -pri->intersect.distant;
	Color trans = Color(exp(absor.r) , exp(absor.g) , exp(absor.b));
	return rcol * trans * pri->GetMaterial()->refr;
}

//光线追踪的主要迭代函数，衔接不同层次的光线
Color Raytracer::RayTracing(Vector3 ray_O , Vector3 ray_V , int dep , int* hash) {
	if (dep > MAX_RAYTRACING_DEP) return Color();

	Color ret;
	Primitive* nearest_primitive = scene.FindNearestPrimitive(ray_O , ray_V);
	Light* nearest_light = scene.FindNearestLight(ray_O , ray_V);

	if (nearest_light != NULL && (nearest_primitive == NULL || nearest_light->intersect_distant < nearest_primitive->intersect.distant)) {
		if (hash != NULL) *hash = (*hash + nearest_light->GetSample()) % HASH_MOD;
		ret += nearest_light->GetColor();
	}
	
	if (nearest_primitive != NULL) {
		if (hash != NULL) *hash = (*hash + nearest_primitive->GetSample()) % HASH_MOD;
		Primitive* primitive = nearest_primitive->PrimitiveCopy();
		if (primitive->GetMaterial()->diff > EPS || primitive->GetMaterial()->spec > EPS) ret += CalnDiffusion(primitive , hash);
		if (primitive->GetMaterial()->refl > EPS) ret += CalnReflection(primitive , ray_V , dep , hash);
		if (primitive->GetMaterial()->refr > EPS) ret += CalnRefraction(primitive , ray_V , dep , hash);
		delete primitive;
	}

	if (hash != NULL) *hash = (*hash * HASH_FAC) % HASH_MOD;
	ret.Valid();
	return ret;
}
//开始渲染，根据输入txt文件中的数据，渲染出图像导出在输出bmp文件中
void Raytracer::Run() {
	Camera* camera = scene.GetCamera();
	scene.CreateScene(input);

	Vector3 ray_O = camera->GetO();
	int H = camera->GetH() , W = camera->GetW();
	int** sample = new int*[H];
	for (int i = 0 ; i < H ; i++) {
		sample[i] = new int[W];
		for (int j = 0 ; j < W ; j++) {
			sample[i][j] = 0;
		}
	}

	for (int i = 0 ; i < H ; ++i)
		for (int j = 0 ; j < W ; j++) {
			Vector3 ray_V = camera->Emit(i , j);
			Color color = RayTracing(ray_O , ray_V , 4 , &sample[i][j]);
			camera->SetColor(i , j , color);
		}

		for (int i = 0 ; i < H ; ++i)
			for (int j = 0 ; j < W ; j++) {
				if ((i == 0 || sample[i][j] == sample[i - 1][j]) && (i == H - 1 || sample[i][j] == sample[i + 1][j]) &&
					 (j == 0 || sample[i][j] == sample[i][j - 1]) && (j == W - 1 || sample[i][j] == sample[i][j + 1])) continue;

				Color color;
				for (int r = -1 ; r <= 1 ; r++)
					for (int c = -1 ; c <= 1 ; c++) {
						Vector3 ray_V = camera->Emit(i + (double) r / 3 , j + (double) c / 3);
						color += RayTracing(ray_O , ray_V , 4 , NULL) / 9;//antiAlias去锯齿
					}
				camera->SetColor(i , j , color);
			}
	
	for (int i = 0 ; i < H ; i++)
		delete[] sample[i];
	delete[] sample;

	Bmp* bmp = new Bmp(H , W);
	camera->Output(bmp);
	bmp->Output(output);
	delete bmp;
}

