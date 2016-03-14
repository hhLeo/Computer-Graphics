#include "func.h"

int main()
{
	//kitten();
	
	double t1 = clock();
	Raytracer* raytracer = new Raytracer;
	raytracer->SetInput("input.txt");
	raytracer->SetOutput("output.bmp");
	raytracer->Run();
	double t2 = clock();
	cout << "time : " << (t2 - t1) / 1000.0 << endl;
	//while(1);
	return 0;
}