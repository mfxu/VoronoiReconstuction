#pragma once
#include "Point3D.h"
#include <vector>

#include <stdlib.h> 

using namespace std;
class iso_surface
{
public:
	iso_surface(int vo)
	{
		resolution = vo;
	}
	void Region_devide();
	void Tag_vertices(char* rbf);

	~iso_surface();
	CPoint3D*** d_vertices;
	int ***	d_vertices_tag;
private:
	int resolution;
	CPoint3D bounding_box[2];
};

