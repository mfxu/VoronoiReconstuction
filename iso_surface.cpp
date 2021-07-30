
#include "iso_surface.h"




iso_surface::~iso_surface()
{
}

void iso_surface::Region_devide() 
{
	double x_axis = bounding_box[1].x - bounding_box[0].x;
	double y_axis = bounding_box[1].y - bounding_box[0].y;
	double z_axis = bounding_box[1].z - bounding_box[0].z;

	double step = x_axis / resolution ;
	int x_v = floor(x_axis / step), y_v= floor(y_axis/step), z_v=floor(z_axis/step);

	**d_vertices = new CPoint3D[x_v];
	**d_vertices_tag=new int[x_v];
	for (int x = 0; x < x_v; ++x)
	{
		*d_vertices[x] = new CPoint3D[y_v];
		*d_vertices_tag[x] = new int[y_v];
		for (int y = 0; y < y_v; ++y)
		{
			d_vertices[x][y] = new CPoint3D[z_v];
			d_vertices_tag[x][y] = new int[z_v];
			for (int z = 0; z < z_v; ++z)
			{
				d_vertices[x][y][z] = CPoint3D(bounding_box[x].x*x, bounding_box[x].y*y, bounding_box[x].z*z);
				d_vertices_tag[x][y][z] = 0;
			}
		}	
	}
		
}

void iso_surface::Tag_vertices(char* rbf)
{


}