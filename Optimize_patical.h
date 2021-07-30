#ifndef __OPTIMIZE_PATICAL_H__
#define __OPTIMIZE_PATICAL_H__

#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
//#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/impl/kdtree.hpp>
#include <pcl/surface/mls.h>
#include <ANN/ANN.h>
#include "Point3D.h"
#include "ReadMOdel\RichModel.h"
//#define NEAREST


#define SQR(a) ((a)*(a))
#pragma once
#include "LBFGS\Non_Linear_Optimization.h"
struct sample_neighbor
{
	double nei_dis[6];
	int		nei_inx[6];
};
class Optimize_patical :
	public Non_Linear_Optimization
{
public:
	Optimize_patical();
	~Optimize_patical();
	vector<CPoint3D> mls_pedal;
	void init_solver(double offset_ratio, double sigma_coefficient,  int sam_num, char* file_name) ;
	void set_sigma(double sigma_coefficient);
	lbfgsfloatval_t Optimize_patical::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	int compute_dis_neighbors();
	int compute_dis_neighbors_from_file();
	double mlsTime;

private:	
	int iter_time;
	double sigma;
	int samples_num = 0;
	int points_num = 0;
	int kernel_num;
	vector<CPoint3D> rbf_center;
	vector<double> center_weight;
	vector<double> samples_value;
	vector<double> points_value;
	vector<double> off_samples_value;
	vector<CPoint3D> samples_set;
	vector<CPoint3D> previous_samples;
	vector<double>	avg_dis;
	
	vector<sample_neighbor>	neighbor;
	ANNkd_tree* points_kdtree;
	ANNpoint			queryPt;
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances

	pcl::PointCloud<pcl::PointXYZ>::Ptr pcd_samples;
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcd_orig_surface;
	pcl::PointCloud<pcl::PointNormal> mls_projected_points;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree;
	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
	double default_mls_radius = 1.5,dis_between_neighbor;
	double off_radius , diagonal, mesh_area, max_axis, min_axis;
	int nearest_index;
	CPoint3D bounding_box[2], point_cloud_center,axis_length;
	CRichModel *mesh;
	vector<sample_neighbor> neigbors;
	vector < vector<int> > samples_cluster;
	int largest_cluster_index;

private:
	double Objection_function();
	void update_partical_location();
	void Output_updated_samples(char* new_samples_set, CRichModel* reference_mesh);
	void analye_point_cloud();
	void produce_samples(double distortion);
	void init_variables(vector<CPoint3D> init_samples);
	ANNkd_tree* GetKDTree(int dim, int pointNum, vector <CPoint3D> &points)const;
	double GetclosetNeighorPointsDistance(ANNkd_tree* kdTree, CPoint3D queryPoint, int num);
	void Optimize_patical::Offset_samples_from_implicit_surface();
	double compute_sigma(vector<CPoint3D> samples);
	void Optimize_patical::mls_projection();


	inline CPoint3D random_distorb(double dis)
	{
		CPoint3D udis;
		udis.x = (rand() + 0.0 )/ RAND_MAX-05; 
		udis.y = (rand() + 0.0) / RAND_MAX-0.5;
		udis.z = (rand() + 0.0) / RAND_MAX-0.5;
		double dadtag = (rand() + 0.0) / RAND_MAX;
		int tag = 0;
		if (dadtag >= 0.5)
			tag = -1;
		else
			tag = 1;
		udis.Normalize();
		udis = udis*dis*tag;
		return udis;
	}
};


#endif // !OPTIMIZE_PATICAL_H__
