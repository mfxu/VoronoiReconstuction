
#include "Optimize_patical.h"
//#define FIX_SAMPLE

Optimize_patical::Optimize_patical() {
	queryPt = annAllocPt(3);
	nnIdx = new ANNidx[Neighbor_NUM];						// allocate near neigh indices
	dists = new ANNdist[Neighbor_NUM];						// allocate near neighbor dists
	default_mls_radius = 1.4;
	mlsTime = 0.0;
}
Optimize_patical::~Optimize_patical() {
	delete[] nnIdx;							// clean things up
	delete[] dists;
	delete points_kdtree;
	annClose();
}
void Optimize_patical::init_solver(  double offset_ratio ,  double sigma_coefficient, int sam_num, char* file_name)
{
	mesh = new CRichModel(file_name);
	mesh->LoadModel();
	samples_num = sam_num;
	points_num = mesh->m_Verts.size();
	points_kdtree = GetKDTree(3, points_num, mesh->m_Verts);
	analye_point_cloud();
	off_radius = diagonal*offset_ratio;
	produce_samples(off_radius);
	dis_between_neighbor = compute_sigma(samples_set);
	init_variables(samples_set);
	
	avg_dis.clear();
	avg_dis.resize(samples_num);
	neighbor.clear();
	neighbor.resize(samples_num);
	sigma = dis_between_neighbor*sigma_coefficient;
	//sigma = 0.25* sqrtf(mesh_area / samples_num);
}
void Optimize_patical::set_sigma(double sigma_coefficient)
{
	dis_between_neighbor = compute_sigma(samples_set);
	sigma = dis_between_neighbor*sigma_coefficient;
	//sigma = 0.25* sqrtf(mesh_area / samples_num);
	cout << "init offset_radius: " << off_radius << "\tdis_for_points: " << (dis_between_neighbor) << "\tsigma:" << sigma << endl;
}
void Optimize_patical::analye_point_cloud()
{
	bounding_box[0].x = bounding_box[1].x = mesh->m_Verts[0].x;
	bounding_box[0].y = bounding_box[1].y = mesh->m_Verts[0].y;
	bounding_box[0].z = bounding_box[1].z = mesh->m_Verts[0].z;

	for (int pi = 0; pi < points_num; pi++)
	{
		if (mesh->m_Verts[pi].x < bounding_box[0].x)
			bounding_box[0].x = mesh->m_Verts[pi].x;
		if (mesh->m_Verts[pi].x > bounding_box[1].x)
			bounding_box[1].x = mesh->m_Verts[pi].x;

		if (mesh->m_Verts[pi].y < bounding_box[0].y)
			bounding_box[0].y = mesh->m_Verts[pi].y;
		if (mesh->m_Verts[pi].y > bounding_box[1].y)
			bounding_box[1].y = mesh->m_Verts[pi].y;

		if (mesh->m_Verts[pi].z < bounding_box[0].z)
			bounding_box[0].z = mesh->m_Verts[pi].z;
		if (mesh->m_Verts[pi].z > bounding_box[1].z)
			bounding_box[1].z = mesh->m_Verts[pi].z;
	}
	diagonal = (bounding_box[1] - bounding_box[0]).Len();
	axis_length = bounding_box[1] - bounding_box[0];
	max_axis = axis_length.x > axis_length.y ? axis_length.x : axis_length.y;
	max_axis = (max_axis > axis_length.z ? max_axis : axis_length.z);
	min_axis = axis_length.x < axis_length.y ? axis_length.x : axis_length.y;// 
	min_axis = (min_axis < axis_length.z ? min_axis : axis_length.z);
	point_cloud_center = (bounding_box[1] + bounding_box[0]) / 2.0;
	mesh_area = 0;
	int e1, e2;
	CPoint3D edge1, edge2;
	double angle = 0;
	if (mesh->m_Faces.size() != 0)
	{
		for (int fi = 0; fi < mesh->m_Faces.size(); fi++)
		{
			e1 = mesh->m_Faces[fi].verts[0];
			e2 = mesh->m_Faces[fi].verts[1];
			edge1 = mesh->m_Verts[e2] - mesh->m_Verts[e1];
			e2 = mesh->m_Faces[fi].verts[2];
			edge2 = mesh->m_Verts[e2] - mesh->m_Verts[e1];
			angle = acosf(Dot(edge1, edge2) / (edge1.Len()*edge2.Len()));
			mesh_area += 0.5* sin(angle)*edge1.Len()*edge2.Len();
		}
	}
	else
	{
		mesh_area = 2.0*( axis_length.x*axis_length.y + axis_length.y*axis_length.z + axis_length.x*axis_length.z);
	}
}
void Optimize_patical::produce_samples(double distortion)
{
#ifdef FIX_SAMPLE
	CPoint3D onep;
	float x, y, z;
	int pindex;
	char buf[200];
	FILE* init_samples;
	fopen_s(&init_samples, "result/example/init_patical.obj","r");
	
	fgets(buf, 200, init_samples);
	while( !feof(init_samples) )
	{
		sscanf_s(buf, "v %f %f %f", &x, &y, &z);
		onep.x = x;	onep.y = y;		onep.z = z;
		samples_set.push_back(onep);
		fgets(buf, 200, init_samples);
	}
	fclose(init_samples);
#else
	CPoint3D onep;
	srand(time(0));
	CPoint3D randirect;
	int pindex;
	ofstream spher_ini("result/init_patical.obj");
	for (int si = 0; si < samples_num; si++)
	{
		pindex = abs((int)((rand() + 0.0) / RAND_MAX * mesh->m_Verts.size() - 1));
		do {
			randirect.x = (rand() + 0.0) / RAND_MAX - 0.5;
			randirect.y = (rand() + 0.0) / RAND_MAX - 0.5;
			randirect.z = (rand() + 0.0) / RAND_MAX - 0.5;
		} while (randirect.Len() == 0.0);
		randirect.Normalize();

		onep = mesh->m_Verts[pindex] + randirect * distortion;
		spher_ini << "vn " << onep.x << ' ' << onep.y << ' ' << onep.z << endl;
		spher_ini << "v " << onep.x << ' ' << onep.y << ' ' << onep.z << endl;
		samples_set.push_back(onep);
	}
	spher_ini.close();
#endif
	

}
double Optimize_patical::compute_sigma(vector<CPoint3D> samples)
{
	double temp_dis;
	ANNkd_tree* skdtree = GetKDTree(3, samples_num, samples);
	ANNpoint			squeryPt = annAllocPt(3);
	ANNidxArray			snnIdx = new ANNidx[7];			// near neighbor indices
	ANNdistArray		sdists = new ANNdist[7]; 			// near neighbor distances

	sample_neighbor temp_neihbor; 
	double points_dis= 0;
	double avg_sigma;
	for (int i = 0; i < samples_num; i++)
	{
		temp_dis = 0;
		squeryPt[0] = samples[i].x;	squeryPt[1] = samples[i].y;	squeryPt[2] = samples[i].z;
		skdtree->annkSearch(squeryPt, 7, snnIdx, sdists);
		for (int ni = 0; ni < 6; ni++)
		{
			temp_neihbor.nei_inx[ni] = snnIdx[ni + 1];
			temp_neihbor.nei_dis[ni] = sdists[ni + 1];
			temp_dis += sqrtf(sdists[ni + 1]);
		}
		avg_sigma = temp_dis / 6.0;
		points_dis += avg_sigma;
	}
	points_dis /= (samples_num + 0.0);

	delete[] snnIdx;							// clean things up
	delete[] sdists;
	delete skdtree;
	annClose();
	return points_dis;
}
int    Optimize_patical::compute_dis_neighbors_from_file()
{
	int number = 3;
	ifstream usam("result/optimized_result.obj");
	char buf[256];
	vector<CPoint3D> uniformed_samples;
	float x, y, z;
	CPoint3D ons;
	int spn = 0;
	while (!usam.eof())
	{
		usam.getline(buf, 200);
		int num =sscanf_s(buf, "v %f %f %f", &x, &y, &z);
		if ( num ==3 )
		{
			ons.x = x; ons.y = y; ons.z = z;
			uniformed_samples.push_back(ons);
		}
	}
	spn = uniformed_samples.size();
	usam.close();
	//create kdTree
	ANNpointArray		dataPts;				// data uniformed_samples
	dataPts = annAllocPts(spn, 3);			// allocate data uniformed_samples
	for (int i = 0; i < spn; i++)
	{
		dataPts[i][0] = uniformed_samples[i].x;
		dataPts[i][1] = uniformed_samples[i].y;
		dataPts[i][2] = uniformed_samples[i].z;
	}
	ANNkd_tree* kdtree = new ANNkd_tree(dataPts, spn, 3);
	ANNpoint			queryPt;
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		nndists;					// near neighbor distances
	queryPt = annAllocPt(3);
	nnIdx = new ANNidx[number + 1];
	nndists = new ANNdist[number + 1];
	neigbors.clear();
	neigbors.resize(spn);
	for (int si = 0; si < spn; si++)
	{
		queryPt[0] = uniformed_samples[si].x;
		queryPt[1] = uniformed_samples[si].y;
		queryPt[2] = uniformed_samples[si].z;
		kdtree->annkSearch(queryPt, number + 1, nnIdx, nndists, 0.0);
		for (int ni = 0; ni < number; ni++)
		{
			neigbors[si].nei_dis[ni] = nndists[ni + 1];
			neigbors[si].nei_inx[ni] = nnIdx[ni + 1];
		}
	}
	// cluster the particles on the offset surface.
	vector<bool> checked;	checked.clear();	checked.resize(spn, false);
	int numOfCluster = 0;
	int nextID;
	while (true)
	{
		bool found = false;
		int topID = 0;
		for (topID; topID < spn; topID++)
		{
			if (!checked[topID])
			{
				found = true;
				break;
			}
		}
		if (!found)
			break;

		vector<int> cluster;
		vector<int> cluster_set;
		cluster.push_back(topID);
		cluster_set.push_back(topID);
		checked[topID] = true;
		while (!cluster.empty())
		{
			nextID = cluster.back();
			cluster.pop_back();
			int nid;
			for (int ni = 0; ni < 6; ++ni)
			{
				nid = neigbors[nextID].nei_inx[ni];
				if (checked[nid])
					continue;
				cluster.push_back(nid);
				cluster_set.push_back(nid);
				checked[nid] = true;
			}
		}
		samples_cluster.push_back(cluster_set);
		++numOfCluster;
	}

	largest_cluster_index = 0;
	nextID = 0;
	cout << "components inforamtion:\n";
	for (int ni = 0; ni < samples_cluster.size(); ni++)
	{
		cout << "\tpoints in cluster " << ni << ":\t is \t" << samples_cluster[ni].size() << endl;
		if (samples_cluster[ni].size() > nextID)
		{
			nextID = samples_cluster[ni].size();
			largest_cluster_index = ni;
		}
	}
	ofstream outo("result/outer.obj");
	for (int vi = 0; vi < samples_cluster[largest_cluster_index].size(); vi++)
	{
		nextID = samples_cluster[largest_cluster_index][vi];
		outo << "v " << uniformed_samples[nextID].x << ' ' << uniformed_samples[nextID].y << ' ' << uniformed_samples[nextID].z << endl;
	}
	outo.close();
	ofstream outi("result/inner.obj");
	for (int many = 0; many < samples_cluster.size(); ++many)
	{ 
		if (many != largest_cluster_index&& samples_cluster[many].size()>500)
		{
			outi << "#cluster " << many << std::endl;
			for (int vi = 0; vi < samples_cluster[many].size(); vi++)
			{
				nextID = samples_cluster[many][vi];
				outi << "v " << uniformed_samples[nextID].x << ' ' << uniformed_samples[nextID].y << ' ' << uniformed_samples[nextID].z << endl;
			}
		}
	}
	outi.close();
	

	delete[] nnIdx;							// clean things up
	delete[] nndists;
	delete kdtree;
	annClose();
	return numOfCluster;
}

void Optimize_patical::init_variables(vector<CPoint3D> init_samples)
{
	// Fill in the cloud data
	tree = pcl::search::KdTree<pcl::PointXYZ>::Ptr(new pcl::search::KdTree<pcl::PointXYZ>);
	// init samples to pcd format
	pcd_samples = pcl::PointCloud<pcl::PointXYZ>::Ptr( new pcl::PointCloud<pcl::PointXYZ>());
	pcd_samples->width = init_samples.size();
	pcd_samples->height = 1;
	pcd_samples->is_dense = false;
	pcd_samples->points.resize(pcd_samples->width * pcd_samples->height);
	//init point surface to pcd format
	pcd_orig_surface = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>());
	pcd_orig_surface->width = mesh->m_Verts.size();
	pcd_orig_surface->height = 1;
	pcd_orig_surface->is_dense = false;
	pcd_orig_surface->points.resize(pcd_orig_surface->width * pcd_orig_surface->height);
	for (int pi = 0; pi < pcd_orig_surface->width; pi++)
	{
		pcd_orig_surface->points[pi].x = mesh->m_Verts[pi].x;
		pcd_orig_surface->points[pi].y = mesh->m_Verts[pi].y;
		pcd_orig_surface->points[pi].z = mesh->m_Verts[pi].z;
	}
	//init pcl mls process
	mls.setComputeNormals(true);
	mls.setSearchMethod(tree);
	mls.setPolynomialOrder(3);
	//mls.setDistinctCloud(pcd_samples);
	
	mls.setUpsamplingMethod(mls.DISTINCT_CLOUD);
	mls.setInputCloud(pcd_orig_surface);
	m_variables.clear();
	m_variables.resize(samples_num * 3, 0);
	mls_pedal.clear();
	mls_pedal.resize(samples_num);
	previous_samples.clear();
	previous_samples.resize(samples_num);
	for (int si = 0; si < samples_num; ++si)
	{
		m_variables[si * 3]     = init_samples[si].x;
		m_variables[si * 3+1] = init_samples[si].y;
		m_variables[si * 3+2] = init_samples[si].z;
		previous_samples[si]  = init_samples[si];
	}
}
lbfgsfloatval_t Optimize_patical::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
	
	char findex[10];
	_itoa_s(iter_time, findex, 9);
	//char before[256] = "result/particle/optimized_";
	//strcat_s(before, findex);
	//strcat_s(before, ".obj");
	//Output_updated_samples(before, mesh);
	//cout << "Offset surface..." << endl;

	Offset_samples_from_implicit_surface();
	Offset_samples_from_implicit_surface();

	//char after[256] = "result/particle/offset_";
	//strcat_s(after, findex);
	//strcat_s(after, ".obj");
	//Output_updated_samples(after, mesh);
	
	cout << iter_time++ << endl;
	double function_value = 0;// Objection_function();
	CPoint3D center;
	double dis;
	/******************************************
	*******  patical system   gradient *******
	******************************************/
	clock_t ts = clock();
	double evalue = 0;
	for (int sk = 0; sk < samples_num; sk++)
	{
		g[sk * 3 + 0] = 0;
		g[sk * 3 + 1] = 0;
		g[sk * 3 + 2] = 0;
		CPoint3D spoint = CPoint3D(m_variables[sk * 3], m_variables[sk * 3 + 1], m_variables[sk * 3 + 2]);
		center.x = center.y = center.z = 0.0;
		for (int ssi = 0; ssi< samples_num; ++ssi)
		{
			CPoint3D tp = CPoint3D(m_variables[ssi * 3], m_variables[ssi * 3 + 1], m_variables[ssi * 3 + 2]);
			dis = (spoint - tp).Len();
			if ( ssi != sk && (dis < sigma * 50))
			{
				evalue = exp(-(dis*dis) / (4 * sigma*sigma));
				center += evalue*(tp - spoint) / (2 * sigma*sigma);
				if (ssi > sk )
					function_value += evalue;
			}
		}
		g[ sk * 3 + 0] += center.x;
		g[ sk * 3 + 1] += center.y;
		g[ sk * 3 + 2] += center.z;
	}
	clock_t te = clock();
	cout << "function value is: "<< function_value<<"; time in computing gradients: "<< te-ts<<"ms.\n";
	return function_value;
}
double Optimize_patical::Objection_function()
{
	double f = 0;
	double inner_sum = 0;
	double dis;
	for (int sn = 0; sn < samples_num; ++sn)
	{
		CPoint3D sp(m_variables[sn * 3], m_variables[sn * 3 + 1], m_variables[sn * 3 + 2]);
		for (int kn = sn + 1; kn < samples_num; ++kn)
		{
			CPoint3D kp(m_variables[kn * 3], m_variables[kn * 3 + 1], m_variables[kn * 3 + 2]);
			dis = (sp - kp).Len();
			if (dis<sigma *50)
				inner_sum += exp( -(dis*dis)/(4*sigma*sigma)) ;
		}
	}
	f += LAMBDA3*inner_sum;
	return f;
}
void Optimize_patical::update_partical_location()
{

	auto updated_samples = GetVariables();
	samples_set.clear();
	samples_set.resize(samples_num);
	ofstream out_samples("result/optimized_result.obj");
	for (int si = 0; si < samples_num; si++)
	{
		samples_set[si].x = updated_samples[si * 3 + 0];
		samples_set[si].y = updated_samples[si * 3 + 1];
		samples_set[si].z = updated_samples[si * 3 + 2];
		out_samples << "v " << updated_samples[si * 3 + 0] << ' ' << updated_samples[si * 3 + 1] << ' ' << updated_samples[si * 3 + 2] << endl;
	}
	out_samples.close();
}
int    Optimize_patical::compute_dis_neighbors()
{
	update_partical_location();
	//create kdTree
	ANNpointArray		dataPts;				// data samples_set
	dataPts = annAllocPts(samples_num, 3);			// allocate data samples_set
	for (int i = 0; i < samples_num; i++)
	{
		dataPts[i][0] = samples_set[i].x;
		dataPts[i][1] = samples_set[i].y;
		dataPts[i][2] = samples_set[i].z;
	}
	ANNkd_tree* kdtree = new ANNkd_tree(dataPts, samples_num, 3);
	ANNpoint			queryPt;
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		nndists;					// near neighbor distances
	queryPt = annAllocPt(3);
	nnIdx = new ANNidx[Neighbor_NUM + 1];
	nndists = new ANNdist[Neighbor_NUM + 1];
	neigbors.clear();
	neigbors.resize(samples_num);
	for (int si = 0; si < samples_num; si++)
	{
		queryPt[0] = samples_set[si].x;
		queryPt[1] = samples_set[si].y;
		queryPt[2] = samples_set[si].z;
		kdtree->annkSearch(queryPt, Neighbor_NUM + 1, nnIdx, nndists, 0.0);
		for (int ni = 0; ni < Neighbor_NUM; ni++)
		{
			neigbors[si].nei_dis[ni] = nndists[ni + 1];
			neigbors[si].nei_inx[ni] = nnIdx[ni + 1];
		}
	}
	// cluster the particles on the offset surface.
	vector<bool> checked;	checked.clear();	checked.resize(samples_num, false);
	int numOfCluster = 0;
	int nextID;
	while (true)
	{
		bool found = false;
		int topID = 0;
		for (topID; topID < samples_num; topID++)
		{
			if (!checked[topID])
			{
				found = true;
				break;
			}
		}
		if (!found)
			break;

		vector<int> cluster;
		vector<int> cluster_set;
		cluster.push_back(topID);
		cluster_set.push_back(topID);
		checked[topID] = true;
		while (!cluster.empty())
		{
			nextID = cluster.back();
			cluster.pop_back();
			int nid;
			for (int ni = 0; ni < 6; ++ni)
			{
				nid = neigbors[nextID].nei_inx[ni];
				if (checked[nid])
					continue;
				cluster.push_back(nid);
				cluster_set.push_back(nid);
				checked[nid] = true;
			}
		}
		samples_cluster.push_back(cluster_set);
		++numOfCluster;
	}

	largest_cluster_index = 0;
	nextID = 0;
	cout << "components inforamtion:\n";
	for (int ni = 0; ni < samples_cluster.size(); ni++)
	{
		cout << "\tpoints in cluster " << ni << ":\t is \t" << samples_cluster[ni].size() << endl;
		if (samples_cluster[ni].size() > nextID)
		{
			nextID = samples_cluster[ni].size();
			largest_cluster_index = ni;
		}
	}
	ofstream outo("result/largest_outer_offset.obj");
	for (int vi = 0; vi < samples_cluster[largest_cluster_index].size(); vi++)
	{
		nextID = samples_cluster[largest_cluster_index][vi];
		outo << "v " << samples_set[nextID].x << ' ' << samples_set[nextID].y << ' ' << samples_set[nextID].z << endl;
	}
	outo.close();

	ofstream outi("result/inner.obj");
	for (int vi = 0; vi < samples_cluster.size(); vi++)
	{
		if (vi != largest_cluster_index && samples_cluster[vi].size() > 500)
		{
			for (int vj = 0; vj < samples_cluster[vi].size(); vj++)
			{
				nextID = samples_cluster[vi][vj];
				outi << "v " << samples_set[nextID].x << ' ' << samples_set[nextID].y << ' ' << samples_set[nextID].z << endl;
			}
		}
	}
	outi.close();

	delete[] nnIdx;							// clean things up
	delete[] nndists;
	delete kdtree;
	annClose();
	return numOfCluster;
}

void Optimize_patical::Output_updated_samples(char* new_samples_set, CRichModel* reference_mesh)
{
	ofstream out(new_samples_set, ios::out | ios::trunc);
	if (out.fail())
		throw "fail to read file";
	for (int i = 0; i < samples_num; i++)
	{
		out << "v  " << m_variables[i * 3] << ' ';
		out << m_variables[i * 3 + 1] << ' ';
		out << m_variables[i * 3 + 2] << " 1.0 0.0 0.0\n";
	}

	/*for (int i = 0; i < mesh->m_Verts.size(); i++)
	{
		out << "v  " << mesh->m_Verts[i].x << ' ';
		out << mesh->m_Verts[i].y << ' ';
		out << mesh->m_Verts[i].z << " 0.0 1.0 0.0\n";
	}

	for (int i = 0; i < mesh->m_Faces.size(); i++)
	{
		out << "f  " << mesh->m_Faces[i].verts[0] + 1 + samples_num << ' ';
		out << mesh->m_Faces[i].verts[1] + 1 + samples_num << ' ';
		out << mesh->m_Faces[i].verts[2] + 1 + samples_num << '\n';
	}*/
	
	out.close();
}
void Optimize_patical::Offset_samples_from_implicit_surface()
{
#ifdef NEAREST
	clock_t ts = clock();
	//mls_projection();

	double nearest_dis_to_surface;
	CPoint3D pedal, direction, offsetP;
	ofstream outc("result/pedal.obj");
	ofstream larf("result/one_offset.obj");
	for (int si = 0; si < samples_num; si++)
	{
		CPoint3D sp(m_variables[si * 3], m_variables[si * 3 + 1], m_variables[si * 3 + 2]);
		queryPt[0] = m_variables[si * 3];
		queryPt[1] = m_variables[si * 3+1];
		queryPt[2] = m_variables[si * 3+2];
		points_kdtree->annkSearch(queryPt, 1, nnIdx, dists);
		mls_pedal[si] = mesh->m_Verts[nnIdx[0]];
		pedal = mls_pedal[si];
		direction = sp - pedal;
		if (direction.Len() == 0.0)
			cout << "special procession\n";
		//if (nearest_dis_to_surface > off_radius*1.5)//使用粒子点到模型的距离
		if ((sp - pedal).Len() > off_radius*1.5)//
		{
			CPoint3D ranP = random_distorb(off_radius*0.1);
			previous_samples[si] = previous_samples[si] + ranP;
			m_variables[si * 3 + 0] = previous_samples[si].x;
			m_variables[si * 3 + 1] = previous_samples[si].y;
			m_variables[si * 3 + 2] = previous_samples[si].z;
		}
		else {
			direction.Normalize();
			offsetP = pedal + direction*off_radius;
			previous_samples[si] = offsetP;
			m_variables[si * 3 + 0] = offsetP.x;
			m_variables[si * 3 + 1] = offsetP.y;
			m_variables[si * 3 + 2] = offsetP.z;
		}
		previous_samples[si].x = m_variables[si * 3 + 0];
		previous_samples[si].y = m_variables[si * 3 + 1];
		previous_samples[si].z = m_variables[si * 3 + 2];

		samples_set[si].x = m_variables[si * 3 + 0];
		samples_set[si].y = m_variables[si * 3 + 1];
		samples_set[si].z = m_variables[si * 3 + 2];
		larf << "v " << samples_set[si].x << ' ' << samples_set[si].y << ' ' << samples_set[si].z << '\n';
		outc << "v " << mls_pedal[si].x << ' ' << mls_pedal[si].y << ' ' << mls_pedal[si].z << '\n';
	}
	outc.close();
	larf.close();
	clock_t te = clock();
	cout << "offsetting in :" << te - ts << "ms" << endl;
#else 
	clock_t ts = clock();
	clock_t mlst = clock();
	mls_projection();
	clock_t mlse = clock();
	mlsTime += (mlse - mlst);
	cout << "once mls time:" << mlse - mlst << "ms"<<endl;
	double nearest_dis_to_surface;
	CPoint3D pedal, direction, offsetP;
	//ofstream outc("result/pedal.obj");
	//ofstream larf("result/one_offset.obj");
	for (int si = 0; si < samples_num; si++)
	{
		CPoint3D sp(m_variables[si * 3], m_variables[si * 3 + 1], m_variables[si * 3 + 2]);

		nearest_dis_to_surface = GetclosetNeighorPointsDistance(points_kdtree, sp, 1);
		pedal = mls_pedal[si];
		//outc << "v " << pedal.x << ' ' << pedal.y << ' ' << pedal.z << '\n';
		direction = sp - pedal;
		if (direction.Len() == 0.0)
			cout << "special procession\n";
		//if (nearest_dis_to_surface > off_radius*1.5)//使用粒子点到模型的距离
		//if ((sp - pedal).Len() > off_radius*1.5)//使用粒子点到垂足的距离
		if ((sp - pedal).Len() > off_radius*1.5 || nearest_dis_to_surface >off_radius*1.80)//使用粒子点到垂足的距离
		{
			CPoint3D ranP = random_distorb(off_radius*0.1);
			previous_samples[si] = previous_samples[si] + ranP;
			double lengh_tag = GetclosetNeighorPointsDistance(points_kdtree, previous_samples[si], 1);
			if (lengh_tag > off_radius*1.80)
			{
				direction = previous_samples[si] - mesh->m_Verts[nearest_index];
				direction.Normalize();
				previous_samples[si] = mesh->m_Verts[nearest_index] + direction * off_radius;
			}

			m_variables[si * 3 + 0] = previous_samples[si].x;
			m_variables[si * 3 + 1] = previous_samples[si].y;
			m_variables[si * 3 + 2] = previous_samples[si].z;
		}
		else {
			direction.Normalize();
			offsetP = pedal + direction*off_radius;
			previous_samples[si] = offsetP;
			m_variables[si * 3 + 0] = offsetP.x;
			m_variables[si * 3 + 1] = offsetP.y;
			m_variables[si * 3 + 2] = offsetP.z;
		}
		previous_samples[si].x = m_variables[si * 3 + 0];
		previous_samples[si].y = m_variables[si * 3 + 1];
		previous_samples[si].z = m_variables[si * 3 + 2];

		samples_set[si].x = m_variables[si * 3 + 0];
		samples_set[si].y = m_variables[si * 3 + 1];
		samples_set[si].z = m_variables[si * 3 + 2];
		//larf << "v " << samples_set[si].x << ' ' << samples_set[si].y << ' ' << samples_set[si].z << '\n';
	}
	//outc.close();
	//larf.close();
	clock_t te = clock();
	cout << "once optimization in :" << te - ts << "ms" << endl;
	

#endif

}
void Optimize_patical::mls_projection()
{
	double pedal_to_surface, sample_to_surface;
	CPoint3D sp;
	//ofstream outs("result/one_update_particles.obj");
	for (int pi = 0; pi < pcd_samples->width; pi++)
	{
		pcd_samples->points[pi].x = m_variables[pi * 3];
		pcd_samples->points[pi].y = m_variables[pi * 3 + 1];
		pcd_samples->points[pi].z = m_variables[pi * 3 + 2];
		//outs << "v " << pcd_samples->points[pi].x << ' ' << pcd_samples->points[pi].y << ' ' << pcd_samples->points[pi].z << '\n';
	}
	//outs.close();

	mls.setDistinctCloud(pcd_samples);
	default_mls_radius = 1.2;
	//default_mls_radius = 1.4;
	do {
		default_mls_radius = default_mls_radius + 0.1;
		mls.setSearchRadius(off_radius*default_mls_radius);
		mls.process(mls_projected_points);
		//cout << "pedal num:"<<mls_projected_points.width << endl;
	} while (mls_projected_points.width < samples_num);

	//ofstream outc1("result/mls_pedal.obj");
	//ofstream mls_outc("result/mls_problem.obj");
	mls_pedal.clear();
	mls_pedal.resize(samples_num);
	for (int pi = 0; pi < samples_num; pi++)
	{
		//queryPt[0] = mls_projected_points.points[pi].x;
		//queryPt[1] = mls_projected_points.points[pi].y;
		//queryPt[2] = mls_projected_points.points[pi].z;
		queryPt[0] = m_variables[3*pi];
		queryPt[1] = m_variables[3 * pi+1];
		queryPt[2] = m_variables[3 * pi+2];
		//outc1 << "v " << mls_projected_points.points[pi].x << ' ' << mls_projected_points.points[pi].y << ' ' << mls_projected_points.points[pi].z << endl;

		points_kdtree->annkSearch(queryPt, 2, nnIdx, dists, 0.0);  //particle to the surface distance .
		CPoint3D pedal2sample = CPoint3D(queryPt[0], queryPt[1], queryPt[2]) - CPoint3D(mls_projected_points.points[pi].x, mls_projected_points.points[pi].y, mls_projected_points.points[pi].z);
		pedal2sample.Normalize();
		sample_to_surface = sqrtf(dists[0]) ;
		queryPt[0] = mls_projected_points.points[pi].x;
		queryPt[1] = mls_projected_points.points[pi].y;
		queryPt[2] = mls_projected_points.points[pi].z;
		points_kdtree->annkSearch(queryPt, 2, nnIdx, dists, 0.0);  //particle to the surface distance .
		CPoint3D nearest2sample = CPoint3D(queryPt[0], queryPt[1], queryPt[2]) - CPoint3D(mesh->m_Verts[nnIdx[0]].x, mesh->m_Verts[nnIdx[0]].y, mesh->m_Verts[nnIdx[0]].z);
		nearest2sample.Normalize();
		double angle = pedal2sample.x*nearest2sample.x + pedal2sample.y*nearest2sample.y + pedal2sample.z*nearest2sample.z;
		pedal_to_surface = sqrtf(dists[0]);
		// the pedal should on the surface and the sample should not be far away from the surface , and the line between pedal and the corresponding points are 
		if (sample_to_surface  >  off_radius*1.3 && pedal_to_surface > off_radius*0.3 )
		{
//			cout << pi << '\t';
			mls_pedal[pi].x = mesh->m_Verts[nnIdx[0]].x;
			mls_pedal[pi].y = mesh->m_Verts[nnIdx[0]].y;
			mls_pedal[pi].z = mesh->m_Verts[nnIdx[0]].z;
		}
		else
		{
			mls_pedal[pi].x = mls_projected_points.points[pi].x;
			mls_pedal[pi].y = mls_projected_points.points[pi].y;
			mls_pedal[pi].z = mls_projected_points.points[pi].z;
		}
	}
//	cout << endl;
//	outc1.close();
}
ANNkd_tree* Optimize_patical::GetKDTree(int dim, int sn, vector <CPoint3D> &points)const
{
	//create kdTree
	ANNpointArray		dataPts;				// data points
	dataPts = annAllocPts(sn, dim);			// allocate data points
	for (int i = 0; i < sn; i++)
	{
		dataPts[i][0] = points[i].x;
		dataPts[i][1] = points[i].y;
		dataPts[i][2] = points[i].z;
	}
	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, sn, dim);
	return kdTree;
}
double Optimize_patical::GetclosetNeighorPointsDistance(ANNkd_tree* kdTree, CPoint3D queryPoint, int num)
{
	queryPt[0] = queryPoint.x;
	queryPt[1] = queryPoint.y;
	queryPt[2] = queryPoint.z;
	//cerr << "kdtree's point : " << kdTree->nPoints() << endl;
	points_kdtree->annkSearch(queryPt, num, nnIdx, dists, 0.0);
	nearest_index = nnIdx[0];
	return (sqrtf(dists[0]));
}
