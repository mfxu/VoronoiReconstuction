#include "Optimize_patical.h"
//#include "LBFGS_Solver.h"
#include <cmath>
using namespace std;

//		#define KNOT
//		#define BLADE
//		#define  PELVIS
//		#define VENUS
//		#define FISH
//		#define SHEET
//		#define HEPTOROID
//		#define SCISSOR
//		#define pulley
//		#define TWOCHAIR
//		#define  SMOOTH_FEATURE
//		#define FLOWER
//		#define OILPUMP
//		#define INSTANCE_EAP
//		#define PLANT
//		#define LEAF1
//		#define oiliverHand
		#define elk
#ifdef SPHERE
char mesh_reference[50] = "result/sphere.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_1460.txt";
#endif // SPHERE



#ifdef KITTEN
char mesh_reference[50] = "data/kitten.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_kitten5773.txt";
#endif // EIGHT

#ifdef TUOYUAN
char mesh_reference[50] = "data/sphere_tuoyuan_dense.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_tuoyuan.txt";
#endif // EIGHT

#ifdef TORUS
char mesh_reference[50] = "data/torus.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_torus.txt";
#endif // TORUS

#ifdef SHEET
char mesh_reference[50] = "data/sheet_point.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_tuoyuan.txt";
#endif // TORUS

#ifdef  KNOT
char mesh_reference[50] = "data/Blade.obj";
#endif //  KNOTdd 

#ifdef SCISSOR
char mesh_reference[50] = "result/scissors_mesh.obj";
#endif // SCISSOR

#ifdef  HALF
char mesh_reference[50] = "test_data/venus_half.obj";
#endif 

#ifdef BLADE
char mesh_reference[50] = "test_data/blade.obj";
char acturally_rbf_parameter[50] = "data/rbf_para_tuoyuan.txt";
#endif // SHADE

#ifdef FISH
char mesh_reference[50] = "data/fish_model.obj";
#endif // fish

#ifdef PELVIS
char mesh_reference[50] = "test_data/pelvis.obj";
#endif // PELVIS

#ifdef pulley
char mesh_reference[50] = "test_data/pulley_surface_noise.obj";
#endif // PELVIS

#ifdef HEPTOROID
char mesh_reference[50] = "test_data/heptoroid72323.obj";
#endif // HEPTOROID

#ifdef VENUS
char mesh_reference[50] = "test_data/venus-half-noisy.obj";
#endif // VENUS

#ifdef SMOOTH_FEATURE
char mesh_reference[50] = "result/smooth_feature.obj";
#endif // VENUS

#ifdef FLOWER
//char mesh_reference[50] = "flower/flower_start.obj";
//char mesh_reference[50] = "flower/flower82.obj";
char mesh_reference[50] = "flower/flower23.obj";
#endif // FLOWER
#ifdef OILPUMP
char mesh_reference[100] = "test_data/oilpump_denoise_liujian_unconsistent.obj";
#endif // OILPUMP
#ifdef INSTANCE_EAP
//char mesh_reference[100] = "test_data/outline_input.obj";
char mesh_reference[100] = "test_data/sphere15k_noise.obj";
#endif // INSTANCE_EAP
#ifdef PLANT
char mesh_reference[100] = "test_data/laurierentire_noise_surface.obj";
#endif
#ifdef TWOCHAIR
char mesh_reference[100] = "test_data/one_chair_surface.obj";
#endif

#ifdef LEAF1
char mesh_reference[100] = "test_data/leaf/leaf1_point.obj";
#endif

#ifdef HOLE
char mesh_reference[100] = "test_data/leaf/smoothEightPoint .obj";
#endif

#ifdef oiliverHand
char mesh_reference[100] = "test_data/leaf/hand-olivier.obj";
#endif

#ifdef elk
char mesh_reference[100] = "test_data/leaf/casting.obj";
#endif

int main(int argc, char *argv[])
{
	Optimize_patical  func;
	/*******************************************
	first parameter:			offset_ratio
	second							sigma_ratio
	third								particles num
	fourth							the input point cloud.
	*******************************************/
	//func.init_solver(0.012,  0.1,  30000, mesh_reference); //sheet
	//func.init_solver(0.015,  0.25, 80000, mesh_reference);//fish
	//func.compute_dis_neighbors_from_file();
	double sig_para =0.5;
	// func.init_solver(0.007, sig_para, 50000, mesh_reference);  //for pulley  
	//func.init_solver(0.01, sig_para, 30000, mesh_reference);	//flower 30
	//func.init_solver(0.01, sig_para, 22000, mesh_reference);	//flower 3
	//func.compute_dis_neighbors_from_file();
	int succed = 0;
	func.init_solver(0.02, sig_para, 22000, mesh_reference);	//example22k;		  sig_para = 0.3;

	//for (sig_para = 0.3; sig_para<5; sig_para = sig_para +0.2)
	{
		for (int i = 0; i<5; ++i)
		{
			func.set_sigma(sig_para);
			succed = func.run();
			//if (succed == -1000)
			//	break;
		}
	}
		
	func.compute_dis_neighbors();
	func.SaveEnergySequence("r-", "enery_sequnce.m");
	func.SaveNormSequence("r-", "norm_sequnce.m");
	//func.compute_dis_neighbors_from_file();
	cout << "MLS time  is " << func.mlsTime << endl<<endl;
	system("pause");
	return 0;
}

