//$Id: surface_integral.h,v 1.9 2013/08/20 11:49:05 zjcao Exp $
#ifndef SURFACE_INTEGRAL_H
#define SURFACE_INTEGRAL_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream>
#include <math.h>
#endif

#include "cgh.h"
#include "ShellPatch.h"
#include "NullShellPatch.h"
#include "NullShellPatch2.h"
#include "var.h"
#include "monitor.h"
#include "bh_interp.h" // BH_OPT: 插值几何预计算缓存

class surface_integral
{

private:
	int Symmetry, factor;
	int N_theta, N_phi; // Number of points in Theta & Phi directions
	double dphi, dcostheta;
	double *arcostheta, *wtcostheta;
	int n_tot; // size of arrays

	double *nx_g, *ny_g, *nz_g; // global list of unit normals
	int myrank, cpusize;

	// BH_OPT: 各探测器球面的插值几何缓存（prepare_benchmark 在计时区外填充）
	int n_cache;
	double *cache_rex;
	bh_detector_cache *cache;

	// BH_OPT2: 球谐(Wigner-d)分解查表——这些量只依赖 (theta 下标 i, phi 下标 j,
	// l,m 序号)，与最外层 n 循环无关，却被放在最内层重算。预先算好存表。
	double *tab_thetap; // [N_theta][2][NN]：2 对应 costheta 取 +arcostheta[i] / -
	double *tab_cos;	// [lpsy][2*maxl+1][N_phi]
	double *tab_sin;	// 同上
	int tab_NN, tab_maxl, tab_spinw;
	void build_decomp_tables(int spinw, int maxl, int NN);

public:
	surface_integral(int iSymmetry);
	~surface_integral();
	void begin_wave_benchmark(int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
							  int field_generation);
	void end_wave_benchmark();
	void prepare_benchmark(int decn, double maxrex, double drex, int lev,
						   cgh *GH); // BH_OPT: 预计算各探测器球面插值几何（计时区外）
	// BH_OPT5: 把「含点块」的最新场数据广播给所有 rank。必须在计时区外调用
	// （紧挨着 Surface_Only_Benchmark 的 MPI_Barrier 之前），每轮一次即可，
	// 3 个探测器共用同一份场数据。
	void sync_block_copies(var *Rpsi4, var *Ipsi4);
	static void report_prof(); // BH_PROF: 临时分解基准探针（提交前移除）

	void surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor); // NN is the length of RP and IP
									  // this routine can only deal with the symmetry of Psi4
	void surf_Wave(double rex, int lev, ShellPatch *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, NullShellPatch *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, NullShellPatch2 *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, ShellPatch *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor); // NN is the length of RP and IP
	void surf_Wave(double rex, int lev, cgh *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor,
				   void (*funcs)(double &, double &, double &,
								 double &, double &, double &, double &, double &, double &, double &,
								 double &, double &, double &, double &, double &, double &,
								 double &, double &)); // NN is the length of RP and IP
	void surf_Wave(double rex, int lev, ShellPatch *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor,
				   void (*funcs)(double &, double &, double &,
								 double &, double &, double &, double &, double &, double &, double &,
								 double &, double &, double &, double &, double &, double &,
								 double &, double &)); // NN is the length of RP and IP
	void surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
					   double *Rout, monitor *Monitor);
	void surf_MassPAng(double rex, int lev, ShellPatch *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
					   double *Rout, monitor *Monitor);
	void surf_Wave(double rex, cgh *GH, ShellPatch *SH,
				   var *chi, var *trK,
				   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
				   var *chix, var *chiy, var *chiz,
				   var *trKx, var *trKy, var *trKz,
				   var *Axxx, var *Axxy, var *Axxz,
				   var *Axyx, var *Axyy, var *Axyz,
				   var *Axzx, var *Axzy, var *Axzz,
				   var *Ayyx, var *Ayyy, var *Ayyz,
				   var *Ayzx, var *Ayzy, var *Ayzz,
				   var *Azzx, var *Azzy, var *Azzz,
				   var *Gamxxx, var *Gamxxy, var *Gamxxz, var *Gamxyy, var *Gamxyz, var *Gamxzz,
				   var *Gamyxx, var *Gamyxy, var *Gamyxz, var *Gamyyy, var *Gamyyz, var *Gamyzz,
				   var *Gamzxx, var *Gamzxy, var *Gamzxz, var *Gamzyy, var *Gamzyz, var *Gamzzz,
				   var *Rxx, var *Rxy, var *Rxz, var *Ryy, var *Ryz, var *Rzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	bool SR_Interp_Points(MyList<var> *VarList, cgh *GH, ShellPatch *SH,
						  int NN, double **XX, double *Shellf);

	void surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
					   double *Rout, monitor *Monitor, MPI_Comm Comm_here);
	void surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor, MPI_Comm Comm_here);
};

extern void bh_prof_mpatch_report(double out[2]); // BH_PROF: 临时分解基准探针（提交前移除）
#endif /* SURFACE_INTEGRAL_H */
