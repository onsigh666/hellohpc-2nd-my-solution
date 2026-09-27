
#ifndef PATCH_H
#define PATCH_H

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "var.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width
#include "bh_interp.h" // BH_OPT: 插值几何预计算缓存

class Patch
{

public:
   int lev;
   int shape[dim];
   double bbox[2 * dim]; // this bbox includes buffer points
   MyList<Block> *blb, *ble;
   int lli[dim], uui[dim]; // denote the buffer points on each boundary

public:
   Patch() {};
   Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry);

   ~Patch();

   void checkPatch(bool buflog);
   void checkPatch(bool buflog, const int out_rank);
   void checkBlock();
   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX,
                      double *Shellf, int Symmetry,
                      const bh_point_geom *geom = 0,    // BH_OPT: geom 非空时走融合内核
                      int n_lo = -1, int n_hi = -1,     // BH_OPT3: 本 rank 只需 [n_lo,n_hi]，-1=全部
                      const int *mine = 0, int n_mine = -1, // BH_OPT4: 只扫这些点，NULL/负=全部
                      const bh_block_copy *bc = 0);     // BH_OPT5: 含点块的数据副本
   void prepare_point_geom(int NN, double **XX, int Symmetry,
                           bh_point_geom *geom,               // BH_OPT: 计时区外预计算每点插值几何
                           bh_block_copy **pbc, int *pn_bc);  // BH_OPT5: 输出「含点块」副本数组
   bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                         double *Shellf, int Symmetry);
   double getdX(int dir);

   void Find_Maximum(MyList<var> *VarList, double *XX,
                     double *Shellf);

   bool Find_Point(double *XX);

   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX,
                      double *Shellf, int Symmetry, MPI_Comm Comm_here);
   bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                         double *Shellf, int Symmetry, MPI_Comm Comm_here);
   void Find_Maximum(MyList<var> *VarList, double *XX,
                     double *Shellf, MPI_Comm Comm_here);
};

#endif /* PATCH_H */
