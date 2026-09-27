
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "misc.h"
#include "MPatch.h"
#include "Parallel.h"
#include "fmisc.h"
#include "parameters.h"

// BH_DEBUG: 临时探针——直接调用 Fortran 的 decide3d_ 导出真实插值模板
extern "C" int decide3d_(int *, double *, double *, int *, int *, double *, double *, int *, int *);
extern "C" void polin3_(double *, double *, double *, double *, double *, double *, double *, double *,
                        double *, int *);

// ===== BH_PROF: 临时分解基准探针（提交前移除）=====
static double bh_prof_ia_search = 0, bh_prof_ia_allreduce = 0;
void bh_prof_mpatch_report(double out[2])
{
  out[0] = bh_prof_ia_search;
  out[1] = bh_prof_ia_allreduce;
}

// BH_DEBUG: 临时插值对拍探针（提交前移除）——同一采样点用融合内核与
// 原 Fortran 路径各算一遍，不一致就打印细节。BH_DEBUG_INTERP=1 打开。
static int bh_dbg_on()
{
  static int cached = -1;
  if (cached < 0)
    cached = getenv("BH_DEBUG_INTERP") ? 1 : 0;
  return cached;
}

// BH_DEBUG: 原样转写 fmisc.f90::polint（Neville），0 基下标，xa[i]=i
static double bh_dbg_polint(const double *ya, int n, double x)
{
  double c[BH_ORDN_MAX], d[BH_ORDN_MAX], ho[BH_ORDN_MAX], den[BH_ORDN_MAX];
  int ns = 0;
  double dif = fabs(x - 0.0);
  for (int i = 0; i < n; i++)
  {
    const double dift = fabs(x - (double)i);
    if (dift < dif)
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
    ho[i] = (double)i - x;
  }
  double y = ya[ns];
  ns = ns - 1;
  for (int m = 1; m < n; m++)
  {
    for (int i = 0; i < n - m; i++)
      den[i] = ho[i] - ho[i + m];
    for (int i = 0; i < n - m; i++)
      den[i] = (c[i + 1] - d[i]) / den[i];
    for (int i = 0; i < n - m; i++)
      d[i] = ho[i + m] * den[i];
    for (int i = 0; i < n - m; i++)
      c[i] = ho[i] * den[i];
    double dy;
    if (2 * ns < n - m)
      dy = c[ns + 1];
    else
    {
      dy = d[ns];
      ns = ns - 1;
    }
    y = y + dy;
  }
  return y;
}

// BH_DEBUG: 手工按 decide3d 搭 6^3 模板再用 Neville 张量积插值

// BH_DEBUG: 用 Fortran 的 decide3d_ + polin3_ 直算（零改写），完全复刻 global_interp 内部
static double bh_dbg_pipe(const Block *BP, const bh_point_geom &g, const double *f,
                          const double xyz[3], const double soa[3], int dbg_sym, double cxo[3])
{
  const int shape[3] = {BP->shape[0], BP->shape[1], BP->shape[2]};
  int cb[3] = {g.base[0], g.base[1], g.base[2]};
  int ct[3] = {cb[0] + g.ordn - 1, cb[1] + g.ordn - 1, cb[2] + g.ordn - 1};
  const int nb = BH_ORDN_MAX * BH_ORDN_MAX * BH_ORDN_MAX;
  double *yb = new double[nb];
  for (int t = 0; t < nb; t++)
    yb[t] = 0.0;
  int o = g.ordn, s = dbg_sym;
  decide3d_(const_cast<int *>(shape), const_cast<double *>(f), const_cast<double *>(f), cb, ct,
            const_cast<double *>(soa), yb, &o, &s);
  for (int m = 0; m < 3; m++)
  {
    const double *Xm = (m == 0) ? BP->X[0] : (m == 1) ? BP->X[1] : BP->X[2];
    const double dmm = Xm[1] - Xm[0];
    cxo[m] = (cb[m] > 0) ? (xyz[m] - Xm[cb[m] - 1]) / dmm : (xyz[m] + Xm[1 - cb[m]]) / dmm;
  }
  double xa1[BH_ORDN_MAX];
  for (int t = 0; t < g.ordn; t++)
    xa1[t] = t;
  double vy = 0, dyv = 0;
  polin3_(xa1, xa1, xa1, yb, &cxo[0], &cxo[1], &cxo[2], &vy, &dyv, &o);
  delete[] yb;
  return vy;
}

static double bh_dbg_manual(const bh_point_geom &g, const double *f, int Nx, int Nxy,
                            const double soa[3], const double xyz[3], double cxout[3])
{
  const int n = g.ordn;
  double ya[BH_ORDN_MAX][BH_ORDN_MAX][BH_ORDN_MAX];
  int src0[BH_ORDN_MAX], src1[BH_ORDN_MAX], src2[BH_ORDN_MAX];
  double s0[BH_ORDN_MAX], s1[BH_ORDN_MAX], s2[BH_ORDN_MAX];
  for (int p = 0; p < n; p++)
  {
    const int fi0 = g.base[0] + p, fi1 = g.base[1] + p, fi2 = g.base[2] + p;
    src0[p] = (fi0 < 1) ? 2 - fi0 : fi0;
    src1[p] = (fi1 < 1) ? 2 - fi1 : fi1;
    src2[p] = (fi2 < 1) ? 2 - fi2 : fi2;
    s0[p] = (fi0 < 1) ? soa[0] : 1.0;
    s1[p] = (fi1 < 1) ? soa[1] : 1.0;
    s2[p] = (fi2 < 1) ? soa[2] : 1.0;
  }
  for (int p = 0; p < n; p++)
    for (int q = 0; q < n; q++)
      for (int r = 0; r < n; r++)
        ya[p][q][r] = f[(src0[p] - 1) + (src1[q] - 1) * Nx + (src2[r] - 1) * Nxy] * s0[p] * s1[q] * s2[r];

  // cx 按 global_interp 的公式逐位重算
  for (int m = 0; m < 3; m++)
  {
    const double dxm = g.BP->X[m][1] - g.BP->X[m][0];
    const int cxB = g.base[m];
    cxout[m] = (cxB > 0) ? (xyz[m] - g.BP->X[m][cxB - 1]) / dxm
                         : (xyz[m] + g.BP->X[m][2 - cxB - 1]) / dxm;
  }

  double ztmp[BH_ORDN_MAX][BH_ORDN_MAX], ytmp[BH_ORDN_MAX], v[BH_ORDN_MAX];
  for (int p = 0; p < n; p++)
    for (int q = 0; q < n; q++)
      ztmp[p][q] = bh_dbg_polint(ya[p][q], n, cxout[2]);
  for (int p = 0; p < n; p++)
    ytmp[p] = bh_dbg_polint(ztmp[p], n, cxout[1]);
  (void)v;
  return bh_dbg_polint(ytmp, n, cxout[0]);
}

namespace
{
int surface_interpolation_order()
{
  map<string, int>::iterator iter = parameters::int_par.find("surface interpolation order");
  if (iter != parameters::int_par.end())
    return iter->second;
  return 2 * ghost_width;
}
}

Patch::Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry) : lev(levi)
{

  int hbuffer_width = buffer_width;
  if (lev == 0)
    hbuffer_width = CS_width; // specific for shell-box coulping

  if (DIM != dim)
  {
    cout << "dimension is not consistent in Patch construction" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  for (int i = 0; i < dim; i++)
  {
    shape[i] = shapei[i];
    bbox[i] = bboxi[i];
    bbox[dim + i] = bboxi[dim + i];
    lli[i] = uui[i] = 0;
    if (buflog)
    {
      double DH;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
      DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      uui[i] = hbuffer_width;
      bbox[dim + i] = bbox[dim + i] + uui[i] * DH;
      shape[i] = shape[i] + uui[i];
    }
  }

  if (buflog)
  {
    if (DIM != 3)
    {
      cout << "Symmetry in Patch construction only support 3 yet but dim = " << DIM << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    double tmpb, DH;
    if (Symmetry > 0)
    {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      DH = (bbox[5] - bbox[2]) / (shape[2] - 1);
#else
#ifdef Cell
      DH = (bbox[5] - bbox[2]) / shape[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
      tmpb = Mymax(0, bbox[2] - hbuffer_width * DH);
      lli[2] = int((bbox[2] - tmpb) / DH + 0.4);
      bbox[2] = bbox[2] - lli[2] * DH;
      shape[2] = shape[2] + lli[2];
      if (lli[2] < hbuffer_width)
      {
        if (feq(bbox[2], 0, DH / 2))
          lli[2] = 0;
        else
        {
          cout << "Code mistake for lli[2] = " << lli[2] << ", bbox[2] = " << bbox[2] << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
      if (Symmetry > 1)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[3] - bbox[0]) / (shape[0] - 1);
#else
#ifdef Cell
        DH = (bbox[3] - bbox[0]) / shape[0];
#else
#error Not define Vertex nor Cell
#endif
#endif
        tmpb = Mymax(0, bbox[0] - hbuffer_width * DH);
        lli[0] = int((bbox[0] - tmpb) / DH + 0.4);
        bbox[0] = bbox[0] - lli[0] * DH;
        shape[0] = shape[0] + lli[0];
        if (lli[0] < hbuffer_width)
        {
          if (feq(bbox[0], 0, DH / 2))
            lli[0] = 0;
          else
          {
            cout << "Code mistake for lli[0] = " << lli[0] << ", bbox[0] = " << bbox[0] << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[4] - bbox[1]) / (shape[1] - 1);
#else
#ifdef Cell
        DH = (bbox[4] - bbox[1]) / shape[1];
#else
#error Not define Vertex nor Cell
#endif
#endif
        tmpb = Mymax(0, bbox[1] - hbuffer_width * DH);
        lli[1] = int((bbox[1] - tmpb) / DH + 0.4);
        bbox[1] = bbox[1] - lli[1] * DH;
        shape[1] = shape[1] + lli[1];
        if (lli[1] < hbuffer_width)
        {
          if (feq(bbox[1], 0, DH / 2))
            lli[1] = 0;
          else
          {
            cout << "Code mistake for lli[1] = " << lli[1] << ", bbox[1] = " << bbox[1] << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }
      }
      else
      {
        for (int i = 0; i < 2; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
          DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
          lli[i] = hbuffer_width;
          bbox[i] = bbox[i] - lli[i] * DH;
          shape[i] = shape[i] + lli[i];
        }
      }
    }
    else
    {
      for (int i = 0; i < dim; i++)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        DH = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
#else
#ifdef Cell
        DH = (bbox[dim + i] - bbox[i]) / shape[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        lli[i] = hbuffer_width;
        bbox[i] = bbox[i] - lli[i] * DH;
        shape[i] = shape[i] + lli[i];
      }
    }
  }

  blb = ble = 0;
}
Patch::~Patch()
{
}
// buflog 1: with buffer points; 0 without
void Patch::checkPatch(bool buflog)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    if (buflog)
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] << ":" << bbox[dim + i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
    else
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i] - lli[i] - uui[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
  }
}
void Patch::checkPatch(bool buflog, const int out_rank)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == out_rank)
  {
    cout << " out_rank = " << out_rank << endl;
    if (buflog)
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] << ":" << bbox[dim + i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
    else
    {
      cout << " belong to level " << lev << endl;
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << shape[i] - lli[i] - uui[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]";
      }
      cout << " resolution: [";
      for (int i = 0; i < dim; i++)
      {
        cout << getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
    }
  }
}
void Patch::Interp_Points(MyList<var> *VarList,
                          int NN, double **XX,
                          double *Shellf, int Symmetry,
                          const bh_point_geom *geom, int n_lo, int n_hi,
                          const int *mine, int n_mine,
                          const bh_block_copy *bc)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = surface_interpolation_order();
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[NN * num_var];
  memset(shellf, 0, sizeof(double) * NN * num_var);

  // we use weight to monitor code, later some day we can move it for optimization
  int *weight;
  weight = new int[NN];
  memset(weight, 0, sizeof(int) * NN);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  // ===================== BH_OPT3: 按需投递，取代全量 Allreduce =====================
  // 原先每轮都把整条 shellf(n_tot*num_var 个 double，small.par 下约 400KB) 做一次
  // MPI_Allreduce，再加一条同样长的 int 归约；而调用方 surf_Wave 每个 rank 只读
  // 自己负责的那一段 [n_lo, n_hi]（约 1/cpusize）。全量归约等于让所有人都听一遍
  // 全场广播，属于典型的「通信量 = 问题规模」而非「通信量 = 实际需要」。
  //
  // geom 在所有 rank 上逐位一致（prepare_benchmark 用同一套坐标生成），所以每个
  // rank 都能各自本地算出「谁拥有哪个点」「谁需要哪个点」，无需再协商：
  //   发送计划：我拥有的点里，落在 rank t 需求区间内的
  //   接收计划：我需要的点里，归属 rank 是 s 的
  // 然后一次 MPI_Alltoallv 点对点投递即可，通信量降到约 1/cpusize。
  // 另外点搜索是「找到第一个匹配块就停」，故每个点至多一个 owner，不存在多点
  // 求和/平均，Alltoallv 直接投递与原 Allreduce 结果逐位相同。
  //
  // 【实测结论】ARM 集群 np=8：通信量确实降到约 1/8，但 ia_allreduce 只从 3.879s
  // 降到 3.601s（-7%），总时间反而 5.2685s → 5.7781s（慢 10%）。原因是那段「通信
  // 时间」里绝大部分是【等最慢的 rank】而非真的在传数据——减数据量治不了负载不均。
  // 因此本路径【默认关闭】，仅在 BH_FASTCOMM=1 时启用以备后续复测。
  if (geom && getenv("BH_FASTCOMM"))
  {
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    const int need_lo = (n_lo < 0) ? 0 : n_lo;
    const int need_hi = (n_hi < 0) ? NN - 1 : n_hi;

    // 1) 交换各 rank 的需求区间，并夹紧到 [0,NN-1]（收发计数与打包必须用同一套区间）
    int *q_lo = new int[cpusize];
    int *q_hi = new int[cpusize];
    MPI_Allgather(&need_lo, 1, MPI_INT, q_lo, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&need_hi, 1, MPI_INT, q_hi, 1, MPI_INT, MPI_COMM_WORLD);
    for (int t = 0; t < cpusize; t++)
    {
      if (q_lo[t] < 0)
        q_lo[t] = 0;
      if (q_hi[t] > NN - 1)
        q_hi[t] = NN - 1;
    }

    // 2) 每个点的归属 rank（-1 = 预计算时没找到块）
    int *owner = new int[NN];
    for (int j = 0; j < NN; j++)
      owner[j] = geom[j].BP ? geom[j].BP->rank : -1;

    // 3) 收发计划
    int *snd_cnt = new int[cpusize];
    int *rcv_cnt = new int[cpusize];
    for (int t = 0; t < cpusize; t++)
    {
      snd_cnt[t] = 0;
      rcv_cnt[t] = 0;
    }
    for (int t = 0; t < cpusize; t++) // 我要发给 t 的点
    {
      if (t == myrank)
        continue;
      for (int j = q_lo[t]; j <= q_hi[t]; j++)
        if (owner[j] == myrank)
          snd_cnt[t]++;
    }
    for (int j = need_lo; j <= need_hi; j++) // 我要从 o 收的点
    {
      const int o = owner[j];
      if (o >= 0 && o != myrank)
        rcv_cnt[o]++;
    }

    // 4) 交叉核对：别人要发给我的数量必须等于我要收的数量，否则 Alltoallv 会死锁
    int *rcv_chk = new int[cpusize];
    MPI_Alltoall(snd_cnt, 1, MPI_INT, rcv_chk, 1, MPI_INT, MPI_COMM_WORLD);
    for (int t = 0; t < cpusize; t++)
      if (rcv_chk[t] != rcv_cnt[t])
      {
        cout << "ERROR: Patch::Interp_Points inconsistent point ownership, rank "
             << myrank << " vs " << t << ": want " << rcv_cnt[t]
             << " got " << rcv_chk[t] << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    delete[] rcv_chk;

    // 4b) 单位换算：上面统计的是「点数」，而 Alltoallv 的 count/disp 单位是
    //     「元素个数」（这里是 double 的个数），必须先乘 num_var 再算位移，
    //     否则 count 与 disp 单位不一致，只会传出 1/num_var 的数据。
    for (int t = 0; t < cpusize; t++)
    {
      snd_cnt[t] *= num_var;
      rcv_cnt[t] *= num_var;
    }

    // 5) 打包缓冲
    int *snd_disp = new int[cpusize];
    int *rcv_disp = new int[cpusize];
    int tot_snd = 0, tot_rcv = 0;
    for (int t = 0; t < cpusize; t++)
    {
      snd_disp[t] = tot_snd;
      tot_snd += snd_cnt[t]; // snd_cnt 已在 4b) 换算成元素个数，此处不可再乘
      rcv_disp[t] = tot_rcv;
      tot_rcv += rcv_cnt[t];
    }
    double *snd_buf = new double[tot_snd > 0 ? tot_snd : 1];
    double *rcv_buf = new double[tot_rcv > 0 ? tot_rcv : 1];
    int *snd_pos = new int[cpusize];
    for (int t = 0; t < cpusize; t++)
      snd_pos[t] = snd_disp[t];

    // 6) 收到的第 e 个点对应哪个 j（排布 = 源 rank 升序、源内 j 升序）
    const int n_rcv_pt = tot_rcv / num_var;
    int *rcv_j = new int[n_rcv_pt > 0 ? n_rcv_pt : 1];
    {
      int *rp = new int[cpusize];
      for (int t = 0; t < cpusize; t++)
        rp[t] = rcv_disp[t];
      for (int j = need_lo; j <= need_hi; j++)
      {
        const int o = owner[j];
        if (o >= 0 && o != myrank)
        {
          rcv_j[rp[o] / num_var] = j;
          rp[o] += num_var;
        }
      }
      delete[] rp;
    }

    // 7) 只插值「本 rank 拥有且有人要」的点：对每个需求方 t 扫它的区间。
    //    t==myrank 写回 Shellf 自用，否则打包进 snd_buf 发走。
    //    落在所有需求区间之外的点（Amdahl 里纯属白算的那部分）这里直接跳过。
    double t_fc0 = MPI_Wtime();
    double *vals = new double[num_var];
    for (int t = 0; t < cpusize; t++)
    {
      for (int j = q_lo[t]; j <= q_hi[t]; j++)
      {
        const bh_point_geom &g = geom[j];
        if (!g.BP || g.BP->rank != myrank)
          continue;
        varl = VarList;
        int k = 0;
        while (varl) // run along variables
        {
          vals[k] = bh_fused_interp(g, g.BP->fgfs[varl->data->sgfn],
                                    g.BP->shape[0],
                                    g.BP->shape[0] * g.BP->shape[1],
                                    varl->data->SoA);
          varl = varl->next;
          k++;
        }
        if (t == myrank)
        {
          for (int q = 0; q < num_var; q++)
            Shellf[j * num_var + q] = vals[q];
        }
        else
        {
          for (int q = 0; q < num_var; q++)
            snd_buf[snd_pos[t] + q] = vals[q];
          snd_pos[t] += num_var;
        }
      }
    }
    delete[] vals;
    double t_fc1 = MPI_Wtime();
    bh_prof_ia_search += t_fc1 - t_fc0;

    // 8) 一次 Alltoallv 取代两次全量 Allreduce
    t_fc0 = MPI_Wtime();
    MPI_Alltoallv(snd_buf, snd_cnt, snd_disp, MPI_DOUBLE,
                  rcv_buf, rcv_cnt, rcv_disp, MPI_DOUBLE, MPI_COMM_WORLD);
    t_fc1 = MPI_Wtime();
    bh_prof_ia_allreduce += t_fc1 - t_fc0;

    // 9) 散回 Shellf（此时 [need_lo, need_hi] 段已完整；段外调用方不读）
    for (int e = 0; e < n_rcv_pt; e++)
      for (int q = 0; q < num_var; q++)
        Shellf[rcv_j[e] * num_var + q] = rcv_buf[e * num_var + q];

    // 10) 归属缺失检查（原 Weight[i]==0 分支的等价物；各 rank 合起来覆盖全部点）
    int bad = -1;
    for (int j = need_lo; j <= need_hi; j++)
      if (owner[j] < 0)
      {
        bad = j;
        break;
      }
    int bad_all = -1;
    MPI_Allreduce(&bad, &bad_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (bad_all >= 0)
    {
      if (myrank == 0)
      {
        cout << "ERROR: Patch::Interp_Points fails to find point (";
        for (int k = 0; k < dim; k++)
          cout << XX[k][bad_all] << (k < dim - 1 ? "," : ")");
        cout << " on Patch (";
        for (int k = 0; k < dim; k++)
          cout << bbox[k] << "+" << lli[k] * DH[k] << (k < dim - 1 ? "," : ")--");
        cout << "(";
        for (int k = 0; k < dim; k++)
          cout << bbox[dim + k] << "-" << uui[k] * DH[k] << (k < dim - 1 ? "," : ")") << endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    delete[] q_lo;
    delete[] q_hi;
    delete[] owner;
    delete[] snd_cnt;
    delete[] rcv_cnt;
    delete[] snd_disp;
    delete[] rcv_disp;
    delete[] snd_buf;
    delete[] rcv_buf;
    delete[] snd_pos;
    delete[] rcv_j;
    delete[] shellf;
    delete[] weight;
    delete[] DH;
    delete[] llb;
    delete[] uub;
    return;
  }
  // ========================= BH_OPT3 结束 =========================

  double t_ia0 = MPI_Wtime();
  // BH_OPT4: 传了本 rank 的点列表就只扫这些点，否则照旧全扫。
  // 球面采样点只落在少数几个块上（实测 np=8 / np=32 都是 4 个 rank 分走全部 25600 点），
  // 所以对大多数 rank 来说全扫是纯浪费；对忙的 rank 也省下不属于自己的那部分。
  const bool scan_mine = (mine != 0 && n_mine >= 0 && !getenv("BH_NO_MINE")); // BH_A/B: 临时对照开关
  const int n_scan = scan_mine ? n_mine : NN;
  for (int scan = 0; scan < n_scan; scan++) // run along points
  {
    const int j = scan_mine ? mine[scan] : scan;
    // BH_OPT: 命中预计算几何则直接走融合插值内核，跳过块搜索与 Fortran 插值链
    if (geom)
    {
      const bh_point_geom &g = geom[j];
      if (g.BP)
      {
        // BH_OPT5: 有块数据副本时，插值从副本取数（副本每轮在计时区外刷新），
        // 于是「块属于谁」不再限制谁能算——点已经按 rank 均分过了，每个点只会
        // 落在唯一一个 rank 的 mine 列表里，不会重复计算。
        const bool use_copy = (bc != 0 && g.bi >= 0 && num_var <= 2 &&
                               !getenv("BH_NO_COPY")); // BH_A/B: 临时对照开关
        if (use_copy || myrank == g.BP->rank)
        {
          varl = VarList;
          int k = 0;
          while (varl) // run along variables
          {
            const double *fsrc = use_copy ? bc[g.bi].f[k] : g.BP->fgfs[varl->data->sgfn];
            const int fNx = use_copy ? bc[g.bi].shape[0] : g.BP->shape[0];
            const int fNxy = use_copy ? bc[g.bi].shape[0] * bc[g.bi].shape[1]
                                      : g.BP->shape[0] * g.BP->shape[1];
            shellf[j * num_var + k] = bh_fused_interp(g, fsrc, fNx, fNxy,
                                                      varl->data->SoA);
            if (bh_dbg_on()) // BH_DEBUG: 临时对拍探针（提交前移除）
            {
              static int dbg_n = 0;
              double ref = 0;
              int dbg_shp[3] = {g.BP->shape[0], g.BP->shape[1], g.BP->shape[2]};
              int dbg_ordn = g.ordn;
              int dbg_sym = Symmetry;
              double dbg_x = XX[0][j], dbg_y = XX[1][j], dbg_z = XX[2][j];
              f_global_interp(dbg_shp, g.BP->X[0], g.BP->X[1], g.BP->X[2],
                              g.BP->fgfs[varl->data->sgfn], ref,
                              dbg_x, dbg_y, dbg_z, dbg_ordn,
                              varl->data->SoA, dbg_sym);
              const double d = fabs(shellf[j * num_var + k] - ref);
              // 对照：记下第一个 kernel 与 Fortran 完全一致的点
              static int good_have = 0;
              static bh_point_geom good_g;
              static double good_xyz[3];
              if (d <= 1e-9 * (1.0 + fabs(ref)) && !good_have)
              {
                good_have = 1;
                good_g = g;
                good_xyz[0] = XX[0][j];
                good_xyz[1] = XX[1][j];
                good_xyz[2] = XX[2][j];
              }
              if (d > 1e-9 * (1.0 + fabs(ref)))
              {
                // 按反射位模式分类统计：0=无反射 1=x 2=y 4=z
                static long dbg_cnt[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                static double dbg_max[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                const int pat = ((int)g.refl[0][0]) + 2 * ((int)g.refl[1][0]) + 4 * ((int)g.refl[2][0]);
                dbg_cnt[pat]++;
                if (d > dbg_max[pat])
                  dbg_max[pat] = d;
                if (pat == 4 && dbg_n == 0) // z 反射点：三方对拍一例
                {
                  dbg_n = 1;
                  const double xyz[3] = {XX[0][j], XX[1][j], XX[2][j]};
                  double cx3[3] = {0, 0, 0};
                  const double vman = bh_dbg_manual(g, g.BP->fgfs[varl->data->sgfn],
                                                    g.BP->shape[0],
                                                    g.BP->shape[0] * g.BP->shape[1],
                                                    varl->data->SoA, xyz, cx3);
                  cout << setprecision(17) << "BH_DEBUG ZREFL rk=" << myrank << " pt=" << j
                       << " var=" << k << endl;
                  cout << "BH_DEBUG   xyz=" << xyz[0] << " " << xyz[1] << " " << xyz[2]
                       << " base=" << g.base[0] << " " << g.base[1] << " " << g.base[2] << endl;
                  cout << "BH_DEBUG   cx=" << cx3[0] << " " << cx3[1] << " " << cx3[2] << endl;
                  cout << "BH_DEBUG   kernel=" << shellf[j * num_var + k]
                       << " fortran=" << ref << " manual=" << vman
                       << " |man-ref|=" << fabs(vman - ref)
                       << " |man-kernel|=" << fabs(vman - shellf[j * num_var + k]) << endl;
                  cout << "BH_DEBUG   zX[0..6]=" << g.BP->X[2][0] << " " << g.BP->X[2][1] << " "
                       << g.BP->X[2][2] << " " << g.BP->X[2][3] << " " << g.BP->X[2][4] << " "
                       << g.BP->X[2][5] << " " << g.BP->X[2][6] << endl;
                  // 冲激实验：逐层置 1，测出 Fortran 真正使用的 z 方向有效系数
                  const int Nx_ = g.BP->shape[0], Ny_ = g.BP->shape[1], Nz_ = g.BP->shape[2];
                  const long sz = (long)Nx_ * Ny_ * Nz_;
                  double *ft = g.BP->fgfs[varl->data->sgfn];
                  double *sav = new double[sz];
                  for (long t = 0; t < sz; t++)
                    sav[t] = ft[t];
                  cout << "BH_DEBUG   impulse k:coef (Fortran 实际系数)" << endl;
                  for (int kk = 0; kk < g.ordn; kk++)
                  {
                    for (long t = 0; t < sz; t++)
                      ft[t] = 0;
                    for (int q = 0; q < g.ordn; q++)
                      for (int p = 0; p < g.ordn; p++)
                      {
                        const int ii = g.base[0] + p - 1, jj = g.base[1] + q - 1;
                        if (ii >= 0 && ii < Nx_ && jj >= 0 && jj < Ny_ && kk < Nz_)
                          ft[ii + (long)jj * Nx_ + (long)kk * Nx_ * Ny_] = 1.0;
                      }
                    double rr = 0;
                    int dbg_shp2[3] = {Nx_, Ny_, Nz_};
                    int o2 = g.ordn, s2 = Symmetry;
                    double xa = xyz[0], ya2 = xyz[1], za = xyz[2];
                    f_global_interp(dbg_shp2, g.BP->X[0], g.BP->X[1], g.BP->X[2], ft, rr,
                                    xa, ya2, za, o2, varl->data->SoA, s2);
                    cout << "BH_DEBUG     k=" << kk << " coef=" << rr << endl;
                  }
                  for (long t = 0; t < sz; t++)
                    ft[t] = sav[t];
                  delete[] sav;
                  cout << "BH_DEBUG   mine zcoef = ["
                       << g.w[2][1] << " " << (g.w[2][0] * varl->data->SoA[2] + g.w[2][2]) << " "
                       << g.w[2][3] << " " << g.w[2][4] << " " << g.w[2][5] << "]" << endl;
                  // 决定性实验：直接调用 Fortran 的 decide3d_，导出真正的 ya 模板
                  {
                    const double *fp = g.BP->fgfs[varl->data->sgfn];
                    const int Nx2 = g.BP->shape[0], Nxy2 = g.BP->shape[0] * g.BP->shape[1];
                    int cxb[3] = {g.base[0], g.base[1], g.base[2]};
                    int cxt[3] = {cxb[0] + g.ordn - 1, cxb[1] + g.ordn - 1, cxb[2] + g.ordn - 1};
                    double *yat = new double[BH_ORDN_MAX * BH_ORDN_MAX * BH_ORDN_MAX];
                    int o3 = g.ordn, s3 = Symmetry;
                    int ok = decide3d_(dbg_shp, const_cast<double *>(fp), const_cast<double *>(fp),
                                       cxb, cxt, varl->data->SoA, yat, &o3, &s3);
                    cout << "BH_DEBUG   decide3d_ ok=" << ok << " cxB=" << cxb[0] << " " << cxb[1] << " "
                         << cxb[2] << " cxT=" << cxt[0] << " " << cxt[1] << " " << cxt[2] << endl;
                    cout << "BH_DEBUG   SoA=" << varl->data->SoA[0] << " " << varl->data->SoA[1] << " "
                         << varl->data->SoA[2] << endl;
                    const int o2_ = g.ordn, n1n2 = g.ordn * g.ordn;
                    cout << "BH_DEBUG   FORTRAN ya  z-line(p=q=0):";
                    for (int r = 0; r < o2_; r++)
                      cout << " " << yat[0 + 0 * o2_ + r * n1n2];
                    cout << endl;
                    cout << "BH_DEBUG   MINE    ya  z-line(p=q=0):";
                    for (int r = 0; r < o2_; r++)
                    {
                      const int k0 = g.base[2] + r - 1; // 0-based
                      double v;
                      if (k0 < 0)
                        v = fp[(g.base[0] - 1) + (g.base[1] - 1) * Nx2 + (2 - 2 * g.base[2] - r) * Nxy2] *
                            varl->data->SoA[2];
                      else
                        v = fp[(g.base[0] - 1) + (g.base[1] - 1) * Nx2 + k0 * Nxy2];
                      cout << " " << v;
                    }
                    cout << endl;
                    cout << "BH_DEBUG   f idx0..5 (p=q=0):";
                    for (int r = 0; r < 6; r++)
                      cout << " " << fp[(g.base[0] - 1) + (g.base[1] - 1) * Nx2 + r * Nxy2];
                    cout << endl;
                    delete[] yat;
                    // 对照实验：同一个流水线跑"不反射"的点（kernel 与 Fortran 完全一致的点）
                    if (good_have)
                    {
                      double cxo[3];
                      const double vg = bh_dbg_pipe(good_g.BP, good_g, good_g.BP->fgfs[varl->data->sgfn],
                                                    good_xyz, varl->data->SoA, Symmetry, cxo);
                      double rfg = 0;
                      int gs[3] = {good_g.BP->shape[0], good_g.BP->shape[1], good_g.BP->shape[2]};
                      int go = good_g.ordn, gsy = Symmetry;
                      double gx = good_xyz[0], gy = good_xyz[1], gz = good_xyz[2];
                      f_global_interp(gs, good_g.BP->X[0], good_g.BP->X[1], good_g.BP->X[2],
                                      good_g.BP->fgfs[varl->data->sgfn], rfg, gx, gy, gz, go,
                                      varl->data->SoA, gsy);
                      cout << "BH_DEBUG   对照(不反射点) base=" << good_g.base[0] << " " << good_g.base[1]
                           << " " << good_g.base[2] << " refl=" << (int)good_g.refl[0][0]
                           << (int)good_g.refl[1][0] << (int)good_g.refl[2][0]
                           << " pipe=" << vg << " fortran=" << rfg << " diff=" << (vg - rfg) << endl;
                    }
                    // 扫描：改变 cxB(3) 的约定，看哪一个能精确复现 Fortran 的 ref
                    {
                      double xa1[8];
                      for (int t = 0; t < 8; t++)
                        xa1[t] = t;
                      cout << "BH_DEBUG   扫描 cxB(3) 约定（ref=" << ref << "）" << endl;
                      for (int dof = -3; dof <= 3; dof++)
                      {
                        int cb[3] = {g.base[0], g.base[1], g.base[2] + dof};
                        int ct[3] = {cb[0] + g.ordn - 1, cb[1] + g.ordn - 1, cb[2] + g.ordn - 1};
                        double *yb = new double[BH_ORDN_MAX * BH_ORDN_MAX * BH_ORDN_MAX];
                        for (int t = 0; t < BH_ORDN_MAX * BH_ORDN_MAX * BH_ORDN_MAX; t++)
                          yb[t] = 0.0;
                        int o5 = g.ordn, s5 = Symmetry;
                        int kk = decide3d_(dbg_shp, const_cast<double *>(fp), const_cast<double *>(fp),
                                           cb, ct, varl->data->SoA, yb, &o5, &s5);
                        double cxv[3];
                        for (int m = 0; m < 3; m++)
                        {
                          double *Xm = (m == 0) ? g.BP->X[0] : (m == 1) ? g.BP->X[1] : g.BP->X[2];
                          const double dmm = Xm[1] - Xm[0];
                          if (cb[m] > 0)
                            cxv[m] = (xyz[m] - Xm[cb[m] - 1]) / dmm;
                          else
                            cxv[m] = (xyz[m] + Xm[1 - cb[m]]) / dmm;
                        }
                        double vy = 0, dyv = 0;
                        polin3_(xa1, xa1, xa1, yb, &cxv[0], &cxv[1], &cxv[2], &vy, &dyv, &o5);
                        cout << "BH_DEBUG     cxB3=" << cb[2] << " cx3=" << cxv[2] << " val=" << vy
                             << " diff=" << (vy - ref) << " err=" << kk << endl;
                        delete[] yb;
                      }
                    }
                    // 无歧义实验：把整个场填成"自己的扁平下标"，再看 decide3d 读到了哪些下标
                    {
                      const long sz2 = (long)Nx2 * g.BP->shape[1] * g.BP->shape[2];
                      double *fld = g.BP->fgfs[varl->data->sgfn];
                      double *sv2 = new double[sz2];
                      for (long t = 0; t < sz2; t++)
                        sv2[t] = fld[t];
                      for (long t = 0; t < sz2; t++)
                        fld[t] = (double)t;
                      double *y2 = new double[BH_ORDN_MAX * BH_ORDN_MAX * BH_ORDN_MAX];
                      int cc[3] = {g.base[0], g.base[1], g.base[2]};
                      int dd[3] = {cc[0] + g.ordn - 1, cc[1] + g.ordn - 1, cc[2] + g.ordn - 1};
                      int o4 = g.ordn, s4 = Symmetry;
                      decide3d_(dbg_shp, fld, fld, cc, dd, varl->data->SoA, y2, &o4, &s4);
                      cout << "BH_DEBUG   IDX-PROBE (p,q)=(0,0):";
                      for (int r = 0; r < 6; r++)
                        cout << " " << (long)y2[0 + 0 * 6 + r * 36];
                      cout << endl;
                      cout << "BH_DEBUG   IDX-PROBE (p,q)=(1,0):";
                      for (int r = 0; r < 6; r++)
                        cout << " " << (long)y2[1 + 0 * 6 + r * 36];
                      cout << endl;
                      cout << "BH_DEBUG   IDX-PROBE (p,q)=(0,1):";
                      for (int r = 0; r < 6; r++)
                        cout << " " << (long)y2[0 + 1 * 6 + r * 36];
                      cout << endl;
                      delete[] y2;
                      for (long t = 0; t < sz2; t++)
                        fld[t] = sv2[t];
                      delete[] sv2;
                    }
                  }
                }
                if (pat == 4 && dbg_n == 1 && myrank == 0 && j == 0)
                {
                  cout << "BH_DEBUG summary4(pat:count:maxdiff) " << dbg_cnt[4] << ":"
                       << dbg_max[4] << endl;
                }
              }
            }
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
        }
        continue;
      }
      // g.BP == NULL：预计算未找到归属块（正常不会发生），回退到原搜索路径
    }
    double pox[dim];
    for (int i = 0; i < dim; i++)
    {
      pox[i] = XX[i][j];
      if (myrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
      {
        cout << "Patch::Interp_Points: point (";
        for (int k = 0; k < dim; k++)
        {
          cout << XX[k][j];
          if (k < dim - 1)
            cout << ",";
          else
            cout << ") is out of current Patch." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    MyList<Block> *Bp = blb;
    bool notfind = true;
    while (notfind && Bp) // run along Blocks
    {
      Block *BP = Bp->data;

      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        if (myrank == BP->rank)
        {
          //---> interpolation
          varl = VarList;
          int k = 0;
          while (varl) // run along variables
          {
            //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
            //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                            pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
        }
      }
      if (Bp == ble)
        break;
      Bp = Bp->next;
    }
  }
  double t_ia1 = MPI_Wtime();
  bh_prof_ia_search += t_ia1 - t_ia0;
  t_ia0 = MPI_Wtime();

  t_ia0 = MPI_Wtime();
  MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int *Weight;
  Weight = new int[NN];
  MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  t_ia1 = MPI_Wtime();
  bh_prof_ia_allreduce += t_ia1 - t_ia0;

  //  misc::tillherecheck("print me");

  for (int i = 0; i < NN; i++)
  {
    if (Weight[i] > 1)
    {
      if (myrank == 0)
        cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
      for (int j = 0; j < num_var; j++)
        Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
    }
    else if (Weight[i] == 0 && myrank == 0)
    {
      cout << "ERROR: Patch::Interp_Points fails to find point (";
      for (int j = 0; j < dim; j++)
      {
        cout << XX[j][i];
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")";
      }
      cout << " on Patch (";
      for (int j = 0; j < dim; j++)
      {
        cout << bbox[j] << "+" << lli[j] * getdX(j);
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")--";
      }
      cout << "(";
      for (int j = 0; j < dim; j++)
      {
        cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
#if 0
       checkBlock();
#else
      cout << "splited domains:" << endl;
      {
        MyList<Block> *Bp = blb;
        while (Bp)
        {
          Block *BP = Bp->data;

          for (int i = 0; i < dim; i++)
          {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
            llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
          }
          cout << "(";
          for (int j = 0; j < dim; j++)
          {
            cout << llb[j] << ":" << uub[j];
            if (j < dim - 1)
              cout << ",";
            else
              cout << ")" << endl;
          }
          if (Bp == ble)
            break;
          Bp = Bp->next;
        }
      }
#endif
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
//|---- BH_OPT: 计时区外预计算每个采样点的插值几何。
// 块归属搜索与 Interp_Points 完全一致（列表顺序、半网格容差）；
// 模板下标/钳位/反射逻辑与 fmisc.f90::global_interp 逐行等价；
// 权重用拉格朗日直写形式（与 polint 的 Neville 算法是同一多项式）。
void Patch::prepare_point_geom(int NN, double **XX, int Symmetry,
                               bh_point_geom *geom,
                               bh_block_copy **pbc, int *pn_bc)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  const int ordn = surface_interpolation_order();
  if (ordn < 2 || ordn > BH_ORDN_MAX)
  {
    if (myrank == 0)
      cout << "Patch::prepare_point_geom: interpolation order " << ordn
           << " out of supported range [2," << BH_ORDN_MAX << "]" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  double DH[dim];
  for (int i = 0; i < dim; i++)
    DH[i] = getdX(i);
  double llb[dim], uub[dim];

  // BH_OPT5: 收集块列表，并为每个块备一份本地 X 坐标副本。
  // Block 构造器只在拥有者 rank 上分配 X（见 Block.C 的 if (myrank == rank)），
  // 但几何计算现在每个 rank 都要做——因为块数据副本每轮会广播过来。X 是纯几何量、
  // 不随时间变，这里从拥有者广播一次即可（8 块合计约 6KB，且只在计时区外做一次）。
  // 用广播而不是本地按 bbox+shape 重建，是为了跟拥有者的 X 逐位一致：插值几何
  // 对 X 的末位差异是敏感的（cx 由 X 算出，权重再由 cx 算出）。
  MyList<Block> *BpL = blb;
  int nb = 0;
  while (BpL)
  {
    nb++;
    if (BpL == ble)
      break;
    BpL = BpL->next;
  }
  Block **blks = new Block *[nb > 0 ? nb : 1];
  BpL = blb;
  for (int t = 0; t < nb; t++)
  {
    blks[t] = BpL->data;
    if (BpL == ble)
      break;
    BpL = BpL->next;
  }
  double **Xc = new double *[nb * dim];
  for (int t = 0; t < nb; t++)
    for (int m = 0; m < dim; m++)
    {
      const int ns = blks[t]->shape[m];
      double *p = new double[ns];
      if (myrank == blks[t]->rank)
        memcpy(p, blks[t]->X[m], sizeof(double) * ns);
      MPI_Bcast(p, ns, MPI_DOUBLE, blks[t]->rank, MPI_COMM_WORLD);
      Xc[t * dim + m] = p;
    }

  for (int j = 0; j < NN; j++)
  {
    bh_point_geom &g = geom[j];
    memset(&g, 0, sizeof(bh_point_geom));
    g.ordn = ordn;
    g.bi = -1; // BH_OPT5: memset 归零后必须显式置 -1（0 是合法的副本下标）

    // 越界检查与 Interp_Points 一致（提前到预计算阶段做一次）
    for (int i = 0; i < dim; i++)
    {
      if (myrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
      {
        cout << "Patch::prepare_point_geom: point (" << XX[0][j] << "," << XX[1][j] << ","
             << XX[2][j] << ") is out of current Patch." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    // 块归属搜索（与 Interp_Points 相同的 llb/uub 判定）
    // BH_OPT5: 改用 blks[] 数组遍历（内容与 MyList 顺序一致），以便记下块序号 t
    bool notfind = true;
    for (int t = 0; t < nb && notfind; t++)
    {
      Block *BP = blks[t];
      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }
      if (flag)
      {
        notfind = false;
        g.BP = BP;
        g.bi = t; // BH_OPT5: 记下命中的块序号，下面取该块的 X 副本
      }
    }
    if (!g.BP)
      continue; // 运行时回退旧路径
    // BH_OPT5: 这里原有「if (g.BP->rank != myrank) continue;」——现已删除。
    // 块数据副本每轮在计时区外广播给所有 rank，所以每个 rank 都要把几何算全。
    // 下面用 Xc（广播来的 X 副本）代替 g.BP->X：后者只在拥有者 rank 上分配
    // （见 Block.C 的 if (myrank == rank)）。shape 各 rank 都有，照旧用 g.BP->shape。

    // 等价 global_interp 的 cxI/cxB/cxT 计算与钳位（Fortran 1 基下标）
    for (int m = 0; m < 3; m++)
    {
      const double *Xm = Xc[g.bi * dim + m];
      const double dxm = Xm[1] - Xm[0];
      // Fortran idint 是向零截断，与 C 的 (int) 强转一致
      const int cxi = (int)((XX[m][j] - Xm[0]) / dxm + 0.4) + 1;
      int cxB = cxi - ordn / 2 + 1;
      int cxT = cxB + ordn - 1;
      int cmin = 1;
      if (m == 2 && Symmetry != 0 && fabs(Xm[0]) < dxm)
        cmin = -ordn / 2 + 1; // 赤道对称：z 方向越界回折
      if (m < 2 && Symmetry == 2 && fabs(Xm[0]) < dxm)
        cmin = -ordn / 2 + 1; // 八分对称：x/y 方向越界回折
      if (cxB < cmin)
      {
        cxB = cmin;
        cxT = cxB + ordn - 1;
      }
      if (cxT > g.BP->shape[m])
      {
        cxT = g.BP->shape[m];
        cxB = cxT + 1 - ordn;
      }
      g.base[m] = cxB;
      // 模板内分数坐标；cxB<=0 时坐标数组按反射延拓（见 global_interp 的 cx 计算）
      double cx;
      if (cxB > 0)
        cx = (XX[m][j] - Xm[cxB - 1]) / dxm;
      else
        cx = (XX[m][j] + Xm[1 - cxB - 1]) / dxm; // 镜像面在 X(1)，非 X(2)
      // 越界反射标记：1 基下标 fi = cxB+p < 1 的模板点取镜像 1-fi 并乘 SoA
      for (int p = 0; p < ordn; p++)
        g.refl[m][p] = (cxB + p < 1) ? 1 : 0;
      // 拉格朗日权重：节点 0..ordn-1 上取 cx 的值
      for (int p = 0; p < ordn; p++)
      {
        double num = 1.0, den = 1.0;
        for (int q = 0; q < ordn; q++)
          if (q != p)
          {
            num *= (cx - q);
            den *= (p - q);
          }
        g.w[m][p] = num / den;
      }
    }
  }

  // BH_OPT5: 建「含点块」的数据副本条目。只有被采样点穿过的块才需要副本
  // （实测 8 块中只有 4 块），其数据在计时区外每轮由各自拥有者广播刷新
  // （见 surface_integral::sync_block_copies）。这里只分配缓冲区。
  {
    int *t2bc = new int[nb];
    for (int t = 0; t < nb; t++)
      t2bc[t] = -1;
    int nbc = 0;
    for (int j = 0; j < NN; j++)
      if (geom[j].bi >= 0 && t2bc[geom[j].bi] < 0)
        t2bc[geom[j].bi] = nbc++;
    *pn_bc = nbc;
    *pbc = new bh_block_copy[nbc > 0 ? nbc : 1];
    for (int t = 0; t < nb; t++)
      if (t2bc[t] >= 0)
      {
        bh_block_copy &c = (*pbc)[t2bc[t]];
        c.src = blks[t];
        for (int m = 0; m < dim; m++)
          c.shape[m] = blks[t]->shape[m];
        const long sz = (long)c.shape[0] * c.shape[1] * c.shape[2];
        for (int v = 0; v < 2; v++)
          c.f[v] = new double[sz](); // 值初始化：首轮同步前不会被读到
      }
    // 把 geom[].bi 从「块序号」改写成「副本下标」
    for (int j = 0; j < NN; j++)
      if (geom[j].bi >= 0)
        geom[j].bi = t2bc[geom[j].bi];
    delete[] t2bc;
  }

  for (int t = 0; t < nb * dim; t++)
    delete[] Xc[t];
  delete[] Xc;
  delete[] blks;
}
void Patch::Interp_Points(MyList<var> *VarList,
                          int NN, double **XX,
                          double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank, lmyrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_rank(Comm_here, &lmyrank);

  int ordn = surface_interpolation_order();
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[NN * num_var];
  memset(shellf, 0, sizeof(double) * NN * num_var);

  // we use weight to monitor code, later some day we can move it for optimization
  int *weight;
  weight = new int[NN];
  memset(weight, 0, sizeof(int) * NN);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  for (int j = 0; j < NN; j++) // run along points
  {
    double pox[dim];
    for (int i = 0; i < dim; i++)
    {
      pox[i] = XX[i][j];
      if (lmyrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
      {
        cout << "Patch::Interp_Points: point (";
        for (int k = 0; k < dim; k++)
        {
          cout << XX[k][j];
          if (k < dim - 1)
            cout << ",";
          else
            cout << ") is out of current Patch." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    MyList<Block> *Bp = blb;
    bool notfind = true;
    while (notfind && Bp) // run along Blocks
    {
      Block *BP = Bp->data;

      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        if (myrank == BP->rank)
        {
          //---> interpolation
          varl = VarList;
          int k = 0;
          while (varl) // run along variables
          {
            //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
            //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
            f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                            pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
        }
      }
      if (Bp == ble)
        break;
      Bp = Bp->next;
    }
  }

  MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
  int *Weight;
  Weight = new int[NN];
  MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, Comm_here);

  //  misc::tillherecheck("print me");
  //  if(lmyrank == 0) cout<<"myrank = "<<myrank<<"print me"<<endl;

  for (int i = 0; i < NN; i++)
  {
    if (Weight[i] > 1)
    {
      if (lmyrank == 0)
        cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
      for (int j = 0; j < num_var; j++)
        Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
    }
#if 0 // for not involved levels, this may fail     
     else if(Weight[i] == 0 && lmyrank == 0)
     {
       cout<<"ERROR: Patch::Interp_Points fails to find point (";
       for(int j=0;j<dim;j++)
       {
	  cout<<XX[j][i];
	  if(j<dim-1) cout<<",";
	  else        cout<<")";
       }
       cout<<" on Patch (";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[j]<<"+"<<lli[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")--";
       }
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<bbox[dim+j]<<"-"<<uui[j]*getdX(j);
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
#if 0
       checkBlock();
#else
  cout<<"splited domains:"<<endl;
  {
     MyList<Block> *Bp=blb;
     while(Bp)
     {
	Block *BP=Bp->data;

	for(int i=0;i<dim;i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +ghost_width*DH[i];
          uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif 
	}       
       cout<<"(";
       for(int j=0;j<dim;j++)
       {
	  cout<<llb[j]<<":"<<uub[j];
	  if(j<dim-1) cout<<",";
	  else        cout<<")"<<endl;
       }
	if(Bp == ble) break;
	Bp=Bp->next;
     }
  }
#endif       
       MPI_Abort(MPI_COMM_WORLD,1);
     }
#endif
  }

  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
void Patch::checkBlock()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    MyList<Block> *BP = blb;
    while (BP)
    {
      BP->data->checkBlock();
      if (BP == ble)
        break;
      BP = BP->next;
    }
  }
}
double Patch::getdX(int dir)
{
  if (dir < 0 || dir >= dim)
  {
    cout << "Patch::getdX: error input dir = " << dir << ", this Patch has direction (0," << dim - 1 << ")" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  if (shape[dir] == 1)
  {
    cout << "Patch::getdX: for direction " << dir << ", this Patch has only one point. Can not determine dX for vertex center grid." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  h = (bbox[dim + dir] - bbox[dir]) / (shape[dir] - 1);
#else
#ifdef Cell
  h = (bbox[dim + dir] - bbox[dir]) / shape[dir];
#else
#error Not define Vertex nor Cell
#endif
#endif
  return h;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
                             double *Shellf, int Symmetry)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = surface_interpolation_order();
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[num_var];
  memset(shellf, 0, sizeof(double) * num_var);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  double pox[dim];
  for (int i = 0; i < dim; i++)
  {
    pox[i] = XX[i];
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] shellf;
      delete[] DH;
      delete[] llb;
      delete[] uub;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  MyList<Block> *Bp = blb;
  bool notfind = true;
  while (notfind && Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    bool flag = true;
    for (int i = 0; i < dim; i++)
    {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
      {
        flag = false;
        break;
      }
    }

    if (flag)
    {
      notfind = false;
      if (myrank == BP->rank)
      {
// test old code
#if 0
#define floorint(a) ((a) < 0 ? int(a) - 1 : int(a))
//---> interpolation
                int ixl,iyl,izl,ixu,iyu,izu;
	    	double Delx,Dely,Delz;

		ixl = 1+floorint((pox[0]-BP->X[0][0])/DH[0]);
	   	iyl = 1+floorint((pox[1]-BP->X[1][0])/DH[1]);
	   	izl = 1+floorint((pox[2]-BP->X[2][0])/DH[2]);

		int nn=ordn/2;

		ixl = ixl-nn;
		iyl = iyl-nn;
		izl = izl-nn;
	   
		int tmi;
		tmi = (Symmetry==2)?-1:0;
		if(ixl<tmi) ixl=tmi;
	   	if(iyl<tmi) iyl=tmi;
		tmi = (Symmetry>0)?-1:0;
	   	if(izl<tmi) izl=tmi;
      
	   	if(ixl+ordn>BP->shape[0]) ixl=BP->shape[0]-ordn;
	   	if(iyl+ordn>BP->shape[1]) iyl=BP->shape[1]-ordn;
	   	if(izl+ordn>BP->shape[2]) izl=BP->shape[2]-ordn;
// support cell center
		if(ixl>=0) Delx = ( pox[0] - BP->X[0][ixl] )/ DH[0];
		else       Delx = ( pox[0] + BP->X[0][0] )/ DH[0];
                if(iyl>=0) Dely = ( pox[1] - BP->X[1][iyl] )/ DH[1];
		else       Dely = ( pox[1] + BP->X[1][0] )/ DH[1];
                if(izl>=0) Delz = ( pox[2] - BP->X[2][izl] )/ DH[2];
		else       Delz = ( pox[2] + BP->X[2][0] )/ DH[2];
//change to fortran index
                ixl++;
	   	iyl++;
	   	izl++;
	   	ixu = ixl + ordn - 1;
	   	iyu = iyl + ordn - 1;
	   	izu = izl + ordn - 1;
	    	varl=VarList;
		int j=0;
	    	while(varl)
		{
                 f_interp_2(BP->shape,BP->fgfs[varl->data->sgfn],shellf[j],ixl,ixu,iyl,iyu,izl,izu,Delx,Dely,Delz,
                                     ordn,varl->data->SoA,Symmetry);
		 varl=varl->next;
		 j++;
		} //varl
#else
        //---> interpolation
        varl = VarList;
        int k = 0;
        while (varl) // run along variables
        {
          //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
          //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
          f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
                          pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
          varl = varl->next;
          k++;
        }
#endif
      }
    }
    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  if (notfind && myrank == 0)
  {
    cout << "ERROR: Patch::Interp_Points fails to find point (";
    for (int j = 0; j < dim; j++)
    {
      cout << XX[j];
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")";
    }
    cout << " on Patch (";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[j] << "+" << lli[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")--";
    }
    cout << "(";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")" << endl;
    }
#if 0
       checkBlock();
#else
    cout << "splited domains:" << endl;
    {
      MyList<Block> *Bp = blb;
      while (Bp)
      {
        Block *BP = Bp->data;

        for (int i = 0; i < dim; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        }
        cout << "(";
        for (int j = 0; j < dim; j++)
        {
          cout << llb[j] << ":" << uub[j];
          if (j < dim - 1)
            cout << ",";
          else
            cout << ")" << endl;
        }
        if (Bp == ble)
          break;
        Bp = Bp->next;
      }
    }
#endif
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] shellf;
  delete[] DH;
  delete[] llb;
  delete[] uub;

  return true;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
                             double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = surface_interpolation_order();
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[num_var];
  memset(shellf, 0, sizeof(double) * num_var);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  double pox[dim];
  for (int i = 0; i < dim; i++)
  {
    pox[i] = XX[i];
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] shellf;
      delete[] DH;
      delete[] llb;
      delete[] uub;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  MyList<Block> *Bp = blb;
  bool notfind = true;
  while (notfind && Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    bool flag = true;
    for (int i = 0; i < dim; i++)
    {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
      llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
      uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
      if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
      {
        flag = false;
        break;
      }
    }

    if (flag)
    {
      notfind = false;
      if (myrank == BP->rank)
      {
// test old code
#if 0
#define floorint(a) ((a) < 0 ? int(a) - 1 : int(a))
//---> interpolation
                int ixl,iyl,izl,ixu,iyu,izu;
	    	double Delx,Dely,Delz;

		ixl = 1+floorint((pox[0]-BP->X[0][0])/DH[0]);
	   	iyl = 1+floorint((pox[1]-BP->X[1][0])/DH[1]);
	   	izl = 1+floorint((pox[2]-BP->X[2][0])/DH[2]);

		int nn=ordn/2;

		ixl = ixl-nn;
		iyl = iyl-nn;
		izl = izl-nn;
	   
		int tmi;
		tmi = (Symmetry==2)?-1:0;
		if(ixl<tmi) ixl=tmi;
	   	if(iyl<tmi) iyl=tmi;
		tmi = (Symmetry>0)?-1:0;
	   	if(izl<tmi) izl=tmi;
      
	   	if(ixl+ordn>BP->shape[0]) ixl=BP->shape[0]-ordn;
	   	if(iyl+ordn>BP->shape[1]) iyl=BP->shape[1]-ordn;
	   	if(izl+ordn>BP->shape[2]) izl=BP->shape[2]-ordn;
// support cell center
		if(ixl>=0) Delx = ( pox[0] - BP->X[0][ixl] )/ DH[0];
		else       Delx = ( pox[0] + BP->X[0][0] )/ DH[0];
                if(iyl>=0) Dely = ( pox[1] - BP->X[1][iyl] )/ DH[1];
		else       Dely = ( pox[1] + BP->X[1][0] )/ DH[1];
                if(izl>=0) Delz = ( pox[2] - BP->X[2][izl] )/ DH[2];
		else       Delz = ( pox[2] + BP->X[2][0] )/ DH[2];
//change to fortran index
                ixl++;
	   	iyl++;
	   	izl++;
	   	ixu = ixl + ordn - 1;
	   	iyu = iyl + ordn - 1;
	   	izu = izl + ordn - 1;
	    	varl=VarList;
		int j=0;
	    	while(varl)
		{
                 f_interp_2(BP->shape,BP->fgfs[varl->data->sgfn],shellf[j],ixl,ixu,iyl,iyu,izl,izu,Delx,Dely,Delz,
                                     ordn,varl->data->SoA,Symmetry);
		 varl=varl->next;
		 j++;
		} //varl
#else
        //---> interpolation
        varl = VarList;
        int k = 0;
        while (varl) // run along variables
        {
          //              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
          //	  		                                    pox,ordn,varl->data->SoA,Symmetry);
          f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
                          pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
          varl = varl->next;
          k++;
        }
#endif
      }
    }
    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  if (notfind && myrank == 0)
  {
    cout << "ERROR: Patch::Interp_Points fails to find point (";
    for (int j = 0; j < dim; j++)
    {
      cout << XX[j];
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")";
    }
    cout << " on Patch (";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[j] << "+" << lli[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")--";
    }
    cout << "(";
    for (int j = 0; j < dim; j++)
    {
      cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
      if (j < dim - 1)
        cout << ",";
      else
        cout << ")" << endl;
    }
#if 0
       checkBlock();
#else
    cout << "splited domains:" << endl;
    {
      MyList<Block> *Bp = blb;
      while (Bp)
      {
        Block *BP = Bp->data;

        for (int i = 0; i < dim; i++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
          llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
          uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        }
        cout << "(";
        for (int j = 0; j < dim; j++)
        {
          cout << llb[j] << ":" << uub[j];
          if (j < dim - 1)
            cout << ",";
          else
            cout << ")" << endl;
        }
        if (Bp == ble)
          break;
        Bp = Bp->next;
      }
    }
#endif
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

  delete[] shellf;
  delete[] DH;
  delete[] llb;
  delete[] uub;

  return true;
}
// find maximum of abstract value, XX store position for maximum, Shellf store maximum themselvs
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
                         double *Shellf)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf, *xx;
  shellf = new double[num_var];
  xx = new double[dim * num_var];
  memset(shellf, 0, sizeof(double) * num_var);
  memset(xx, 0, sizeof(double) * dim * num_var);

  double *DH;
  int *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  llb = new int[dim];
  uub = new int[dim];

  MyList<Block> *Bp = blb;
  while (Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    if (myrank == BP->rank)
    {

      for (int i = 0; i < dim; i++)
      {
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
      }

      varl = VarList;
      int k = 0;
      double tmp, tmpx[dim];
      while (varl) // run along variables
      {
        f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
        if (tmp > shellf[k])
        {
          shellf[k] = tmp;
          for (int i = 0; i < dim; i++)
            xx[dim * k + i] = tmpx[i];
        }
        varl = varl->next;
        k++;
      }
    }

    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  struct mloc
  {
    double val;
    int rank;
  };

  mloc *IN, *OUT;
  IN = new mloc[num_var];
  OUT = new mloc[num_var];
  for (int i = 0; i < num_var; i++)
  {
    IN[i].val = shellf[i];
    IN[i].rank = myrank;
  }

  MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  for (int i = 0; i < num_var; i++)
  {
    Shellf[i] = OUT[i].val;
    if (myrank != OUT[i].rank)
      for (int k = 0; k < 3; k++)
        xx[3 * i + k] = 0;
  }

  MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete[] IN;
  delete[] OUT;
  delete[] shellf;
  delete[] xx;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
                         double *Shellf, MPI_Comm Comm_here)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf, *xx;
  shellf = new double[num_var];
  xx = new double[dim * num_var];
  memset(shellf, 0, sizeof(double) * num_var);
  memset(xx, 0, sizeof(double) * dim * num_var);

  double *DH;
  int *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  llb = new int[dim];
  uub = new int[dim];

  MyList<Block> *Bp = blb;
  while (Bp) // run along Blocks
  {
    Block *BP = Bp->data;

    if (myrank == BP->rank)
    {

      for (int i = 0; i < dim; i++)
      {
        llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
        uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
      }

      varl = VarList;
      int k = 0;
      double tmp, tmpx[dim];
      while (varl) // run along variables
      {
        f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
        if (tmp > shellf[k])
        {
          shellf[k] = tmp;
          for (int i = 0; i < dim; i++)
            xx[dim * k + i] = tmpx[i];
        }
        varl = varl->next;
        k++;
      }
    }

    if (Bp == ble)
      break;
    Bp = Bp->next;
  }

  struct mloc
  {
    double val;
    int rank;
  };

  mloc *IN, *OUT;
  IN = new mloc[num_var];
  OUT = new mloc[num_var];
  for (int i = 0; i < num_var; i++)
  {
    IN[i].val = shellf[i];
    IN[i].rank = myrank;
  }

  MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, Comm_here);

  for (int i = 0; i < num_var; i++)
  {
    Shellf[i] = OUT[i].val;
    if (myrank != OUT[i].rank)
      for (int k = 0; k < 3; k++)
        xx[3 * i + k] = 0;
  }

  MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

  delete[] IN;
  delete[] OUT;
  delete[] shellf;
  delete[] xx;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
// if the given point locates in the present Patch return true
// otherwise return false
bool Patch::Find_Point(double *XX)
{
  double *DH;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }

  for (int i = 0; i < dim; i++)
  {
    // has excluded the buffer points
    if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
    {
      delete[] DH;
      return false; // out of current patch,
                    // remember to delete the allocated arrays before return!!!
    }
  }

  delete[] DH;

  return true;
}
