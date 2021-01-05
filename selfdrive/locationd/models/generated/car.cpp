
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7305527196165684333) {
   out_7305527196165684333[0] = delta_x[0] + nom_x[0];
   out_7305527196165684333[1] = delta_x[1] + nom_x[1];
   out_7305527196165684333[2] = delta_x[2] + nom_x[2];
   out_7305527196165684333[3] = delta_x[3] + nom_x[3];
   out_7305527196165684333[4] = delta_x[4] + nom_x[4];
   out_7305527196165684333[5] = delta_x[5] + nom_x[5];
   out_7305527196165684333[6] = delta_x[6] + nom_x[6];
   out_7305527196165684333[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5181902846881929823) {
   out_5181902846881929823[0] = -nom_x[0] + true_x[0];
   out_5181902846881929823[1] = -nom_x[1] + true_x[1];
   out_5181902846881929823[2] = -nom_x[2] + true_x[2];
   out_5181902846881929823[3] = -nom_x[3] + true_x[3];
   out_5181902846881929823[4] = -nom_x[4] + true_x[4];
   out_5181902846881929823[5] = -nom_x[5] + true_x[5];
   out_5181902846881929823[6] = -nom_x[6] + true_x[6];
   out_5181902846881929823[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2569559215721594454) {
   out_2569559215721594454[0] = 1.0;
   out_2569559215721594454[1] = 0.0;
   out_2569559215721594454[2] = 0.0;
   out_2569559215721594454[3] = 0.0;
   out_2569559215721594454[4] = 0.0;
   out_2569559215721594454[5] = 0.0;
   out_2569559215721594454[6] = 0.0;
   out_2569559215721594454[7] = 0.0;
   out_2569559215721594454[8] = 0.0;
   out_2569559215721594454[9] = 1.0;
   out_2569559215721594454[10] = 0.0;
   out_2569559215721594454[11] = 0.0;
   out_2569559215721594454[12] = 0.0;
   out_2569559215721594454[13] = 0.0;
   out_2569559215721594454[14] = 0.0;
   out_2569559215721594454[15] = 0.0;
   out_2569559215721594454[16] = 0.0;
   out_2569559215721594454[17] = 0.0;
   out_2569559215721594454[18] = 1.0;
   out_2569559215721594454[19] = 0.0;
   out_2569559215721594454[20] = 0.0;
   out_2569559215721594454[21] = 0.0;
   out_2569559215721594454[22] = 0.0;
   out_2569559215721594454[23] = 0.0;
   out_2569559215721594454[24] = 0.0;
   out_2569559215721594454[25] = 0.0;
   out_2569559215721594454[26] = 0.0;
   out_2569559215721594454[27] = 1.0;
   out_2569559215721594454[28] = 0.0;
   out_2569559215721594454[29] = 0.0;
   out_2569559215721594454[30] = 0.0;
   out_2569559215721594454[31] = 0.0;
   out_2569559215721594454[32] = 0.0;
   out_2569559215721594454[33] = 0.0;
   out_2569559215721594454[34] = 0.0;
   out_2569559215721594454[35] = 0.0;
   out_2569559215721594454[36] = 1.0;
   out_2569559215721594454[37] = 0.0;
   out_2569559215721594454[38] = 0.0;
   out_2569559215721594454[39] = 0.0;
   out_2569559215721594454[40] = 0.0;
   out_2569559215721594454[41] = 0.0;
   out_2569559215721594454[42] = 0.0;
   out_2569559215721594454[43] = 0.0;
   out_2569559215721594454[44] = 0.0;
   out_2569559215721594454[45] = 1.0;
   out_2569559215721594454[46] = 0.0;
   out_2569559215721594454[47] = 0.0;
   out_2569559215721594454[48] = 0.0;
   out_2569559215721594454[49] = 0.0;
   out_2569559215721594454[50] = 0.0;
   out_2569559215721594454[51] = 0.0;
   out_2569559215721594454[52] = 0.0;
   out_2569559215721594454[53] = 0.0;
   out_2569559215721594454[54] = 1.0;
   out_2569559215721594454[55] = 0.0;
   out_2569559215721594454[56] = 0.0;
   out_2569559215721594454[57] = 0.0;
   out_2569559215721594454[58] = 0.0;
   out_2569559215721594454[59] = 0.0;
   out_2569559215721594454[60] = 0.0;
   out_2569559215721594454[61] = 0.0;
   out_2569559215721594454[62] = 0.0;
   out_2569559215721594454[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8098244067719744350) {
   out_8098244067719744350[0] = state[0];
   out_8098244067719744350[1] = state[1];
   out_8098244067719744350[2] = state[2];
   out_8098244067719744350[3] = state[3];
   out_8098244067719744350[4] = state[4];
   out_8098244067719744350[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8098244067719744350[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8098244067719744350[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8290869482667701516) {
   out_8290869482667701516[0] = 1;
   out_8290869482667701516[1] = 0;
   out_8290869482667701516[2] = 0;
   out_8290869482667701516[3] = 0;
   out_8290869482667701516[4] = 0;
   out_8290869482667701516[5] = 0;
   out_8290869482667701516[6] = 0;
   out_8290869482667701516[7] = 0;
   out_8290869482667701516[8] = 0;
   out_8290869482667701516[9] = 1;
   out_8290869482667701516[10] = 0;
   out_8290869482667701516[11] = 0;
   out_8290869482667701516[12] = 0;
   out_8290869482667701516[13] = 0;
   out_8290869482667701516[14] = 0;
   out_8290869482667701516[15] = 0;
   out_8290869482667701516[16] = 0;
   out_8290869482667701516[17] = 0;
   out_8290869482667701516[18] = 1;
   out_8290869482667701516[19] = 0;
   out_8290869482667701516[20] = 0;
   out_8290869482667701516[21] = 0;
   out_8290869482667701516[22] = 0;
   out_8290869482667701516[23] = 0;
   out_8290869482667701516[24] = 0;
   out_8290869482667701516[25] = 0;
   out_8290869482667701516[26] = 0;
   out_8290869482667701516[27] = 1;
   out_8290869482667701516[28] = 0;
   out_8290869482667701516[29] = 0;
   out_8290869482667701516[30] = 0;
   out_8290869482667701516[31] = 0;
   out_8290869482667701516[32] = 0;
   out_8290869482667701516[33] = 0;
   out_8290869482667701516[34] = 0;
   out_8290869482667701516[35] = 0;
   out_8290869482667701516[36] = 1;
   out_8290869482667701516[37] = 0;
   out_8290869482667701516[38] = 0;
   out_8290869482667701516[39] = 0;
   out_8290869482667701516[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8290869482667701516[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8290869482667701516[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8290869482667701516[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8290869482667701516[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8290869482667701516[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8290869482667701516[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8290869482667701516[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8290869482667701516[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8290869482667701516[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8290869482667701516[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8290869482667701516[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8290869482667701516[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8290869482667701516[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8290869482667701516[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8290869482667701516[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8290869482667701516[56] = 0;
   out_8290869482667701516[57] = 0;
   out_8290869482667701516[58] = 0;
   out_8290869482667701516[59] = 0;
   out_8290869482667701516[60] = 0;
   out_8290869482667701516[61] = 0;
   out_8290869482667701516[62] = 0;
   out_8290869482667701516[63] = 1;
}
void h_25(double *state, double *unused, double *out_86988055339073109) {
   out_86988055339073109[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6071874413648555756) {
   out_6071874413648555756[0] = 0;
   out_6071874413648555756[1] = 0;
   out_6071874413648555756[2] = 0;
   out_6071874413648555756[3] = 0;
   out_6071874413648555756[4] = 0;
   out_6071874413648555756[5] = 0;
   out_6071874413648555756[6] = 1;
   out_6071874413648555756[7] = 0;
}
void h_24(double *state, double *unused, double *out_1524082693182885610) {
   out_1524082693182885610[0] = state[4];
   out_1524082693182885610[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1603688280036473186) {
   out_1603688280036473186[0] = 0;
   out_1603688280036473186[1] = 0;
   out_1603688280036473186[2] = 0;
   out_1603688280036473186[3] = 0;
   out_1603688280036473186[4] = 1;
   out_1603688280036473186[5] = 0;
   out_1603688280036473186[6] = 0;
   out_1603688280036473186[7] = 0;
   out_1603688280036473186[8] = 0;
   out_1603688280036473186[9] = 0;
   out_1603688280036473186[10] = 0;
   out_1603688280036473186[11] = 0;
   out_1603688280036473186[12] = 0;
   out_1603688280036473186[13] = 1;
   out_1603688280036473186[14] = 0;
   out_1603688280036473186[15] = 0;
}
void h_30(double *state, double *unused, double *out_6857823281689424045) {
   out_6857823281689424045[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3542024497300627524) {
   out_3542024497300627524[0] = 0;
   out_3542024497300627524[1] = 0;
   out_3542024497300627524[2] = 0;
   out_3542024497300627524[3] = 0;
   out_3542024497300627524[4] = 1;
   out_3542024497300627524[5] = 0;
   out_3542024497300627524[6] = 0;
   out_3542024497300627524[7] = 0;
}
void h_26(double *state, double *unused, double *out_8881011243771440403) {
   out_8881011243771440403[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7378110185603396901) {
   out_7378110185603396901[0] = 0;
   out_7378110185603396901[1] = 0;
   out_7378110185603396901[2] = 0;
   out_7378110185603396901[3] = 0;
   out_7378110185603396901[4] = 0;
   out_7378110185603396901[5] = 0;
   out_7378110185603396901[6] = 0;
   out_7378110185603396901[7] = 1;
}
void h_27(double *state, double *unused, double *out_1254731856415214846) {
   out_1254731856415214846[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4829606485137252836) {
   out_4829606485137252836[0] = 0;
   out_4829606485137252836[1] = 0;
   out_4829606485137252836[2] = 0;
   out_4829606485137252836[3] = 1;
   out_4829606485137252836[4] = 0;
   out_4829606485137252836[5] = 0;
   out_4829606485137252836[6] = 0;
   out_4829606485137252836[7] = 0;
}
void h_29(double *state, double *unused, double *out_4153777357591426758) {
   out_4153777357591426758[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7365823266429459477) {
   out_7365823266429459477[0] = 0;
   out_7365823266429459477[1] = 1;
   out_7365823266429459477[2] = 0;
   out_7365823266429459477[3] = 0;
   out_7365823266429459477[4] = 0;
   out_7365823266429459477[5] = 0;
   out_7365823266429459477[6] = 0;
   out_7365823266429459477[7] = 0;
}
void h_28(double *state, double *unused, double *out_7232633107963649866) {
   out_7232633107963649866[0] = state[5];
   out_7232633107963649866[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8621520758272968052) {
   out_8621520758272968052[0] = 0;
   out_8621520758272968052[1] = 0;
   out_8621520758272968052[2] = 0;
   out_8621520758272968052[3] = 0;
   out_8621520758272968052[4] = 0;
   out_8621520758272968052[5] = 1;
   out_8621520758272968052[6] = 0;
   out_8621520758272968052[7] = 0;
   out_8621520758272968052[8] = 0;
   out_8621520758272968052[9] = 0;
   out_8621520758272968052[10] = 0;
   out_8621520758272968052[11] = 0;
   out_8621520758272968052[12] = 0;
   out_8621520758272968052[13] = 0;
   out_8621520758272968052[14] = 1;
   out_8621520758272968052[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
