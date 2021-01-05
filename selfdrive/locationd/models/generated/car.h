/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7305527196165684333);
void inv_err_fun(double *nom_x, double *true_x, double *out_5181902846881929823);
void H_mod_fun(double *state, double *out_2569559215721594454);
void f_fun(double *state, double dt, double *out_8098244067719744350);
void F_fun(double *state, double dt, double *out_8290869482667701516);
void h_25(double *state, double *unused, double *out_86988055339073109);
void H_25(double *state, double *unused, double *out_6071874413648555756);
void h_24(double *state, double *unused, double *out_1524082693182885610);
void H_24(double *state, double *unused, double *out_1603688280036473186);
void h_30(double *state, double *unused, double *out_6857823281689424045);
void H_30(double *state, double *unused, double *out_3542024497300627524);
void h_26(double *state, double *unused, double *out_8881011243771440403);
void H_26(double *state, double *unused, double *out_7378110185603396901);
void h_27(double *state, double *unused, double *out_1254731856415214846);
void H_27(double *state, double *unused, double *out_4829606485137252836);
void h_29(double *state, double *unused, double *out_4153777357591426758);
void H_29(double *state, double *unused, double *out_7365823266429459477);
void h_28(double *state, double *unused, double *out_7232633107963649866);
void H_28(double *state, double *unused, double *out_8621520758272968052);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
