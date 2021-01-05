/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4464007461933898312);
void inv_err_fun(double *nom_x, double *true_x, double *out_4014484379048568467);
void H_mod_fun(double *state, double *out_3971040135259682870);
void f_fun(double *state, double dt, double *out_8955640682202544921);
void F_fun(double *state, double dt, double *out_2770007890841314472);
void h_3(double *state, double *unused, double *out_562273081113642881);
void H_3(double *state, double *unused, double *out_6250517764500677237);
void h_4(double *state, double *unused, double *out_7300734182737055273);
void H_4(double *state, double *unused, double *out_5853402369882574238);
void h_9(double *state, double *unused, double *out_5375458734568416845);
void H_9(double *state, double *unused, double *out_3515132462623138557);
void h_10(double *state, double *unused, double *out_4791808667020168634);
void H_10(double *state, double *unused, double *out_7211772746493218348);
void h_12(double *state, double *unused, double *out_7450697125728691139);
void H_12(double *state, double *unused, double *out_3746846858143364669);
void h_31(double *state, double *unused, double *out_7224047753869266240);
void H_31(double *state, double *unused, double *out_6421936633286784490);
void h_32(double *state, double *unused, double *out_6110400123915083194);
void H_32(double *state, double *unused, double *out_4344055274088874277);
void h_13(double *state, double *unused, double *out_2280115523582535088);
void H_13(double *state, double *unused, double *out_412605668220139245);
void h_14(double *state, double *unused, double *out_5375458734568416845);
void H_14(double *state, double *unused, double *out_3515132462623138557);
void h_19(double *state, double *unused, double *out_8367225692275273032);
void H_19(double *state, double *unused, double *out_7947812313725569738);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);