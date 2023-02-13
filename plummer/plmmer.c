#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#define N 1000 //粒子数
#define R 1.0  //モデル半径
#define G 1.0
#define M 1.0
#define m 1.0 / (double)N

double t = 0.0;
double x;
double y;
double z;
double r_0[N];
double V_e[N];
double r[N][3];
double v[N][3];
double a[N][3];
double Ep[N];
double Ek[N];
double du_dt(double x_i, double x_j, double r_ij_3);
double dv_dt(double y_i, double y_j, double r_ij_3);
double dw_dt(double z_i, double z_j, double r_ij_3);
double Ep_total = 0.0;
double Ek_total = 0.0;
double k;
double eta = 0.1;
double dt;
double t_ff;          //自由落下時間
double eps = R / 100; //ソフトニングパラメータ
double c_virial;      //ビリアル比

double rho_0;
// double V;

int main()
{
  FILE *fp;
  char *fname = "plummer.csv";
  fp = fopen(fname, "w");
  int i = 0;
  int j = 0;
  int k = 0;
  double eps_2 = eps * eps;
  // 座標の初期化
  for (i = 0; i < N; i++)
  {
    // 0~1の乱数作成
    double X1 = (double)rand() / (RAND_MAX + 1.0);
    double X2 = (double)rand() / (RAND_MAX + 1.0);
    double X3 = (double)rand() / (RAND_MAX + 1.0);
    r_0[i] = 1.0 / sqrt(pow(X1, -2.0 / 3.0) - 1.0);
    double zn = (1.0 - 2.0 * X2) * r_0[i];
    double xn = sqrt(r_0[i] * r_0[i] - zn * zn) * cos(2.0 * M_PI * X3);
    double yn = sqrt(r_0[i] * r_0[i] - zn * zn) * sin(2.0 * M_PI * X3);

    r[i][0] = xn;
    r[i][1] = yn;
    r[i][2] = zn;
    // printf("%f\n",r_0[i]);
  }

  //ポテンシャルエネルギーを求める。
  for (i = 0; i < N; i++)
  {
    Ep[i] = 0.0;
  }
  for (i = 0; i < N; i++)
  {
    Ep[i] = (G * M * m) / eps;
    double x_i = r[i][0];
    double y_i = r[i][1];
    double z_i = r[i][2];
    for (j = 0; j < N; j++)
    {

      Ep[i] += -(G * M * m) / sqrt((x_i - r[j][0]) * (x_i - r[j][0]) + (y_i - r[j][1]) * (y_i - r[j][1]) + (z_i - r[j][2]) * (z_i - r[j][2]) + eps_2);
    }
    Ep_total += Ep[i];
  }
  Ep_total = 0.5 * Ep_total;
  //　ビリアル速度の変化を見る(ビリアル速度の係数を変化させる)　
  // double v_virial= sqrt(fabs(Ep_total) /(2*N));//ビリアル速度
  // double v_virial_dash = 0.7*v_virial;

  // 速度の初期化
  i = 0;
  while (i < N)
  {
    // for (i=0;i<N;i++){
    V_e[i] = sqrt(2.0) / pow(1.0 + r_0[i] * r_0[i], 0.25);
    double X4 = (double)rand() / (RAND_MAX + 1.0);
    double X5 = (double)rand() / (RAND_MAX + 1.0);
    if (0.1 * X5 < X4 * X4 * pow(1.0 - X4 * X4, 7.0 / 2.0))
    {
      // printf("%f\n",X4);
      double V = X4 * V_e[i];
      double X6 = (double)rand() / (RAND_MAX + 1.0);
      double X7 = (double)rand() / (RAND_MAX + 1.0);
      double wn = (1.0 - 2.0 * X6) * V;
      double un = sqrt(V * V - wn * wn) * cos(2.0 * M_PI * X7);
      double vn = sqrt(V * V - wn * wn) * sin(2.0 * M_PI * X7);

      v[i][0] = un;
      v[i][1] = vn;
      v[i][2] = wn;

      k++;
      if (k == 100000)
      {
        printf("１0万回");
        break;
      }

      // printf("%f\n",0.1*X5 - X4*X4*pow(1.0-X4*X4,7.0/2.0) );
      // printf("%f\n",v[i][0]);
      i++;
    }
  }

  // for(i=0;i<N;i++){
  //   printf("%f,\n",v[i][0]);
  // }

  for (i = 0; i < N; i++)
  {
    Ek[i] = (0.5) * (pow(v[i][0], 2.0) + pow(v[i][1], 2.0) + pow(v[i][2], 2.0));
    Ek_total += Ek[i];
  }
  // printf("%f\n",Ek_total);

  // 加速度の初期値
  for (i = 0; i < N; i++)
  {
    a[i][0] = 0.0;
    a[i][1] = 0.0;
    a[i][2] = 0.0;
  }

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N - 1; j++)
    {

      double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps_2, 1.5);
      a[i][0] += du_dt(r[i][0], r[j][0], r_ij_3);
      a[i][1] += dv_dt(r[i][1], r[j][1], r_ij_3);
      a[i][2] += dw_dt(r[i][2], r[j][2], r_ij_3);
    }
  }

  // leap frog
  rho_0 = 3.0 * M / (4.0 * M_PI * R * R * R * pow(2.0, 2.5));
  t_ff = sqrt(3.0 * M_PI / (32.0 * G * rho_0)); //自由落下時間
  double v_typical = R / t_ff;
  double eps_pl = R / 10.0;
  dt = eps_pl / v_typical;
  // printf("%f,%f\n",t_ff,dt);
  for (t = 0; t < 10 * t_ff; t += dt)
  {
    fprintf(fp, "%f\n", t);
    for (i = 0; i < N; i++)
    {

      v[i][0] += (dt * a[i][0]) / 2.0;
      v[i][1] += (dt * a[i][1]) / 2.0;
      v[i][2] += (dt * a[i][2]) / 2.0;
      r[i][0] += dt * v[i][0];
      r[i][1] += dt * v[i][1];
      r[i][2] += dt * v[i][2];
    }
    for (i = 0; i < N; i++)
    {
      a[i][0] = 0.0;
      a[i][1] = 0.0;
      a[i][2] = 0.0;
    }
    for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
      {
        double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps_2, 1.5);
        a[i][0] += du_dt(r[i][0], r[j][0], r_ij_3);
        a[i][1] += dv_dt(r[i][1], r[j][1], r_ij_3);
        a[i][2] += dw_dt(r[i][2], r[j][2], r_ij_3);
      }

      v[i][0] += (dt * a[i][0]) / 2.0;
      v[i][1] += (dt * a[i][1]) / 2.0;
      v[i][2] += (dt * a[i][2]) / 2.0;
    }

    // printf("%f\n",v[1][0]);

    Ep_total = 0;
    for (i = 0; i < N; i++)
    {
      Ep[i] = (G * M * m) / 1.0;
      for (j = 0; j < N; j++)
      {
        double r_ij = sqrt((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps_2);
        double r_ij_2 = sqrt((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + 1);
        Ep[i] += -(G * M * m) / r_ij;
      }
      Ep_total += Ep[i];
    }
    Ep_total = 0.5 * Ep_total;

    Ek_total = 0;
    for (i = 0; i < N; i++)
    {
      Ek[i] = (0.5) * (pow(v[i][0], 2.0) + pow(v[i][1], 2.0) + pow(v[i][2], 2.0));
      Ek_total += Ek[i];
    }
    double c_virial = (2 * Ek_total) / fabs(Ep_total);
    double E_total = Ep_total + Ek_total;

    // fprintf(fp,"%f,%f,%f,%f,%f\n",t,Ep_total,Ek_total,E_total,c_virial);
    for (i = 0; i < N; i++)
    {
      fprintf(fp, "%f,%f,%f\n", r[i][0], r[i][1], r[i][2]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "\n");
  }
  fclose(fp);
}

double du_dt(double x_i, double x_j, double r_ij_3)
{
  return -m * (x_i - x_j) / r_ij_3;
}
double dv_dt(double y_i, double y_j, double r_ij_3)
{
  return -m * (y_i - y_j) / r_ij_3;
}
double dw_dt(double z_i, double z_j, double r_ij_3)
{
  return -m * (z_i - z_j) / r_ij_3;
}
