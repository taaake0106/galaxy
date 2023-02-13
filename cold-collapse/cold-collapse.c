#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#define N 1000 // 粒子数
#define R 1.0  // モデルの半径
#define G 1
#define M 1
#define m 1
int test = 0;
double t = 0.0;
double x;
double y;
double z;
double r[N][3];
double v[N][3];
double a[N][3];
double Ep[N];
double Ek[N];
double du_dt(double x_i, double x_j, double r_ij_3);
double dv_dt(double y_i, double y_j, double r_ij_3);
double dw_dt(double z_i, double z_j, double r_ij_3);
double Ep_total = 0;
double Ek_total = 0;
double Em_total = 0;
double t_unit = 14.9181;
double eta = 0.1;
double dt;
double t_ff;
int last_chunk = 0;
double eps = R / 100; // ソフトニングパラメータ

int main()
{
    FILE *fp;
    char *fname = "cold-collapse.dat";
    fp = fopen(fname, "w");
    FILE *fp_enegy;
    char *fname_enegy = "cold-collapse_enegy.dat";
    fp_enegy = fopen(fname_enegy, "w");
    int i = 0;
    int j = 0;
    int k = 0;
    double eps2 = eps * eps;
    // 座標の初期化
    while (i < N)
    {
        // 0~1の乱数作成
        double xn = (double)rand() / RAND_MAX;
        double yn = (double)rand() / RAND_MAX;
        double zn = (double)rand() / RAND_MAX;

        double x = 2.0 * xn - 1.0;
        double y = 2.0 * yn - 1.0;
        double z = 2.0 * zn - 1.0;

        if ((double)R >= pow(x * x + y * y + z * z, 0.5))
        {
            double test = pow(x * x + y * y + z * z, 0.5);
            r[i][0] = x;
            r[i][1] = y;
            r[i][2] = z;
            i++;
        }
        // ストッパー
        k++;
        if (k == 10000)
        {
            printf("10000回に達しました。");
            break;
        }
    }

    // 速度の初期化
    for (i = 0; i < N; i++)
    {
        Ep[i] = 0;
    }

    // Ep_totalからv_virialを求めにいく
    for (i = 0; i < N; i++)
    {
        Ep[i] = (G * m * m) / eps;
        double x_i = r[i][0];
        double y_i = r[i][1];
        double z_i = r[i][2];
        for (j = 0; j < N; j++)
        {
            Ep[i] += -(G * m * m) / sqrt((x_i - r[j][0]) * (x_i - r[j][0]) + (y_i - r[j][1]) * (y_i - r[j][1]) + (z_i - r[j][2]) * (z_i - r[j][2]) + eps2);
        }
        Ep_total += Ep[i];
    }
    Ep_total = Ep_total * 0.5;

    // v_virialの変化も見てみる  物理状態を見るために物理パラメータを変化させる。
    double v_virial = pow(fabs(Ep_total) / (2 * N), 0.5);
    double v_virial_dash = v_virial;

    // 速度を初期化していく
    for (i = 0; i < N; i++)
    {
        double un = (double)rand() / RAND_MAX;
        double vn = (double)rand() / RAND_MAX;
        double wn = (double)rand() / RAND_MAX;
        double A = 2.0 * un - 1.0;
        double B = 2.0 * vn - 1.0;
        double C = 2.0 * wn - 1.0;

        v[i][0] = v_virial_dash * (A / pow(A * A + B * B + C * C, 0.5));
        v[i][1] = v_virial_dash * (B / pow(A * A + B * B + C * C, 0.5));
        v[i][2] = v_virial_dash * (C / pow(A * A + B * B + C * C, 0.5));
    }
    for (i = 0; i < N; i++)
    {

        Ek[i] = (0.5) * (pow(v[i][0], 2.0) + pow(v[i][1], 2.0) + pow(v[i][2], 2.0));
        Ek_total += Ek[i];
    }

    // 加速度の初期値
    for (i = 0; i < N; i++)
    {
        a[i][0] = 0.0;
        a[i][1] = 0.0;
        a[i][2] = 0.0;
    }

    // 加速度の初期化
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps, 1.5);
            a[i][0] += du_dt(r[i][0], r[j][0], r_ij_3);
            a[i][1] += dv_dt(r[i][1], r[j][1], r_ij_3);
            a[i][2] += dw_dt(r[i][2], r[j][2], r_ij_3);
        }
    }

    // 以下リープフロッグ
    dt = (eps / v_virial) * eta;
    t_ff = sqrt((R * R * R * M_PI * M_PI) / (8 * N));
    int step = 0;

    for (t = 0; t < 10 * t_ff; t += dt)
    {
        // まず、vを0.5dtでvを作り、そのvでrを作る
        for (i = 0; i < N; i++)
        {
            for (k = 0; k < 3; k++)
            {
                v[i][k] += (dt * a[i][k]) / 2.0;
                r[i][k] += dt * v[i][k];
            }
        }
        // 次に、aを求め直し、そのaでvをつくる。
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
                double r_ij_3 = pow((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps2, 1.5);
                a[i][0] += du_dt(r[i][0], r[j][0], r_ij_3);
                a[i][1] += dv_dt(r[i][1], r[j][1], r_ij_3);
                a[i][2] += dw_dt(r[i][2], r[j][2], r_ij_3);
            }
            for (k = 0; k < 3; k++)
            {
                v[i][k] += (dt * a[i][k]) / 2.0;
            }
        }
        if (step % 50 == 0)
        {
            // 位置と速度座標への書きこみ
            fprintf(fp, "%f\n", t * t_unit);
            for (i = 0; i < N; i++)
            {
                fprintf(fp, "%f,%f,%f,%f,%f,%f\n", r[i][0], r[i][1], r[i][2], v[i][0], v[i][1], v[i][2]);
            }
            fprintf(fp, "\n\n");

            // エネルギーの書き込み

            // 初期化
            Ep_total = 0;
            Ek_total = 0;
            for (i = 0; i < N; i++)
            {
                Ep[i] = 0;
                Ek[i] = 0;
            }

            // Ep_totalからv_virialを求めにいく
            for (i = 0; i < N; i++)
            {
                Ep[i] = (G * m * m) / eps;
                double x_i = r[i][0];
                double y_i = r[i][1];
                double z_i = r[i][2];
                for (j = 0; j < N; j++)
                {
                    Ep[i] += -(G * m * m) / sqrt((x_i - r[j][0]) * (x_i - r[j][0]) + (y_i - r[j][1]) * (y_i - r[j][1]) + (z_i - r[j][2]) * (z_i - r[j][2]) + eps2);
                }
                Ep_total += Ep[i];
            }
            Ep_total = Ep_total * 0.5;

            for (i = 0; i < N; i++)
            {
                Ek[i] = (0.5) * (pow(v[i][0], 2.0) + pow(v[i][1], 2.0) + pow(v[i][2], 2.0));
                Ek_total += Ek[i];
            }

            double E_unit = 4.966 * pow(10.0, 6.0);

            fprintf(fp_enegy, "%f,%f,%f,%f,%f\n", t * t_unit, Ek_total * E_unit, Ep_total * E_unit, (Ek_total + Ep_total) * E_unit, -Ek_total / Ep_total);
        }
        last_chunk++;
        step++;
    } // リープフロッグ終わり
    fclose(fp);
    fclose(fp_enegy);
    printf("%d", last_chunk);
    printf("v_virial=%f[km s^-1],t_ff=%f[Myr],dt=%f[Myr]", v_virial * 65.589, t_ff * t_unit, dt * t_unit);
}

double du_dt(double x_i, double x_j, double r_ij_3)
{
    return -(x_i - x_j) / r_ij_3;
}
double dv_dt(double y_i, double y_j, double r_ij_3)
{
    return -(y_i - y_j) / r_ij_3;
}
double dw_dt(double z_i, double z_j, double r_ij_3)
{
    return -(z_i - z_j) / r_ij_3;
}
