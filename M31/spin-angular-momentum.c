// BHを重心と仮定し、重心を中心とした密度分布を考える。

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define N 550000           // 粒子数
#define M_sat 1            // 衛星銀河の質量(×10^9 M_sun)
#define G 1.0              // 計算機内の重力定数G’
#define m 0.00001          // 衛星銀河を構成する粒子の質量(10^3 M_sun)
#define m_BH 0.001         // BHの質量(10^8 M_sun)
#define r_c 1.0            // 衛星銀河のコア半径(kpc)
#define gc_margin 0.000001 // 重心座標の相相対誤差
#define eps 0.1            // ソフトニングパラメーター

// stellarのデータ
#define I (550000) // 行の数
#define J (7)     // 列の数

#define N_bulge 100000 // bulge粒子の数
#define N_disk 450000  // disk粒子の数

// BHのデータ
#define K (200) // 行の数 (stellarのデータの数が100まで)
#define L (8)   // 列の数

double r_bulge[N][3]; // 各bulge粒子の座標
double v_bulge[N][3]; // 各bulge粒子の速度

double r_disk[N][3]; // 各disk粒子の座標
double v_disk[N][3]; // 各disk粒子の速度

double x_rela_BH_bulge[N]; // BHから見たbulge粒子の相対座標
double y_rela_BH_bulge[N];
double z_rela_BH_bulge[N];

double v_rela_BH_x_bulge[N]; // BHから見たbulge粒子の相対速度
double v_rela_BH_y_bulge[N];
double v_rela_BH_z_bulge[N];

double x_rela_BH_disk[N]; // BHから見たdisk粒子の相対座標
double y_rela_BH_disk[N];
double z_rela_BH_disk[N];

double v_rela_BH_x_disk[N]; // BHから見たdisk粒子の相対速度
double v_rela_BH_y_disk[N];
double v_rela_BH_z_disk[N];

double L_bulge[100000][3]; //各バルジ粒子の角運動量
double L_disk[450000][3];//各ディスク粒子の角運動量

double L_bulge_tot[3];//バルジ粒子の角運動量
double L_disk_tot[3];//各ディスク粒子の角運動量
double L_tot[3];

// stellarのファイルからデータを読み込む
void read_bound(char filename1[], double mat_bound[I][J])
{
    int i, j;
    FILE *fpstellar;

    fpstellar = fopen(filename1, "r"); /*  読み込みモードでファイルをオープン  */
    if (fpstellar == NULL)
    {
        printf("ファイルを開くことが出来ませんでした．\n");
    }

    for (i = 0; i < I; i++)
    {

        fscanf(fpstellar, "%lf %lf %lf %lf %lf %lf %lf ", &mat_bound[i][0], &mat_bound[i][1], &mat_bound[i][2], &mat_bound[i][3], &mat_bound[i][4], &mat_bound[i][5], &mat_bound[i][6]); /*  1行読む  */
    }

    fclose(fpstellar);
}

// BHのファイルからデータを読み込む 1粒子なので１行に一つの時間
void read_BH(char filename2[], double mat_BH[L])
{
    int i, j;
    FILE *fpBH;

    fpBH = fopen(filename2, "r"); /*  読み込みモードでファイルをオープン  */
    if (fpBH == NULL)
    {
        printf("ファイルを開くことが出来ませんでした．\n");
    }

    fscanf(fpBH, "%lf %lf %lf %lf %lf %lf %lf %lf", &mat_BH[0], &mat_BH[1], &mat_BH[2], &mat_BH[3], &mat_BH[4], &mat_BH[5], &mat_BH[6], &mat_BH[7]); /*  1行読む  */

    fclose(fpBH);
}

void read_num_bound(char filename2[], double mat_bound_particle[2])
{
    int i, j;
    FILE *fpnum;

    fpnum = fopen(filename2, "r"); /*  読み込みモードでファイルをオープン  */
    if (fpnum == NULL)
    {
        printf("ファイルを開くことが出来ませんでした．\n");
    }

    fscanf(fpnum, "%lf %lf ", &mat_bound_particle[0], &mat_bound_particle[1]); /*  1行読む  */

    fclose(fpnum);
}

int main()
{

    // PI=M_PI;
    char f_bound[256];
    char f_BH[256];
    char f_num[256];
    int number;
    int i;
    double data_bound[I][J];
    double data_BH[L];
    double data_num[2];
    FILE *fpwrite1;
    FILE *fpwrite2;
    FILE *fpwrite3;
    char *fname1[256];
    char *fname2[256];
    char *fname3[256];
    // FILE *fpwrite;
    // char *fname = "num_bound_particle_margin=100000.dat";
    // fpwrite = fopen(fname, "w");

    for (number = 0; number < 101; number++)
    {
        if (number % 10 == 0 || number == 94)
        {

            sprintf(fname1, "L_3t=%d.dat", number);
            fpwrite1 = fopen(("%s", fname1), "w");
            // sprintf(fname2, "analy_bound_ptcl_t=%05d.dat", number);
            // fpwrite2 = fopen(("%s", fname2), "w");
            // sprintf(fname3, "analy_unbound_ptcl_t=%05d.dat", number);
            // fpwrite3 = fopen(("%s", fname3), "w");
            // sprintf(filename1,"satellitedata/stellar/data_rini%05d.txt",number);
            sprintf(f_bound, "../stellar.BH.DMH/5analy_bound_ptcl_DMH_t=%05d.dat", number);
            sprintf(f_BH, "../stellarandBH/stellar/BHpos-vel%05d.txt", number);
            sprintf(f_num, "../stellar.BH.DMH/5analy_num_bound_ptcl_DMH_t=%05d.dat", number);
            read_bound(f_bound, data_bound);
            read_BH(f_BH, data_BH);
            read_num_bound(f_num, data_num);
            printf("%f\n", data_num[1]); // N->data_numとなる

            for (i = 0; i < data_num[1]; i++)
            {
                if (data_bound[i][0] < N_bulge)
                {

                    v_bulge[i][0] = 0.0;
                    v_bulge[i][1] = 0.0;
                    v_bulge[i][2] = 0.0;
                }

                if (data_bound[i][0] > N_bulge - 1)
                {

                    v_disk[i][0] = 0.0;
                    v_disk[i][1] = 0.0;
                    v_disk[i][2] = 0.0;
                }
            }

            // 粒子の座標を代入
            for (i = 0; i < data_num[1]; i++)
            {
                if (data_bound[i][0] < N_bulge)
                {
                    r_bulge[i][0] = data_bound[i][1];
                    r_bulge[i][1] = data_bound[i][2];
                    r_bulge[i][2] = data_bound[i][3];
                }

                if (data_bound[i][0] > N_bulge - 1)
                {
                    r_disk[i][0] = data_bound[i][1];
                    r_disk[i][1] = data_bound[i][2];
                    r_disk[i][2] = data_bound[i][3];
                }
            }

            for (i = 0; i < data_num[1]; i++)
            {
                if (data_bound[i][0] < N_bulge)
                {

                    v_bulge[i][0] = 0.0;
                    v_bulge[i][1] = 0.0;
                    v_bulge[i][2] = 0.0;
                }

                if (data_bound[i][0] > N_bulge - 1)
                {

                    v_disk[i][0] = 0.0;
                    v_disk[i][1] = 0.0;
                    v_disk[i][2] = 0.0;
                }
            }

            // stellar粒子の速度を代入
            for (i = 0; i < data_num[1]; i++)
            {
                if (data_bound[i][0] < N_bulge)
                {
                    v_bulge[i][0] = data_bound[i][4];
                    v_bulge[i][1] = data_bound[i][5];
                    v_bulge[i][2] = data_bound[i][6];
                }

                if (data_bound[i][0] > N_bulge - 1)
                {
                    v_disk[i][0] = data_bound[i][4];
                    v_disk[i][1] = data_bound[i][5];
                    v_disk[i][2] = data_bound[i][6];
                }
            }

            double BH_pos_x;
            double BH_pos_y;
            double BH_pos_z;
            BH_pos_x = 0.0;
            BH_pos_y = 0.0;
            BH_pos_z = 0.0;

            // BHの座標を代入
            BH_pos_x = data_BH[2];
            BH_pos_y = data_BH[3];
            BH_pos_z = data_BH[4];

            double BH_vel_x;
            double BH_vel_y;
            double BH_vel_z;
            BH_vel_x = 0.0;
            BH_vel_y = 0.0;
            BH_vel_z = 0.0;

            // BHの速度を代入
            BH_vel_x = data_BH[5];
            BH_vel_y = data_BH[6];
            BH_vel_z = data_BH[7];

            for (i = 0; i < data_num[1]; i++)
            {
                if (data_bound[i][0] < N_bulge)
                {
                    x_rela_BH_bulge[i] = 0.0;
                    y_rela_BH_bulge[i] = 0.0;
                    z_rela_BH_bulge[i] = 0.0;

                    v_rela_BH_x_bulge[i] = 0.0;
                    v_rela_BH_y_bulge[i] = 0.0;
                    v_rela_BH_z_bulge[i] = 0.0;
                }

                if (data_bound[i][0] < N_bulge - 1)
                {
                    x_rela_BH_disk[i] = 0.0;
                    y_rela_BH_disk[i] = 0.0;
                    z_rela_BH_disk[i] = 0.0;

                    v_rela_BH_x_disk[i] = 0.0;
                    v_rela_BH_y_disk[i] = 0.0;
                    v_rela_BH_z_disk[i] = 0.0;
                }
            }

            for (i = 0; i < data_num[1]; i++)
            {

                if (data_bound[i][0] < N_bulge)
                {
                    x_rela_BH_bulge[i] = r_bulge[i][0] - BH_pos_x;
                    y_rela_BH_bulge[i] = r_bulge[i][1] - BH_pos_y;
                    z_rela_BH_bulge[i] = r_bulge[i][2] - BH_pos_z;

                    v_rela_BH_x_bulge[i] = v_bulge[i][0] - BH_vel_x;
                    v_rela_BH_y_bulge[i] = v_bulge[i][1] - BH_vel_y;
                    v_rela_BH_z_bulge[i] = v_bulge[i][2] - BH_vel_z;
                }

                if (data_bound[i][0] > N_bulge - 1)
                {
                    x_rela_BH_disk[i] = r_disk[i][0] - BH_pos_x;
                    y_rela_BH_disk[i] = r_disk[i][1] - BH_pos_y;
                    z_rela_BH_disk[i] = r_disk[i][2] - BH_pos_z;

                    v_rela_BH_x_disk[i] = v_disk[i][0] - BH_vel_x;
                    v_rela_BH_y_disk[i] = v_disk[i][1] - BH_vel_y;
                    v_rela_BH_z_disk[i] = v_disk[i][2] - BH_vel_z;
                }
            }

                for (i = 0; i < data_num[1]; i++)
                {
                    if (data_bound[i][0] < N_bulge)
                    {
                        L_bulge[i][0] = 0.0;
                        L_bulge[i][1] = 0.0;
                        L_bulge[i][2] = 0.0;
                    }

                    if (data_bound[i][0] > N_bulge - 1)
                    {
                        L_disk[i][0] = 0.0;
                        L_disk[i][1] = 0.0;
                        L_disk[i][2] = 0.0;
                    }
                }

                for (i = 0; i < data_num[1]; i++)
                {
                    if (data_bound[i][0] < N_bulge)
                    {
                        L_bulge[i][0] = m * 1000 * (y_rela_BH_bulge[i] * v_rela_BH_z_bulge[i] - z_rela_BH_bulge[i] * v_rela_BH_y_bulge[i]);
                        L_bulge[i][1] = m * 1000 * (z_rela_BH_bulge[i] * v_rela_BH_x_bulge[i] - x_rela_BH_bulge[i] * v_rela_BH_z_bulge[i]);
                        L_bulge[i][2] = m * 1000 * (x_rela_BH_bulge[i] * v_rela_BH_y_bulge[i] - y_rela_BH_bulge[i] * v_rela_BH_x_bulge[i]);
                    }

                    if (data_bound[i][0] > N_bulge - 1)
                    {
                        L_disk[i][0] = m * 1000 * (y_rela_BH_disk[i] * v_rela_BH_z_disk[i] - z_rela_BH_disk[i] * v_rela_BH_y_disk[i]);
                        L_disk[i][1] = m * 1000 * (z_rela_BH_disk[i] * v_rela_BH_x_disk[i] - x_rela_BH_disk[i] * v_rela_BH_z_disk[i]);
                        L_disk[i][2] = m * 1000 * (x_rela_BH_disk[i] * v_rela_BH_y_disk[i] - y_rela_BH_disk[i] * v_rela_BH_x_disk[i]);
                    }
                }

                L_bulge_tot[0] = 0.0;
                L_bulge_tot[1] = 0.0;
                L_bulge_tot[2] = 0.0;
                L_disk_tot[0] = 0.0;
                L_disk_tot[1] = 0.0;
                L_disk_tot[2] = 0.0;
                L_tot[0] = 0.0;
                L_tot[1] = 0.0;
                L_tot[2] = 0.0;

                for (i = 0; i < data_num[1]; i++)
                {
                    if (data_bound[i][0] < N_bulge)
                    {
                        L_bulge_tot[0] += L_bulge[i][0];
                        L_bulge_tot[1] += L_bulge[i][1];
                        L_bulge_tot[2] += L_bulge[i][2];
                    }

                    if (data_bound[i][0] > N_bulge - 1)
                    {
                        L_disk_tot[0] += L_disk[i][0];
                        L_disk_tot[1] += L_disk[i][1];
                        L_disk_tot[2] += L_disk[i][2];
                    }
                }

                L_tot[0] = L_bulge_tot[0] + L_disk_tot[0];
                L_tot[1] = L_bulge_tot[1] + L_disk_tot[1];
                L_tot[2] = L_bulge_tot[2] + L_disk_tot[2];

            fprintf(fpwrite1,"%f %f %f %f %f %f %f %f %f %f\n", number * 23.6, L_bulge_tot[0]*3.0857*pow(10,16)*20.74*pow(10,-20), L_bulge_tot[1]*3.0857*pow(10,16)*20.74*pow(10,-20), L_bulge_tot[2]*3.0857*pow(10,16)*20.74*pow(10,-20), L_disk_tot[0]*3.0857*pow(10,16)*20.74*pow(10,-20), L_disk_tot[1]*3.0857*pow(10,16)*20.74*pow(10,-20), L_disk_tot[2]*3.0857*pow(10,16)*20.74*pow(10,-20), L_tot[0]*3.0857*pow(10,16)*20.74*pow(10,-20), L_tot[1]*3.0857*pow(10,16)*20.74*pow(10,-20), L_tot[2]*3.0857*pow(10,16)*20.74*pow(10,-20));
            //fprintf(fpwrite1,"%f %f %f %f %f %f %f %f %f %f\n", number * 23.6, log10(L_bulge_tot[0]*3.0857*pow(10,16)*20.74), log10(L_bulge_tot[1]*3.0857*pow(10,16)*20.74), log10(L_bulge_tot[2]*3.0857*pow(10,16)*20.74), log10(L_disk_tot[0]*3.0857*pow(10,16)*20.74), log10(L_disk_tot[1]*3.0857*pow(10,16)*20.74), log10(L_disk_tot[2]*3.0857*pow(10,16)*20.74), log10(L_tot[0]*3.0857*pow(10,16)*20.74), log10(L_tot[1]*3.0857*pow(10,16)*20.74), log10(L_tot[2]*3.0857*pow(10,16)*20.74));
            fclose(fpwrite1);
            }
    }
        // fclose(fpwrite2);
        // fclose(fpwrite3);
}//3.0857\times{10}^{19}
