// BHを重心と仮定し、束縛されている粒子を考える。

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
#define m 0.00001          // stellar粒子の質量(10^8 M_sun)
#define m_BH 0.001         // BHの質量(10^8 M_sun)
#define m_DMH 0.0001       // DMH粒子の質量(10^8 M_sun)
#define r_c 1.0            // 衛星銀河のコア半径(kpc)
#define gc_margin 0.000001 // 重心座標の相相対誤差

// stellarのデータ
#define I (550000) // 行の数
#define J (10)     // 列の数

// BHのデータ
#define K (8) // 列の数 (stellarのデータの数が100まで)

// DMHのデータ
#define L (1000000) // 行の数

double r[N][3]; // 各粒子の座標
double v[N][3]; // 各粒子の速度

double r_DMH[L][3]; // DMH粒子の座標
double v_DMH[L][3]; // DMH粒子の速速度

double v_rela_x[N]; // 重心から見た各粒子の速度
double v_rela_y[N];
double v_rela_z[N];

double v_rela_pre_x[N]; // 重心を決めなおす前の相対速度
double v_rela_pre_y[N];
double v_rela_pre_z[N];

double E_kin_rela[N]; // 重心から見た粒子の運動エネルギー
double E_pot[N];      // 粒子のポテンシャルエネルギー

double E_kin_rela_pre[N]; // 重心を決めなおす前の粒子の運動エネルギー
double E_pot_pre[N];      // 重心を決めなおす前の粒子のポテンシャルエネルギー

double E_kin_rela_new[N]; // 重心を決めなおした時の粒子の運動エネルギー
double E_pot_new[N];      // 重心を決めなおした時の粒子のポテンシャルエネルギー

double x_rela_BH[N]; // Bhから見た各粒子の相対座標
double y_rela_BH[N];
double z_rela_BH[N];

double v_rela_BH_x[N]; // Bhから見た各粒子の相対速度
double v_rela_BH_y[N];
double v_rela_BH_z[N];

double E_kin_rela_BH[N]; // BHから見た各粒子の運動エネルギー
double E_pot_BH[N];      // BHのみによる各粒子のポテンシャルエネルギー
double E_pot_DMH[N];     // DMHのみによる各粒子のポテンシャルエネルギー
double E_pot_all[N];     // DMHとBHとstellarによる各粒子のポテンシャルエネルギー

double boundptcl_pos[N][3];   // 束縛されている粒子の座標
double unboundptcl_pos[N][3]; // 束縛されていない粒子の座標
int bound_id[N];

double boundptcl_vel[N][3];   // 束縛されている粒子の速度
double unboundptcl_vel[N][3]; // 束縛されていない粒子の速度

double eps = 0.1; // ソフトニングパラメーター

int t = 0;
int i = 0;
int j = 0;
int k = 0;

int num_bound_particle = 0;    // 束縛されている粒子の数
int num_bound_particle_BH = 0; // BHを重心としたときの束縛されている粒子の数

int num_process = 0; // 試行回数

// stellarのファイルからデータを読み込む
void read_stellar(char filename1[], double mat_stellar[I][J])
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

        fscanf(fpstellar, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &mat_stellar[i][0], &mat_stellar[i][1], &mat_stellar[i][2], &mat_stellar[i][3], &mat_stellar[i][4], &mat_stellar[i][5], &mat_stellar[i][6], &mat_stellar[i][7], &mat_stellar[i][8], &mat_stellar[i][9]); /*  1行読む  */
    }

    fclose(fpstellar);
}

// BHのファイルからデータを読み込む 1粒子なので１行に一つの時間
void read_BH(char filename2[], double mat_BH[K])
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

// DMHのファイルからデータを読み込む
void read_DMH(char filename3[], double mat_DMH[L][8])
{
    int i, j;
    FILE *fpDMH;

    fpDMH = fopen(filename3, "r"); /*  読み込みモードでファイルをオープン  */
    if (fpDMH == NULL)
    {
        printf("ファイルを開くことが出来ませんでした．\n");
    }
    for (i = 0; i < L; i++)
    {
        fscanf(fpDMH, "%lf %lf %lf %lf %lf %lf %lf %lf", &mat_DMH[i][0], &mat_DMH[i][1], &mat_DMH[i][2], &mat_DMH[i][3], &mat_DMH[i][4], &mat_DMH[i][5], &mat_DMH[i][6], &mat_DMH[i][7]); /*  1行読む  */
    }
    fclose(fpDMH);
}

int main()
{
    char f_stellar[256];
    char f_BH[256];
    char f_DMH[256];
    int number;
    int i;
    double data_stellar[I][J];
    double data_BH[K];
    double data_DMH[L][8];
    FILE *fpwrite1;
    FILE *fpwrite2;
    FILE *fpwrite3;
    char *fname1[256];
    char *fname2[256];
    char *fname3[256];
    // FILE *fpwrite;
    // char *fname = "num_bound_particle_margin=100000.dat";
    // fpwrite = fopen(fname, "w");

    for (i = 0; i < N; i++)
    {
        bound_id[i] = 1;
    }
    
    for (number = 0; number < 101; number++)
    {
        if (number % 10 == 0 || number % 94 == 0 )
        {

            sprintf(fname1, "5analy_num_bound_ptcl_DMH_t=%05d.dat", number);
            fpwrite1 = fopen(("%s", fname1), "w");
            sprintf(fname2, "5analy_bound_ptcl_DMH_t=%05d.dat", number);
            fpwrite2 = fopen(("%s", fname2), "w");
            sprintf(fname3, "5analy_unbound_ptcl_DMH_t=%05d.dat", number);
            fpwrite3 = fopen(("%s", fname3), "w");
            // sprintf(filename1,"satellitedata/stellar/data_rini%05d.txt",number);
            sprintf(f_stellar, "../stellarandBH/stellar/data_rini%05d.txt", number);
            sprintf(f_BH, "../stellarandBH/stellar/BHpos-vel%05d.txt", number);
            sprintf(f_DMH, "../DMH/DMHdata_%05d.txt", number);
            read_stellar(f_stellar, data_stellar);
            read_BH(f_BH, data_BH);
            read_DMH(f_DMH, data_DMH);
            // for (i = 0; i < L; i++)
            // {
            //     printf("%f,%f,%f,%f,%f,%f\n", data_DMH[i][2], data_DMH[i][3], data_DMH[i][4], data_DMH[i][5], data_DMH[i][6], data_DMH[i][7]);
            // }

            // printf("%f,%f,%f,%f,%f,%f\n", data_BH[2], data_BH[3], data_BH[4], data_BH[5], data_BH[6], data_BH[7]);

            // 束縛されている粒子の初期化
            num_bound_particle = 0;
            num_bound_particle_BH = 0;

            // 試行回数の初期化
            num_process = 0;

            // rとvの初期化
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < 3; j++)
                {

                    r[i][j] = 0.0;
                    v[i][j] = 0.0;
                }
            }

            // 運動エネルギーとポテンシャルエネルギーの初期化
            for (i = 0; i < N; i++)
            {
                E_kin_rela[i] = 0.0;
                E_pot[i] = 0.0;
            }

            // stellarから読み込んだデータをr[i][0] ~v[i][2] に代入する
            for (i = 0; i < N; i++)
            {
                r[i][0] = data_stellar[i][2];
                r[i][1] = data_stellar[i][3];
                r[i][2] = data_stellar[i][4];
                v[i][0] = data_stellar[i][5];
                v[i][1] = data_stellar[i][6];
                v[i][2] = data_stellar[i][7];
                // printf("%f,%f,%f\n", r[i][0], r[i][1], r[i][2]);
            }

            // DMHについて
            for (i = 0; i < L; i++)
            {
                for (j = 0; j < 3; j++)
                {

                    r_DMH[i][j] = 0.0;
                    v_DMH[i][j] = 0.0;
                }
            }

            // ポテンシャルエネルギーの初期化
            for (i = 0; i < L; i++)
            {
                E_pot_DMH[i] = 0.0;
            }

            // stellarから読み込んだデータをr[i][0] ~v[i][2] に代入する
            for (i = 0; i < L; i++)
            {
                r_DMH[i][0] = data_DMH[i][2];
                r_DMH[i][1] = data_DMH[i][3];
                r_DMH[i][2] = data_DMH[i][4];
                v_DMH[i][0] = data_DMH[i][5];
                v_DMH[i][1] = data_DMH[i][6];
                v_DMH[i][2] = data_DMH[i][7];
                //printf("%f,%f,%f\n", r_DMH[i][0], r_DMH[i][1], r_DMH[i][2]);
            }

            // BHについて
            // BHの座標と速度を宣言
            double BH_pos_x;
            double BH_pos_y;
            double BH_pos_z;
            double BH_vel_x;
            double BH_vel_y;
            double BH_vel_z;
            // BHの座標、速度を初期化する
            BH_pos_x = 0.0;
            BH_pos_y = 0.0;
            BH_pos_z = 0.0;
            BH_vel_x = 0.0;
            BH_vel_y = 0.0;
            BH_vel_z = 0.0;
            // BHから読み込んだデータを代入する
            BH_pos_x = data_BH[2];
            BH_pos_y = data_BH[3];
            BH_pos_z = data_BH[4];
            BH_vel_x = data_BH[5];
            BH_vel_y = data_BH[6];
            BH_vel_z = data_BH[7];
            // printf("%f %f %f %f",r[0][0],v[0][0],BH_pos_x,BH_vel_x);
            //  printf("\n");

            // BHから見た各粒子の相対座標と相対速度を決める
            for (i = 0; i < N; i++)
            {

                x_rela_BH[i] = 0.0;
                y_rela_BH[i] = 0.0;
                z_rela_BH[i] = 0.0;

                v_rela_BH_x[i] = 0.0;
                v_rela_BH_y[i] = 0.0;
                v_rela_BH_z[i] = 0.0;
            }
            for (i = 0; i < N; i++)
            {

                x_rela_BH[i] = r[i][0] - BH_pos_x;
                y_rela_BH[i] = r[i][1] - BH_pos_y;
                z_rela_BH[i] = r[i][2] - BH_pos_z;

                v_rela_BH_x[i] = v[i][0] - BH_vel_x;
                v_rela_BH_y[i] = v[i][1] - BH_vel_y;
                v_rela_BH_z[i] = v[i][2] - BH_vel_z;
            }

            //  BHからみた各粒子の運動エネルギーとポテンシャルエネルギーを求める。
            // 運動エネルギーとポテンシャルエネルギーの初期化
            for (i = 0; i < N; i++)
            {

                E_kin_rela_BH[i] = 0.0;
                E_pot_BH[i] = 0.0;
                E_pot_all[i] = 0.0;
            }

            // 運動エネルギーを求める
            for (i = 0; i < N; i++)
            {
                E_kin_rela_BH[i] = 0.5 * m * (v_rela_BH_x[i] * v_rela_BH_x[i] + v_rela_BH_y[i] * v_rela_BH_y[i] + v_rela_BH_z[i] * v_rela_BH_z[i]);
            }

            // ポテンシャルエネルギーを求める
            // まずstellarのみのポテンシャルエネルギー
            // #pragma omp parallel for
            // for (i = 0; i < N; i++)
            // {
            //     // E_pot[i] = -G * m * m / eps;

            //     E_pot[i]=G*m*m/eps;
            //     for (j = 0; j < N; j++)
            //     {
            //         E_pot[i] += (-G * m * m) / sqrt((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps * eps);
            //     }
            // }

            // ポテンシャルエネルギーを求める
            for (i = 0; i < N; i++)
            {
                E_pot[i] = G * m * m / eps;
                for (j = 0; j < N; j++)
                {
                    E_pot[i] += -G * m * m / sqrt((r[i][0] - r[j][0]) * (r[i][0] - r[j][0]) + (r[i][1] - r[j][1]) * (r[i][1] - r[j][1]) + (r[i][2] - r[j][2]) * (r[i][2] - r[j][2]) + eps * eps);
                }
            }

            // 次にBHのみによるポテンシャルエネルギー
            for (i = 0; i < N; i++)
            {
                E_pot_BH[i] = (-G * m * m_BH / sqrt(x_rela_BH[i] * x_rela_BH[i] + y_rela_BH[i] * y_rela_BH[i] + z_rela_BH[i] * z_rela_BH[i]));
                // + (-G * m * m_BH / eps);
            }

            // DMHのみによるポテンシャルエネルギー
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < L; j++)
                {
                    E_pot_DMH[i] += -G * m * m_DMH / sqrt((r[i][0] - r_DMH[j][0]) * (r[i][0] - r_DMH[j][0]) + (r[i][1] - r_DMH[j][1]) * (r[i][1] - r_DMH[j][1]) + (r[i][2] - r_DMH[j][2]) * (r[i][2] - r_DMH[j][2]));
                }
            }

            // 各粒子のポテンシャルエネルギー
            for (i = 0; i < N; i++)
            {
                E_pot_all[i] = E_pot[i] + E_pot_BH[i] + E_pot_DMH[i];
            }

            // for (i = 0; i < N; i++)
            // {

            //     printf("%f %.20f %.20f %.20f %.20f\n", E_kin_rela_BH[i], E_pot_all[i], E_pot[i], E_pot_BH[i], E_pot_DMH[i]); // ptclout.csv
            // }

            // ここまでで各粒子の運動エネルギーとポテンシャルエネルギーが求まった

            // bound粒子の判定を行う
            // 束縛されている粒子と束縛されていない粒子の座標と速度の初期化
            for (i = 0; i < N; i++)
            {
                for (k = 0; k < 3; k++)
                {

                    boundptcl_pos[i][k] = 0.0;
                    boundptcl_vel[i][k] = 0.0;
                    unboundptcl_pos[i][k] = 0.0;
                    unboundptcl_vel[i][k] = 0.0;
                }
            }

            // 束縛されている粒子と束縛されていない粒子の座標と速度の代入
            for (i = 0; i < N; i++)
            {
                if (bound_id[i] == 1)//以前の計算で束縛されていない粒子を除外する。
                {
                    if (E_kin_rela_BH[i] + E_pot_all[i] < 0.0)
                    {
                        for (k = 0; k < 3; k++)
                        {

                            boundptcl_pos[i][k] = r[i][k];
                            boundptcl_vel[i][k] = v[i][k];
                        }
                        num_bound_particle_BH++;
                    }
                    else
                    {
                        for (k = 0; k < 3; k++)
                        {

                            unboundptcl_pos[i][k] = r[i][k];
                            unboundptcl_vel[i][k] = v[i][k];
                        }
                    }
                }
               
               if(bound_id[i]==0)
                {
                        for (k = 0; k < 3; k++)
                        {

                            unboundptcl_pos[i][k] = r[i][k];
                            unboundptcl_vel[i][k] = v[i][k];
                        } 
                }
            }
            // unboundptclのIDを特定し、以降カウントされないようにする。
 for(i=0;i<N;i++)
 {

     if (E_kin_rela_BH[i] + E_pot_all[i] > 0.0)
     {
         bound_id[i] = 0;
     }
 }

 fprintf(fpwrite1, "%d %d\n", number, num_bound_particle_BH);
 // 必要であれば粒子のＩＤも
 for (i = 0; i < N; i++)
 {

     if (E_kin_rela_BH[i] + E_pot_all[i] < 0.0)
     {
         fprintf(fpwrite2, "%d %f %f %f %f %f %f\n", i, boundptcl_pos[i][0], boundptcl_pos[i][1], boundptcl_pos[i][2], boundptcl_vel[i][0], boundptcl_vel[i][1], boundptcl_vel[i][2]);
     }
 }
 for (i = 0; i < N; i++)
 {
     if (E_kin_rela_BH[i] + E_pot_all[i] > 0.0)
     {
         fprintf(fpwrite3, "%d %f %f %f %f %f %f\n", i, unboundptcl_pos[i][0], unboundptcl_pos[i][1], unboundptcl_pos[i][2], unboundptcl_vel[i][0], unboundptcl_vel[i][1], unboundptcl_vel[i][2]);
     }
 }


 fclose(fpwrite1);
 fclose(fpwrite2);
 fclose(fpwrite3);
        }
    }
}
