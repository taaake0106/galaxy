// BHを重心と仮定し、重心を中心とした束縛された粒子のhalfmassradiusを考える。

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define N 182199           // 粒子数
#define M_sat 1            // 衛星銀河の質量(×10^9 M_sun)
#define G 1.0              // 計算機内の重力定数G’
#define m 0.00001          // 衛星銀河を構成する粒子の質量(10^8 M_sun)
#define m_BH 0.001         // BHの質量(10^8 M_sun)
#define r_c 1.0            // 衛星銀河のコア半径(kpc)
#define gc_margin 0.000001 // 重心座標の相相対誤差
#define eps 0.1            // ソフトニングパラメーター

// stellarのデータ
#define I (182199) // 行の数
#define J (7)      // 列の数

// BHのデータ
#define K (200) // 行の数 (stellarのデータの数が100まで)
#define L (8)   // 列の数

// #define delta_x 0.005
// #define delta_y 0.005

#define eps 0.0000000001

int i;
int j;

double x;
double y;

int count; // 1binあたりの粒子数

double v_z_sum;
double v_z;
double v_z_bin_sum;
double v_z_bin;
double v_rot;

// bound-ptclのデータを読み込む
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

        fscanf(fpstellar, "%lf %lf %lf %lf %lf %lf %lf", &mat_bound[i][0], &mat_bound[i][1], &mat_bound[i][2], &mat_bound[i][3], &mat_bound[i][4], &mat_bound[i][5], &mat_bound[i][6]); /*  1行読む  */
    }

    fclose(fpstellar);
}

int main()
{

    // PI=M_PI;
    char f_bound[256];
    char f_BH[256];
    char f_num[256];
    int number;
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

    for (number = 94; number < 95; number++)
    {
        sprintf(fname1, "2velrot_x-y_t=%05d.dat", number);
        sprintf(fname3, "max.min_t=%05d.dat", number);
        fpwrite1 = fopen(("%s", fname1), "w");
        fpwrite3 = fopen(("%s", fname3), "w");
        sprintf(f_bound, "../../stellar.BH.DMH/5analy_bound_ptcl_DMH_t=%05d.dat", number);
        read_bound(f_bound, data_bound);

        //         for(i=0;i<N;i++){

        // fprintf(fpwrite2,"%f %f %f %f %f %f\n",data_bound[i][1],data_bound[i][2],data_bound[i][3],data_bound[i][4],data_bound[i][5],data_bound[i][6]);
        //         }

double max_x = -100000.0;
double max_y = -100000.0;
double min_x = 100000.0;
double min_y = 100000.0;

for(i=0;i<N;i++){
    
    if(data_bound[i][1] > max_x){
        max_x = data_bound[i][1];
    }

    if(data_bound[i][2] > max_y){
        max_y = data_bound[i][2];
    }

    if(data_bound[i][1] < min_x){
        min_x = data_bound[i][1];
    }

    if(data_bound[i][2] < min_y){
        min_y = data_bound[i][2];
    }

}
fprintf(fpwrite3,"%f %f %f %f",min_x,max_x,min_y,max_y);



min_x = -0.3037;
max_x = -0.2963;

min_y = -0.9037;
max_y = -0.8963;

double delta_x = 0.0;
double delta_y = 0.0;

delta_x = (max_x - min_x)/100;
delta_y = (max_y - min_y)/100;

        v_z_sum = 0.0;
        v_z = 0.0;
        
        for (i = 0; i < N; i++)
        {
            v_z_sum += data_bound[i][6];
        }
        v_z = v_z_sum / N;

        for (x = min_x; x < (max_x + eps); x += delta_x)
        {
            for (y = min_y; y < (max_y + eps); y += delta_y)
            {
                count = 0;
                v_z_bin_sum = 0.0;
                v_z_bin = 0.0;
                for (i = 0; i < N; i++)
                {
                   

                    if (data_bound[i][1] > x && data_bound[i][1] < x + delta_x && data_bound[i][2] > y && data_bound[i][2] < y + delta_y)
                    {
                        count++;
                        v_z_bin_sum += data_bound[i][6];
                    }
                    
                }
                if (count > 0)
                {
                    v_z_bin = v_z_bin_sum / count;
                }
                
                if(log10((count + eps)*1000/(delta_x*delta_y)) > 8.0){
                v_rot = v_z_bin - v_z;
                fprintf(fpwrite1, "%f %f %f\n", x, y, v_rot*20.74);
                }
                // v_rot = v_z_bin - v_z;
                // fprintf(fpwrite1, "%f %f %f\n", x, y, v_rot*20.74);
                
                // fprintf(fpwrite1, "%d %d %f\n", x, y, count*1000);
            }
            fprintf(fpwrite1, "\n");
        }
    }
    fclose(fpwrite1);
    fclose(fpwrite3);
}
