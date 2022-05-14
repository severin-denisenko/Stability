/*
Runge-Kutta 4 Method to solve the trajectory
of N - 1 particles and 1 planet under the gravitational force.

Input and output in meters, kilograms and seconds.

* Initial data: "data.csv"
The 7 columns correspond to: rx, ry, rz, vx, vy, vz, m.
Each row corresponds to the data of each particle.

* Result writes in data: "results_R.txt" i "results_v.txt"
Each row corresponds to an instant of time (n-iteration).
The data is separated by blank spaces and correspond to:
rx1, ry1, rz1, ... , rxN, ryN, rzN   ("results_R.txt")
vx1, vy1, vz1, ... , vxN, vyN, vzN   ("results_V.txt")
where the order of the particles is the same as "data.csv".
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define R_NORM 1.00000e+11
#define V_NORM 1.00000e+04
#define M_NORM 6.67384e-30
#define T_NORM 3.15000e+00
#define G 6.67384e-11

//#define MAKE_OUT 0

double solve(int N, int i, int j, long double *R, long double *V, long double *M);

long double dnorm(long double *r1, long double *r2)
{
    return sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]));
}

int main(int argc, char *argv[])
{
    int N = 3;
    // Read data //

    long double R[3 * N]; // rx1, ry1, rz1, ... , rxN, ryN, rzN
    long double V[3 * N]; // vx1, vy1, vz1, ... , vxN, vyN, vzN
    long double M[N];     // m1, ... , mN

    FILE *input;
    input = fopen("data.csv", "r");
    for (int i = 0; i < N - 1; i++)
    {
        long double read;

        fscanf(input, "%Le,", &read);
        R[3 * i] = read / R_NORM;
        fscanf(input, "%Le,", &read);
        R[3 * i + 1] = read / R_NORM;
        fscanf(input, "%Le,", &read);
        R[3 * i + 2] = read / R_NORM;
        fscanf(input, "%Le,", &read);
        V[3 * i] = read / V_NORM;
        fscanf(input, "%Le,", &read);
        V[3 * i + 1] = read / V_NORM;
        fscanf(input, "%Le,", &read);
        V[3 * i + 2] = read / V_NORM;
        fscanf(input, "%Le\n", &read);
        M[i] = read * M_NORM;
    }
    fclose(input);

    // MODELLING STABILITY //

    int resolution = 100;

    long double **LifeTime = (long double **)malloc(sizeof(long double *) * resolution);
    for (int i = 0; i < resolution; i++)
    {
        LifeTime[i] = (long double *)malloc(sizeof(long double) * resolution);
    }

    for (int i = 0; i < resolution; i++)
    {
        for (int j = 0; j < resolution; j++)
        {
            // Adding N'th element -- planet //

            long double a = 8700 / 2 * 1.495978707e11 * (i + 1) / resolution;
            long double e = 0.5 * (j + 1) / resolution;

            R[3 * 2] = R[3 * 1];
            R[3 * 2 + 1] = R[3 * 1 + 1];
            R[3 * 2 + 2] = R[3 * 1 + 2] - a * (1 - e) / R_NORM;
            V[3 * 2] = V[3 * 1] - sqrt(G * M[1] / a * (1 + e) / (1 - e)) / V_NORM;
            V[3 * 2 + 1] = V[3 * 1 + 1] ;
            V[3 * 2 + 2] = V[3 * 1 + 2];
            M[2] = 2.4298E+10 * M_NORM;

            LifeTime[i][j] = solve(N, i, j, R, V, M);
            printf("done %d, %d.\n", i, j);
        }
    }

    // Output main result //

    FILE *main_result;

    main_result = fopen("result.dat", "w");

    for (int i = 0; i < resolution; i++)
    {
        for (int j = 0; j < resolution; j++)
        {
            fprintf(main_result, " %.5Le ", LifeTime[i][j]);
        }
        fprintf(main_result, "\n");
    }

    for (int i = 0; i < resolution; i++)
    {
        free(LifeTime[i]);
    }
    free(LifeTime);

    fclose(main_result);

    return 0;
}

// Gets i,j -- numbers,  N -- numbers of bodes, vectors of R, V and M
// Returns the T -- lifetime, or 0 if can't determine maximum
double solve(int N, int in, int jn, long double *R_n0, long double *V_n0, long double *M)
{
    // 550000 yr -- Period of Proxima around AB
    long double Number_of_periods = 50;
    long double Total_time = 550000.0 * Number_of_periods * T_NORM;

    // 1 yr -- optimal
    long double Delta_time = T_NORM * 10;

    int Number_of_iterations = Total_time / Delta_time;
    int Output_iteration = Number_of_iterations / 1000;

    #ifdef MAKE_OUT
    FILE *output_r;
    FILE *output_v;

    char filename_r[100];
    char filename_v[100];

    sprintf(filename_r, "result/result_%d_%d_R.txt", in, jn);
    sprintf(filename_v, "result/result_%d_%d_V.txt", in, jn);

    output_r = fopen(filename_r, "w");
    output_v = fopen(filename_v, "w");
    #endif

    // Copying data

    long double R_n1[3 * N];
    long double V_n1[3 * N];

    for (int i = 0; i < N; i++)
    {
        R_n1[i * 3] = R_n0[i * 3];
        R_n1[i * 3 + 1] = R_n0[i * 3 + 1];
        R_n1[i * 3 + 2] = R_n0[i * 3 + 2];
        V_n1[i * 3] = V_n0[i * 3];
        V_n1[i * 3 + 1] = V_n0[i * 3 + 1];
        V_n1[i * 3 + 2] = V_n0[i * 3 + 2];
    }

    // RK4 //

    long double R_n2[3 * N];
    long double V_n2[3 * N];
    long double Ria;
    long double k1R[3 * N], k2R[3 * N], k3R[3 * N], k4R[3 * N];
    long double k1V[3 * N], k2V[3 * N], k3V[3 * N], k4V[3 * N];

    #ifdef MAKE_OUT
    // print data
    fprintf(output_r, "%.10Le", R_n1[0] * R_NORM);
    fprintf(output_v, "%.10Le", V_n1[0] * V_NORM);
    for (int i = 1; i < 3 * N; i++)
    {
        fprintf(output_r, " %.10Le", R_n1[i] * R_NORM);
        fprintf(output_v, " %.10Le", V_n1[i] * V_NORM);
    }
    fprintf(output_r, "\n");
    fprintf(output_v, "\n");
    #endif

    for (int n = 1; n < Number_of_iterations + 1; n++)
    {

        // k1
        for (int i = 0; i < 3 * N; i++)
        {
            k1R[i] = Delta_time * V_n1[i];
        }
        for (int i = 0; i < N; i++)
        {
            k1V[3 * i] = 0;     // x
            k1V[3 * i + 1] = 0; // y
            k1V[3 * i + 2] = 0; // z
            for (int a = 0; a < N; a++)
            {
                if (a != i)
                {
                    Ria = sqrt((R_n1[3 * a] - R_n1[3 * i]) * (R_n1[3 * a] - R_n1[3 * i]) +
                               (R_n1[3 * a + 1] - R_n1[3 * i + 1]) * (R_n1[3 * a + 1] - R_n1[3 * i + 1]) +
                               (R_n1[3 * a + 2] - R_n1[3 * i + 2]) * (R_n1[3 * a + 2] - R_n1[3 * i + 2]));
                    k1V[3 * i] += M[a] * (R_n1[3 * a] - R_n1[3 * i]) / (Ria * Ria * Ria);
                    k1V[3 * i + 1] += M[a] * (R_n1[3 * a + 1] - R_n1[3 * i + 1]) / (Ria * Ria * Ria);
                    k1V[3 * i + 2] += M[a] * (R_n1[3 * a + 2] - R_n1[3 * i + 2]) / (Ria * Ria * Ria);
                }
            }
            k1V[3 * i] = k1V[3 * i] * Delta_time;         // x
            k1V[3 * i + 1] = k1V[3 * i + 1] * Delta_time; // y
            k1V[3 * i + 2] = k1V[3 * i + 2] * Delta_time; // z
        }

        // k2
        for (int i = 0; i < 3 * N; i++)
        {
            k2R[i] = Delta_time * (V_n1[i] + 0.5 * k1V[i]);
        }
        for (int i = 0; i < N; i++)
        {
            k2V[3 * i] = 0;     // x
            k2V[3 * i + 1] = 0; // y
            k2V[3 * i + 2] = 0; // z
            for (int a = 0; a < N; a++)
            {
                if (a != i)
                {
                    Ria = sqrt((R_n1[3 * a] + 0.5 * k1R[3 * a] - R_n1[3 * i] - 0.5 * k1R[3 * i]) * (R_n1[3 * a] + 0.5 * k1R[3 * a] - R_n1[3 * i] - 0.5 * k1R[3 * i]) +
                               (R_n1[3 * a + 1] + 0.5 * k1R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k1R[3 * i + 1]) * (R_n1[3 * a + 1] + 0.5 * k1R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k1R[3 * i + 1]) +
                               (R_n1[3 * a + 2] + 0.5 * k1R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k1R[3 * i + 2]) * (R_n1[3 * a + 2] + 0.5 * k1R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k1R[3 * i + 2]));
                    k2V[3 * i] += M[a] * (R_n1[3 * a] + 0.5 * k1R[3 * a] - R_n1[3 * i] - 0.5 * k1R[3 * i]) / (Ria * Ria * Ria);
                    k2V[3 * i + 1] += M[a] * (R_n1[3 * a + 1] + 0.5 * k1R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k1R[3 * i + 1]) / (Ria * Ria * Ria);
                    k2V[3 * i + 2] += M[a] * (R_n1[3 * a + 2] + 0.5 * k1R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k1R[3 * i + 2]) / (Ria * Ria * Ria);
                }
            }
            k2V[3 * i] = k2V[3 * i] * Delta_time;         // x
            k2V[3 * i + 1] = k2V[3 * i + 1] * Delta_time; // y
            k2V[3 * i + 2] = k2V[3 * i + 2] * Delta_time; // z
        }

        // k3
        for (int i = 0; i < 3 * N; i++)
        {
            k3R[i] = Delta_time * (V_n1[i] + 0.5 * k2V[i]);
        }
        for (int i = 0; i < N; i++)
        {
            k3V[3 * i] = 0;     // x
            k3V[3 * i + 1] = 0; // y
            k3V[3 * i + 2] = 0; // z
            for (int a = 0; a < N; a++)
            {
                if (a != i)
                {
                    Ria = sqrt((R_n1[3 * a] + 0.5 * k2R[3 * a] - R_n1[3 * i] - 0.5 * k2R[3 * i]) * (R_n1[3 * a] + 0.5 * k2R[3 * a] - R_n1[3 * i] - 0.5 * k2R[3 * i]) +
                               (R_n1[3 * a + 1] + 0.5 * k2R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k2R[3 * i + 1]) * (R_n1[3 * a + 1] + 0.5 * k2R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k2R[3 * i + 1]) +
                               (R_n1[3 * a + 2] + 0.5 * k2R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k2R[3 * i + 2]) * (R_n1[3 * a + 2] + 0.5 * k2R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k2R[3 * i + 2]));
                    k3V[3 * i] += M[a] * (R_n1[3 * a] + 0.5 * k2R[3 * a] - R_n1[3 * i] - 0.5 * k2R[3 * i]) / (Ria * Ria * Ria);
                    k3V[3 * i + 1] += M[a] * (R_n1[3 * a + 1] + 0.5 * k2R[3 * a + 1] - R_n1[3 * i + 1] - 0.5 * k2R[3 * i + 1]) / (Ria * Ria * Ria);
                    k3V[3 * i + 2] += M[a] * (R_n1[3 * a + 2] + 0.5 * k2R[3 * a + 2] - R_n1[3 * i + 2] - 0.5 * k2R[3 * i + 2]) / (Ria * Ria * Ria);
                }
            }
            k3V[3 * i] = k3V[3 * i] * Delta_time;         // x
            k3V[3 * i + 1] = k3V[3 * i + 1] * Delta_time; // y
            k3V[3 * i + 2] = k3V[3 * i + 2] * Delta_time; // z
        }

        // k4
        for (int i = 0; i < 3 * N; i++)
        {
            k4R[i] = Delta_time * (V_n1[i] + k3V[i]);
        }
        for (int i = 0; i < N; i++)
        {
            k4V[3 * i] = 0;     // x
            k4V[3 * i + 1] = 0; // y
            k4V[3 * i + 2] = 0; // z
            for (int a = 0; a < N; a++)
            {
                if (a != i)
                {
                    Ria = sqrt((R_n1[3 * a] + k3R[3 * a] - R_n1[3 * i] - k3R[3 * i]) * (R_n1[3 * a] + k3R[3 * a] - R_n1[3 * i] - k3R[3 * i]) +
                               (R_n1[3 * a + 1] + k3R[3 * a + 1] - R_n1[3 * i + 1] - k3R[3 * i + 1]) * (R_n1[3 * a + 1] + k3R[3 * a + 1] - R_n1[3 * i + 1] - k3R[3 * i + 1]) +
                               (R_n1[3 * a + 2] + k3R[3 * a + 2] - R_n1[3 * i + 2] - k3R[3 * i + 2]) * (R_n1[3 * a + 2] + k3R[3 * a + 2] - R_n1[3 * i + 2] - k3R[3 * i + 2]));
                    k4V[3 * i] += M[a] * (R_n1[3 * a] + k3R[3 * a] - R_n1[3 * i] - k3R[3 * i]) / (Ria * Ria * Ria);
                    k4V[3 * i + 1] += M[a] * (R_n1[3 * a + 1] + k3R[3 * a + 1] - R_n1[3 * i + 1] - k3R[3 * i + 1]) / (Ria * Ria * Ria);
                    k4V[3 * i + 2] += M[a] * (R_n1[3 * a + 2] + k3R[3 * a + 2] - R_n1[3 * i + 2] - k3R[3 * i + 2]) / (Ria * Ria * Ria);
                }
            }
            k4V[3 * i] = k4V[3 * i] * Delta_time;         // x
            k4V[3 * i + 1] = k4V[3 * i + 1] * Delta_time; // y
            k4V[3 * i + 2] = k4V[3 * i + 2] * Delta_time; // z
        }

        // R & V
        for (int i = 0; i < 3 * N; i++)
        {
            R_n2[i] = R_n1[i] + (k1R[i] + 2 * k2R[i] + 2 * k3R[i] + k4R[i]) / 6.;
            V_n2[i] = V_n1[i] + (k1V[i] + 2 * k2V[i] + 2 * k3V[i] + k4V[i]) / 6.;
        }

        #ifdef MAKE_OUT
        // Save result
        if (n % Output_iteration == 0)
        {
            fprintf(output_r, "%.10Le", R_n2[0] * R_NORM);
            fprintf(output_v, "%.10Le", V_n2[0] * V_NORM);
            for (int i = 1; i < 3 * N; i++)
            {
                fprintf(output_r, " %.10Le", R_n2[i] * R_NORM);
                fprintf(output_v, " %.10Le", V_n2[i] * V_NORM);
            }
            fprintf(output_r, "\n");
            fprintf(output_v, "\n");
        }
        #endif

        // OLD = NEW
        for (int i = 0; i < 3 * N; i++)
        {
            R_n1[i] = R_n2[i];
            V_n1[i] = V_n2[i];
        }

        if (dnorm(&R_n1[6], &R_n1[3]) * R_NORM > (1.3015e+15 * 10))
        {
            return n * Delta_time / T_NORM;
        }
    }

    #ifdef MAKE_OUT
    fclose(output_r);
    fclose(output_v);
    #endif

    return Total_time / T_NORM;
}
