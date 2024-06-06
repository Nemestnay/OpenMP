#include <iostream>
#include <vector>
#include <omp.h>
#include <chrono>
#include <fstream>

using namespace std;


double minmax(int n, vector<vector<double>> matrix) {
    double result = matrix[0][0];
#pragma omp parallel for reduction(max: result)
    for (int i = 0; i < n; ++i) {
        double minimum = matrix[i][0];
        for (int j = 1; j < n; ++j) {
            minimum = min(minimum, matrix[i][j]);
        }
        result = max(result, minimum);
    }
    return result;
}

double parallel_minmax(vector<vector<double>> matrix, int n, int j) {
    double result = matrix[0][0];

    omp_set_nested(true);
    omp_set_num_threads(j);
    #pragma omp parallel for reduction(max: result) num_threads(1)
    for (int i = 0; i < n; ++i) {
        double minimum = matrix[i][0];

        omp_set_num_threads(j);
        #pragma omp parallel for reduction(min: minimum) num_threads(j)
        for (int j = 1; j < n; ++j) {
            minimum = min(minimum, matrix[i][j]);
        }
        result = max(result, minimum);
    }
    return result;
}

int main() {

    srand(static_cast<unsigned int>(time(0)));

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/10_2.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    ofstream g;
    g.open ("/home/alena/CLionProjects/OpenMP/tabl/10_1.csv");
    g << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    for(int n = 4; n <= 4096; n *= 4) {
        vector<vector<double>> matrix(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix[i].push_back(rand() % 100 + 1);
            }
        }

        for (int j = 1; j <= 16; ++j) {
            omp_set_num_threads(j);

            double t1, t2, dt;
            t1 = omp_get_wtime ();
            double result = parallel_minmax(matrix, n, j);
            t2 = omp_get_wtime ();
            dt = (t2 - t1) * 1000000;

            f << dt;
            if (j<16) {
                f << ',';
            }



            t1 = omp_get_wtime ();
            result = minmax(n, matrix);
            t2 = omp_get_wtime ();
            dt = (t2 - t1) * 1000000;

            g << dt;
            if (j<16) {
                g << ',';
            }
        }
        f << '\n';
        g << '\n';
    }
    f.close();
    g.close();
}