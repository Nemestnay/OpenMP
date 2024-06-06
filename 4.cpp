#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>

using namespace std;

double minmax(int n, double** matrix) {
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

int main() {
   // omp_set_num_threads(3);

    srand(static_cast<unsigned int>(time(0)));

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/4.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    for(int n = 16; n <= 4096 * 4; n*=4) {
        double** matrix = new double*[n];
            for (int i = 0; i < n; ++i) {
                matrix[i] = new double[n];
                for (int j = 0; j < n; ++j) {
                    matrix[i][j] = rand() % 100 + 1;
                }
            }

            for (int j = 1; j <= 16; ++j) {

                omp_set_num_threads(j);

                double t1, t2, dt;
                t1 = omp_get_wtime ();
                double result = minmax(n, matrix);
                t2 = omp_get_wtime ();
                dt = (t2 - t1) * 1000000;

                f << dt;
                if (j<16) {
                    f << ',';
                }
            }
            f << '\n';
        }

    f.close();
}