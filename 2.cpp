#include <iostream>
#include <vector>
#include <omp.h>
#include <random>
#include <fstream>

using namespace std;

double scaliar(vector<double> a, vector<double> b) {
    double result = 0.0;
    #pragma omp parallel for reduction(+: result)
    for (int i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

int main() {

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/2.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";

    for(int i = 16; i <= 16777216; i*=16) {
        vector<double> v1(i);
        vector<double> v2(i);

        for (int k = 0; k < i; ++k) {
            v1[k] = dis(gen);
            v2[k] = dis(gen);
        }

        for (int j = 1; j <= 16; ++j) {
            omp_set_num_threads(j);

            double t1, t2, dt;
            t1 = omp_get_wtime ();
            double result = scaliar(v1, v2);
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

    return 0;
}