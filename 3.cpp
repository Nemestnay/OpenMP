#include <iostream>
#include <omp.h>
#include <fstream>
#include <math.h>

using namespace std;

double f(double x) {
    return sin(x);
}

double integral(int n) {
    double a = 0.0;
    double b = 1;
    double h = (b - a) / n;
    double sum = 0.0;

    #pragma omp parallel for reduction(+: sum)
    for (int i = 0; i < n; ++i) {
        double seredina = a + h * (i - 0.5);
        sum += h * f(seredina);
    }
    return sum;
}

int main() {
    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/3.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    for(int i = 16; i <= 16777216; i*=16) {
        int n = i;
        for (int j = 1; j <= 16; ++j) {
            omp_set_num_threads(j);

            double t1, t2, dt;
            t1 = omp_get_wtime ();
            double result = integral(n);
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
