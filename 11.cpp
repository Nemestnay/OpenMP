#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>

using namespace std;

double sum(vector<vector<double>> matrix, int n) {
    double result = 0.0;
    #pragma omp parallel for reduction(+: result)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < (n-i); ++j) {
            sum += matrix[i][j];
        }
        result += sum;
    }
    return result;
}

double static_sum(vector<vector<double>> matrix, int n) {

    double result = 0.0;
    #pragma omp parallel for reduction(+: result) schedule(static)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < (n-i); ++j) {
            sum += matrix[i][j];
        }
        result += sum;
    }
    return result;
}

double dynamic_sum(vector<vector<double>> matrix, int n) {
    double result = 0.0;
    #pragma omp parallel for reduction(+: result) schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < (n-i); ++j) {
            sum += matrix[i][j];
        }
        result += sum;
    }
    return result;
}

double guided_sum(vector<vector<double>> matrix, int n) {

    double result = 0.0;
    #pragma omp parallel for reduction(+: result) schedule(guided)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < (n-i); ++j) {
            sum += matrix[i][j];
        }
        result += sum;
    }
    return result;
}



double auto_sum(vector<vector<double>> matrix, int n) {

    double total = 0.0;
    #pragma omp parallel for reduction(+: total) schedule(auto)
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < (n-i); ++j) {
            sum += matrix[i][j];
        }
        total += sum;
    }
    return total;
}

int main() {

    srand(static_cast<unsigned int>(time(0)));

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/11.csv");
    int n = 4096;

    vector<vector<double>> matrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < (n-i); ++j) {
            matrix[i].push_back(rand() % 10);
        }
    }

    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);

        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = sum(matrix, n);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;
        f << dt;
        if (j<16) {
            f << ',';
        }

    }
    f << '\n';

    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);

        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = static_sum(matrix, n);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;
        f << dt;
        if (j<16) {
            f << ',';
        }

    }
    f << '\n';

    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);

        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = dynamic_sum(matrix, n);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;
        f << dt;
        if (j<16) {
            f << ',';
        }

    }
    f << '\n';

    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);

        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = guided_sum(matrix, n);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;
        f << dt;
        if (j<16) {
            f << ',';
        }

    }
    f << '\n';

    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);

        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = auto_sum(matrix, n);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;
        f << dt;
        if (j<16) {
            f << ',';
        }

    }
    f << '\n';
    f.close();
    return 0;
}