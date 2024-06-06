#include <iostream>
#include <vector>
#include <omp.h>
#include <fstream>

using namespace std;

double minmax_without_parallel(vector<vector<double>> matrix, int n) {
    double result = matrix[0][0];
    for (int i = 0; i < n; ++i) {
        double minimum = matrix[i][0];
        for (int j = 1; j < n; ++j) {
            minimum = min(minimum, matrix[i][j]);
        }
        result = max(result, minimum);
    }
    return result;
}

double minmax(vector<vector<double>> matrix, int n){
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


double shedule(vector<vector<double>> matrix, int n, string s) {
    double result = matrix[0][0];
    if (s == "static") {
        #pragma omp parallel for reduction(max: result) schedule(static)
        for (int i = 0; i < n; ++i) {
            double minimum = matrix[i][0];
            for (int j = 1; j < n; ++j) {
                minimum = min(minimum, matrix[i][j]);
            }
            result = max(result, minimum);
        }
    }

    if (s == "dynamic") {
    #pragma omp parallel for reduction(max: result) schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            double minimum = matrix[i][0];
            for (int j = 1; j < n; ++j) {
                minimum = min(minimum, matrix[i][j]);
            }
            result = max(result, minimum);
        }
    }

    if (s == "guided") {
    #pragma omp parallel for reduction(max: result) schedule(guided)
        for (int i = 0; i < n; ++i) {
            double minimum = matrix[i][0];
            for (int j = 1; j < n; ++j) {
                minimum = min(minimum, matrix[i][j]);
            }
            result = max(result, minimum);
        }
    }

    if (s == "auto") {
    #pragma omp parallel for reduction(max: result) schedule(auto)
        for (int i = 0; i < n; ++i) {
            double minimum = matrix[i][0];
            for (int j = 1; j < n; ++j) {
                minimum = min(minimum, matrix[i][j]);
            }
            result = max(result, minimum);
        }
    }
    return result;
}


int main() {

    srand(static_cast<unsigned int>(time(0)));

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/5.csv");
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
        double result =minmax_without_parallel(matrix, n);
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
        double result = minmax(matrix, n);
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
        string s = "static";
        double result = shedule(matrix, n, s);
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
        string s = "dynamic";
        double result = shedule(matrix, n, s);
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
        string s = "guided";
        double result = shedule(matrix, n, s);
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
        string s = "auto";
        double result = shedule(matrix, n, s);
        t2 = omp_get_wtime ();
        dt = (t2 - t1) * 1000000;

        f << dt;
        if (j<16) {
            f << ',';
        }
    }
    f << '\n';

    f.close();
}