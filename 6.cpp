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
/*
double f(double x, double y) {
    if (x > y) {return x;}
    else {return y;}
}
double atomic(vector<vector<double>> matrix, int n){
    double result = matrix[0][0];
#pragma omp parallel shared(matrix)
    for (int i = 0; i < n; ++i) {
        double minimum = matrix[i][0];
        for (int j = 1; j < n; ++j) {
            minimum = min(minimum, matrix[i][j]);
        }
#pragma omp atomic
        result = f(result, minimum);
    }
    return result;
}
*/
double crit_minmax(vector<vector<double>> matrix, int n) {
    double result = matrix[0][0];
    omp_lock_t writelock;
    omp_init_lock(&writelock);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double minimum = matrix[i][0];
        for (int j = 1; j < n; ++j) {
            minimum = min(minimum, matrix[i][j]);
        }

        omp_set_lock(&writelock);
        result = max(result, minimum);
        omp_unset_lock(&writelock);
    }

    omp_destroy_lock(&writelock);
    return result;
}


double lock(vector<vector<double>> matrix, int n) {
        double result = matrix[0][0];
        omp_lock_t writelock;
        omp_init_lock(&writelock);

#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double minimum = matrix[i][0];
            for (int j = 1; j < n; ++j) {
                minimum = min(minimum, matrix[i][j]);
            }

            omp_set_lock(&writelock);
            result = max(result, minimum);
            omp_unset_lock(&writelock);
        }

        omp_destroy_lock(&writelock);
        return result;
    }


int main() {

    srand(static_cast<unsigned int>(time(0)));

    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/6.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";

    int n = 2048*2;

    vector<vector<double>> matrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i].push_back(rand() % 10);
        }
    }

    for (int j = 1; j <= 16; ++j) {
        omp_set_num_threads(j);
        double t1, t2, dt;
        t1 = omp_get_wtime ();
        double result = minmax_without_parallel(matrix, n);
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
        double result = crit_minmax(matrix, n);
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
        double result = lock(matrix, n);
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