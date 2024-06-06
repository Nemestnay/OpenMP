#include <iostream>
#include <vector>
#include <omp.h>
#include <random>
#include <fstream>

using namespace std;


double multi_scaliar(const vector<vector<double>>& vectors) {
    int num_vectors = vectors.size();
    int vector_size = vectors[0].size();
    double result = 0.0;
    #pragma omp parallel
    {
        #pragma omp sections reduction(+:result)
        {
            #pragma omp section
            {
                double partial_sum = 0.0;
                for (int i = 0; i < vector_size; ++i) {
                    for (int j = 0; j < num_vectors; ++j) {
                        partial_sum += vectors[j][i];
                    }
                }
                result += partial_sum;
            }
        }
    }
    return result;
}

vector<vector<double>> generateRandomVectors(int n, int vectorSize) {
    random_device rd;
    std::mt19937 gen(rd());

    vector<vector<double>> randomVectors;

    for (int i = 0; i < n; ++i) {
        vector<double> vec;
        for (int j = 0; j < vectorSize; ++j) {
            uniform_real_distribution<double> dis(0.0, 1.0);
            vec.push_back(dis(gen));
        }
        randomVectors.push_back(vec);
    }
    return randomVectors;
}


int main() {
    int vectorSize = 100;
    ofstream f;
    f.open ("/home/alena/CLionProjects/OpenMP/tabl/7.csv");
    f << "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16\n";
    for(int n = 16; n <= 65536; n *= 16) {
        vector<vector<double>> input_vectors = generateRandomVectors(n, vectorSize);
        for (int j = 1; j <= 16; ++j) {
            omp_set_num_threads(j);
            double t1, t2, dt;
            t1 = omp_get_wtime ();
            double result = multi_scaliar(input_vectors);
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
