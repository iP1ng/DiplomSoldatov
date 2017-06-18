
#ifndef DIPLOMSOLDATOV_FUNCTIONS_H
#define DIPLOMSOLDATOV_FUNCTIONS_H

#include <stdlib.h>     /* qsort */


double_t func_calculate_rho(double_t t) {
    return (RHO_0 * exp(- Coeff_K * (ACCELERATION * t * t / 2)));
}

double_t func_calculate_q(double_t t) {
    double_t rho = func_calculate_rho(t);
    double_t V = ACCELERATION * t;

    return (Coeff_S * rho * V * V * V) /
            (2 * Thermal_Conductivity * PI * TRIANGLE_BASE * sqrt(TRIANGLE_BASE * TRIANGLE_BASE + TRIANGLE_HEIGHT * TRIANGLE_HEIGHT));
}

void func_multiply_matrix_and_vector(double_t* result_vector, double_t** matrix, double_t* vect, uint_fast32_t dim) {
    for (uint_fast32_t ii = 0; ii < DIMENSION; ii++)
        result_vector[ii] = 0;
    for (int ix = 0; ix < DIMENSION; ix++) {
        result_vector[ix] = 0;
        for (int jx = 0; jx < DIMENSION; jx++)
            result_vector[ix] += matrix[ix][jx] * vect[jx];
    }
}

void func_substract_two_matrices(double_t ** result_matrix, double_t** matrix1, double_t** matrix2, uint_fast32_t dim, double_t devide_element) {
    for (auto i = 0; i < dim; i++)
        for (auto j = 0; j < dim; j++)
            result_matrix[i][j] = (devide_element * matrix1[i][j]) - matrix2[i][j];

}

void print(vector< vector<double_t> > A) {
    ofstream file_slau("SLAU.dat", std::ofstream::out | std::ofstream::app);
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            file_slau << A[i][j] << "\t";
            if (j == n-1) {
                file_slau << "| ";
            }
        }
        file_slau << "\n";
    }
    cout << endl;
    file_slau.close();
}

#endif //DIPLOMSOLDATOV_FUNCTIONS_H
