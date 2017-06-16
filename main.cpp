
#include "Triangle/triangle.h"
#include "Gauss/gauss.h"
#include <iomanip> // setprecicion for cout Debug


double func_calculate_rho(double t);

double func_calculate_q(double t);

void func_multiply_matrix_and_vector(double* result_vector, double** matrix, double* vect, uint_fast32_t dim);

void func_substract_two_matrices(double** result_matrix, double** matrix1, double** matrix2, uint_fast32_t dim, double devide_element);

void print(vector< vector<double> > A) {
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

using namespace std;

int main() {
    /**
     * Отладочные сообщения для функци построения триангуляции и нумерации узлов
     */
    bool debug = false;
    /**
     * Инстанс класса, отвечающий за триангулацию треугольника
     */
    IsoscelesTriangleGrid triangle((points){A_x, A_y}, (points){B_x, B_y}, (points){C_x, C_y}, STEP_X);
    /**
     * Число точек
     */
    uint_fast32_t n = triangle.GetGreed(debug) +1;
    /**
     * Число узлов
     */
    uint_fast32_t k = triangle.triangles_array.size();
    /**
     * Тепловой поток, заданный на границе
     */
    double q = 0;
    /**
     * Столбец правых частей для элемента
     */
    double *F = new double[DIMENSION];
    /**
     * Значения температуры в узлах
     */
    double *Temp = new double[n];
    /**
     * Cюда записываем те температуры с предыдыщуего слоя которые нас интересуют для данного конечного элемента
     */
    double Fi[3];
    /**
     * Промежуточная переменная для записи произведения матрицы и столбца
     */
    double *Resultic = new double[DIMENSION];
    /**
     * Вектор правых частей итоговой системы
     */
    double *Result = new double[n];
    /**
     * Итоговая матрица жесткости системы
     */
    double **Main_Matrix= new double *[n];
    /**
     * Матрица коэффициентов элемента K
     */
    double **K= new double *[DIMENSION];
    /**
     * Матрица коэффициентов элемента C
     */
    double **C = new double *[DIMENSION];
    /**
     * Массив номеров строк в итоговом векторе правых частей Result,
     * к которому нужно будет прибавить полученные правые части F на элементе
     */
    uint_fast32_t ind[DIMENSION];
    /**
     * Результат вычитания двух матриц. Вспомогательная переменная
     */
    double ** substracted_matrix = new double* [DIMENSION];
    /**
     * Еще один Гаусс
     */
    vector< vector<double> > equation1_A(n, std::vector<double>(n+1));
    Gauss equation;


    /* Debug */
    cout << "Number of dots = " << n << endl;
    cout << "Number of triangles = " << k << endl;
    cout << "Heat flux q = " << q << endl;
    cout << "Space step h = " << STEP_X << endl;
    cout << "Time step tau = " << TAU << endl;
    cout << endl;
    /* End debug */

    /* Задаем начальное условие распределения температуры */
    for (uint_fast32_t i = 0; i < n; i++) {
        Temp[i] = INITIAL_TEMPERATURE;
    }

    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        substracted_matrix[i] = new double[DIMENSION];
    }

    /* Инициализация матрц K, C и вектора правых частей F для элемента */
    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        K[i] = new double[DIMENSION];
        C[i] = new double[DIMENSION];
    }

    /* Ищем решения по временным слоям */
    uint_fast32_t number_of_time_steps = floor(TIME_WHEN_FUEL_ENDS / TAU);
    //for (uint_fast32_t global_tau =0; global_tau  < number_of_time_steps; global_tau++) {
    for (uint_fast32_t global_tau =0; global_tau  < 100; global_tau++) {

        q = -func_calculate_q(global_tau * TAU);
        cout << "q = " << q << endl;
        //q = 100000;

        for (uint_fast32_t i = 0; i < DIMENSION; i++) {
            F[i] = 0;
            Resultic[i] = 0;
            Fi[i] = 0;
            for (uint_fast32_t j = 0; j < DIMENSION; j++) {
                K[i][j] = 0;
                C[i][j] = 0;
                substracted_matrix[i][j] = 0;
            }
        }
        /* Debug */
        cout << "Time: "<< global_tau * TAU << endl;
        cout << "Temperature on step " << global_tau << endl;
        for (uint_fast32_t i = 0; i < n; i++) {
            cout <<setprecision(16)<< Temp[i] << "\t";
        }
        cout << endl;
        /* End debug */

        /* Инициализация матрицы жесткости и вектора правой части итоговой системы */
        for (uint_fast32_t i = 0; i < n; i++) {
            Main_Matrix[i] = new double[n];
            Result[i] = 0;
            for (int j = 0; j < n; j++) {
                Main_Matrix[i][j] = 0;
            }
        }


        double **R1; R1 = new double *[n];
        double *FF; FF  = new double [n];
        for (int kk = 0; kk<n; kk++){
            R1[kk] = new double[n];
            FF[kk] = 0;
            for (int j = 0; j< n; j++)
                R1[kk][j] = 0;
        }

        /* Идем по всем элементам */
        for (uint_fast32_t i = 0; i < k; i++) {
            /* Вычисляем матрицы коэффициентов K, C и вектор правых частей F для элемента k */
            triangle.triangles_array[i].Matrix_K(K);
            triangle.triangles_array[i].Matrix_C(C);
            triangle.triangles_array[i].Column_F(F, q);

            ind[0] = triangle.triangles_array[i].first_point.point_num;
            ind[1] = triangle.triangles_array[i].second_point.point_num;
            ind[2] = triangle.triangles_array[i].third_point.point_num;

            /* В Fi записываем температуры Temp на текущем временном слое тех узлов, которые входят в элемент k */
            for (uint_fast32_t j = 0; j < DIMENSION; j++) {
                Fi[j] = Temp[ind[j]];
            }

            func_substract_two_matrices(substracted_matrix, C, K, 3, 2.0/TAU);

            func_multiply_matrix_and_vector(Resultic, substracted_matrix, Fi, 3);

            /* Заполняем итоговую матрицу жесткости и вектор правой части  полученными на элементе значениями */
            for (uint_fast32_t j = 0; j < 3; j++) {
                //Result[ind[j]] += (Resultic[j]) - 2 * F[j];
                FF[ind[j]] += - 2 * F[j];
                //Result[ind[j]] += F[j];
                for (uint_fast32_t l = 0; l < 3; l++) {
                    Main_Matrix[ind[j]][ind[l]] += (C[j][l] * 2 / TAU) + K[j][l];
                    R1[ind[j]][ind[l]] += substracted_matrix[j][l];
                    //Main_Matrix[ind[j]][ind[l]] += K[j][l];
                }
            }
        }
        for (int ix = 0; ix < n; ix++) {
            Result[ix] = 0;
            for (int jx = 0; jx <n; jx++)
                Result[ix] += R1[ix][jx] * Temp[jx];
        }

        for (int ii = 0; ii<n; ii++){
            Result[ii]+= FF[ii];
            cout<<setprecision(9)<<Result[ii]<<endl;
        }

        for (int p=0; p<n; p++) {
            for (int j=0; j<n; j++) {
                equation1_A[p][j] = Main_Matrix[p][j];
            }
        }
        for (int p=0; p<n; p++) {
            equation1_A[p][n] = Result[p];
        }

        equation.Solve(Temp, equation1_A);
    }
    return 0;
}

double func_calculate_rho(double t)
{
    return (RHO_0 * exp(- Coeff_K * (ACCELERATION * t * t / 2)));
}

double func_calculate_q(double t)
{

    double rho = func_calculate_rho(t);
    double V = 0;

    if (t < TIME_WHEN_FUEL_ENDS)
        V = ACCELERATION * t;

    return ((Coeff_S * rho * V * V * V) /
            (Thermal_Conductivity * PI * TRIANGLE_BASE * sqrt(TRIANGLE_BASE * TRIANGLE_BASE + TRIANGLE_HEIGHT * TRIANGLE_HEIGHT)));
}

void func_multiply_matrix_and_vector(double* result_vector, double** matrix, double* vect, uint_fast32_t dim) {
    for (uint_fast32_t ii = 0; ii < DIMENSION; ii++)
        result_vector[ii] = 0;
    for (int ix = 0; ix < DIMENSION; ix++) {
        result_vector[ix] = 0;
        for (int jx = 0; jx < DIMENSION; jx++)
            result_vector[ix] += matrix[ix][jx] * vect[jx];
    }
}

void func_substract_two_matrices(double ** result_matrix, double** matrix1, double** matrix2, uint_fast32_t dim, double devide_element) {

    for (auto i = 0; i < dim; i++)
        for (auto j = 0; j < dim; j++)
            result_matrix[i][j] = (devide_element * matrix1[i][j]) - matrix2[i][j];

}