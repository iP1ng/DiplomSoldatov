
#include "Triangle/triangle.h"
#include "Gauss/gauss.h"
#include <iomanip> // setprecicion for cout Debug


double_t func_calculate_rho(double_t t);

double_t func_calculate_q(double_t t);

void func_multiply_matrix_and_vector(double_t* result_vector, double_t** matrix, double_t* vect, uint_fast32_t dim);

void func_substract_two_matrices(double_t** result_matrix, double_t** matrix1, double_t** matrix2, uint_fast32_t dim, double_t devide_element);

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
    double_t q = 0;
    /**
     * Столбец правых частей для элемента на текущем слое
     */
    double_t *F = new double_t[DIMENSION];
    /**
     * Столбец правых частей для элемента
     */
    double_t *F_old = new double_t[n];
    /**
     * Значения температуры в узлах
     */
    double_t *Temp = new double_t[n];
    /**
     * Cюда записываем те температуры с предыдыщуего слоя которые нас интересуют для данного конечного элемента
     */
    double_t Fi[3];
    /**
     * Промежуточная переменная для записи произведения матрицы и столбца
     */
    double_t *Resultic = new double_t[DIMENSION];
    /**
     * Вектор правых частей итоговой системы
     */
    double_t *Result = new double_t[n];
    /**
     * Итоговая матрица жесткости системы
     */
    double_t **Main_Matrix= new double_t *[n];
    /**
     * Матрица коэффициентов элемента K
     */
    double_t **K = new double_t *[DIMENSION];
    /**
     * Матрица коэффициентов элемента C
     */
    double_t **C = new double_t *[DIMENSION];
    /**
     * Массив номеров строк в итоговом векторе правых частей Result,
     * к которому нужно будет прибавить полученные правые части F на элементе
     */
    uint_fast32_t ind[DIMENSION];
    /**
     * Результат вычитания двух матриц. Вспомогательная переменная
     */
    double_t ** substracted_matrix = new double_t* [DIMENSION];
    /**
     * Еще один Гаусс
     */
    vector< vector<double_t> > equation1_A(n, std::vector<double_t>(n+1));
    Gauss equation;
    /**
     * Вспомогательные переменные для накапливания матриц поэлементно
     */
    double_t **R1; R1 = new double_t *[n];
    double_t *FF; FF  = new double_t [n];
    for (int kk = 0; kk<n; kk++){
        R1[kk] = new double_t[n];
        FF[kk] = 0;
        for (int j = 0; j< n; j++)
            R1[kk][j] = 0;
    }

    /**
     * Wolfram Mathematica
     */
    ofstream triangles_coordinates("triangles_coordinates.dat", std::ofstream::out);

    triangles_coordinates << n << endl;
    for (auto i = 0; i < k; i ++) {
            triangles_coordinates << triangle.triangles_array[i].first_point.x << ' ';
            triangles_coordinates << triangle.triangles_array[i].first_point.y << ' ';
            triangles_coordinates << endl;
            triangles_coordinates << triangle.triangles_array[i].second_point.x << ' ';
            triangles_coordinates << triangle.triangles_array[i].second_point.y << ' ';
            triangles_coordinates << endl;
            triangles_coordinates << triangle.triangles_array[i].third_point.x << ' ';
            triangles_coordinates << triangle.triangles_array[i].third_point.y << ' ';
            triangles_coordinates << endl;
        }
    triangles_coordinates.close();

    ofstream q_difference("q_difference.dat", std::ofstream::out);
    ofstream rho_difference("rho_difference.dat", std::ofstream::out);
    ofstream form_fuc_coeffs("form_fuc_coeffs.dat", std::ofstream::out);
    form_fuc_coeffs << k * 3 << endl;

    ofstream triangle_temp("triangle_temp.dat", std::ofstream::out);
    triangle_temp << "Number of dots = " << n << endl;
    triangle_temp << "Number of triangles = " << k << endl;
    triangle_temp << "Space step h = " << STEP_X << endl;
    triangle_temp << "Time step tau = " << TAU << endl;
    triangle_temp << "Thermal_Diffusivity = " << Thermal_Diffusivity << endl;
    triangle_temp << "Thermal_Conductivity = " << Thermal_Conductivity << endl;



    /* Debug */
    cout << "Number of dots = " << n << endl;
    cout << "Number of triangles = " << k << endl;
    cout << "Heat flow q = " << q << endl;
    cout << "Space step h = " << STEP_X << endl;
    cout << "Time step tau = " << TAU << endl;
    cout << endl;
    /* End debug */

    /* Задаем начальное условие распределения температуры */
    for (uint_fast32_t i = 0; i < n; i++) {
        Temp[i] = INITIAL_TEMPERATURE;
    }

    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        substracted_matrix[i] = new double_t[DIMENSION];
    }

    /* Инициализация матрц K, C и вектора правых частей F для элемента */
    for (uint_fast32_t i = 0; i < DIMENSION; i++) {
        K[i] = new double_t[DIMENSION];
        C[i] = new double_t[DIMENSION];
    }

    /* Инициализация матрицы жесткости и вектора правой части итоговой системы */
    for (uint_fast32_t i = 0; i < n; i++) {
        Main_Matrix[i] = new double_t[n];
        Result[i] = 0;
        for (int j = 0; j < n; j++) {
            Main_Matrix[i][j] = 0;
        }
    }

    /************************************************************************************/
    /*
     * Ищем решения по временным слоям
     * Пока не начнется плавление наконечника
     */
    /************************************************************************************/
    uint_fast32_t number_of_time_steps = floor(TIME_WHEN_FUEL_ENDS / TAU);
    /* Средняя температура на наконечнике */
    double_t average_temp = 0;
    /* Время, которое уже успела пролететь ракета */
    double_t flight_time = - TAU;
    //double_t rho_for_file = 0;
    //for (uint_fast32_t global_tau =0; global_tau  < number_of_time_steps; global_tau++) {
    //for (uint_fast32_t global_tau =0; global_tau  < DEBUG_TIME_STEPS; global_tau++) {
    while (average_temp < MELTING_TEMPERATURE) {
        average_temp = 0;
        flight_time += TAU;

        for(uint_fast32_t i = 0; i < n; i++)
            average_temp += Temp[i] / n;

        cout << "flight_time = " << flight_time << endl;
        cout << "average_temp = " << average_temp << endl;

        //rho_for_file = func_calculate_rho(global_tau * TAU);
        q = -func_calculate_q(flight_time);
        //rho_difference << global_tau * TAU << ' ' << rho_for_file << endl;
        //q_difference << global_tau * TAU << ' '  << - q << endl;


        //cout << "q = " << q << endl;
        //q = - 10000;

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

        /*
         * Обнуление матрицы жесткости и вектора правой части итоговой системы
         */
        for (uint_fast32_t i = 0; i < n; i++) {
            Result[i] = 0;
            for (int j = 0; j < n; j++) {
                Main_Matrix[i][j] = 0;
            }
        }

        /*
         * Обнуляем вспомогательные переменные для поэлементного накопления
         */
        for (int kk = 0; kk<n; kk++){
            F_old[kk] = FF[kk];
            FF[kk] = 0;
            for (int j = 0; j< n; j++)
                R1[kk][j] = 0;
        }

        /************************************************************************************/
        /*
         * Идем по всем элементам
         */
        /************************************************************************************/
        for (uint_fast32_t i = 0; i < k; i++) {
            /* Вычисляем матрицы коэффициентов K, C и вектор правых частей F для элемента k */
            triangle.triangles_array[i].Matrix_K(K);
            triangle.triangles_array[i].Matrix_C(C);
            triangle.triangles_array[i].Column_F(F, q);

            ind[0] = triangle.triangles_array[i].first_point.point_num;
            ind[1] = triangle.triangles_array[i].second_point.point_num;
            ind[2] = triangle.triangles_array[i].third_point.point_num;

            /*
             * В Fi записываем температуры Temp на текущем временном слое тех узлов, которые входят в элемент k
             */
            for (uint_fast32_t j = 0; j < DIMENSION; j++) {
                Fi[j] = Temp[ind[j]];
            }

            func_substract_two_matrices(substracted_matrix, C, K, 3, 2.0/TAU);

            func_multiply_matrix_and_vector(Resultic, substracted_matrix, Fi, 3);

            /*
             * Заполняем итоговую матрицу жесткости и вектор правой части  полученными на элементе значениями
             * Реализуем схему 11.23 стр. 206 Сегерлинд
             *
             * FF = -2F*
             *
             * Main_Matrix = [K] + 2/tau [C]
             *
             * R1 = 2/tau [C] - [K]
             *
             */
            for (uint_fast32_t j = 0; j < 3; j++) {
                // по формуле 11.21
                FF[ind[j]] += - (F[j]);
                for (uint_fast32_t l = 0; l < 3; l++) {
                    Main_Matrix[ind[j]][ind[l]] += (C[j][l] * 2 / TAU) + K[j][l];
                    R1[ind[j]][ind[l]] += substracted_matrix[j][l];
                }
            }
        }

        /*
         * Заполняыем Result - итоговый вектор правой части уравнения 11.23
         * Он состоит из поэлементно накопленных значений
         * Result = R1 + FF
         */
        for (int ii = 0; ii < n; ii++) {
            Result[ii] = 0;
            for (int jj = 0; jj <n; jj++)
                Result[ii] += R1[ii][jj] * Temp[jj];
        }
        for (int ii = 0; ii<n; ii++){
            // На этом моменте в Result содержится полностью сформированный вектор правой части
            Result[ii]+= (FF[ii] - F_old[ii]);
            //cout<<setprecision(9)<<Result[ii]<<endl;
        }

        //TODO заменить Main_Matrix везде на equation1_A
        for (int p=0; p<n; p++) {
            for (int j=0; j<n; j++) {
                equation1_A[p][j] = Main_Matrix[p][j];
            }
        }
        for (int p=0; p<n; p++) {
            equation1_A[p][n] = Result[p];
        }

        /*
         * Вычисляем температуру Temp на новом слое по времени
         */
        equation.Solve(Temp, equation1_A);

        /* Deebug */
        uint_fast32_t dots_num = floor(1 / STEP_X) + 1;
        uint_fast32_t dots_num_next = dots_num;
        uint_fast32_t dots_shift = 0;

        triangle_temp << endl;
        triangle_temp << "Temperature on time = " << flight_time << endl;
        triangle_temp << "with heat flow q = " << q << endl;
        for (auto strs = 0; strs < dots_num; strs++)
            triangle_temp << "**********";
        triangle_temp << endl;
        for (auto d = 0; d < n; d++) {
            if (d == dots_num) {
                dots_num_next = dots_num_next - 1;
                dots_num = dots_num + dots_num_next;
                dots_shift ++;
                triangle_temp << endl;
                for (auto st = 0; st < dots_shift; st++)
                    triangle_temp << setw(10) << " ";
            }
            triangle_temp << setw(10) << setprecision(8)<< Temp[d];
        }
        /* End debug */
        /*
         * Вычисляем коэффициенты функции формы
         */

        /*double_t N[3];

        for (int j = 0; j<k; j++ ){

            ind[0] = triangle.triangles_array[j].first_point.point_num;
            ind[1] = triangle.triangles_array[j].second_point.point_num;
            ind[2] = triangle.triangles_array[j].third_point.point_num;


            for (uint_fast32_t jj = 0; jj < DIMENSION; jj++)
                Fi[jj] = Temp[ind[jj]];

            triangle.triangles_array[j].GetN(N,Fi);
            for (auto pp = 0; pp < DIMENSION; pp++)
                form_fuc_coeffs << N[pp] << ' ' << endl;
        }*/

    } // Конец цикла по времени

    q_difference.close();
    form_fuc_coeffs.close();
    return 0;
}

double_t func_calculate_rho(double_t t)
{
    return (RHO_0 * exp(- Coeff_K * (ACCELERATION * t * t / 2)));
}

double_t func_calculate_q(double_t t)
{

    double_t rho = func_calculate_rho(t);
    double_t V = 0;

    if (t < TIME_WHEN_FUEL_ENDS)
        V = ACCELERATION * t;

    return ((Coeff_S * rho * V * V * V) /
            (2 * Thermal_Conductivity * PI * TRIANGLE_BASE * sqrt(TRIANGLE_BASE * TRIANGLE_BASE + TRIANGLE_HEIGHT * TRIANGLE_HEIGHT)));
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