
#include "Triangle/triangle.h"
#include "Gauss/gauss.h"
#include "functions.h"
#include <iomanip> // setprecicion for cout Debug

using namespace std;

int main() {
    /************************************************************************************/
    /*
     * Треугольник и триангуляция
     */
    /************************************************************************************/
    /**
     * Отладочные сообщения для функци построения триангуляции и нумерации узлов
     */
    bool debug = false;
    /**
     * Создаем треугольник с координатами (A_x, A_y), (B_x, B_y), (C_x, C_y)
     */
    IsoscelesTriangleGrid triangle((points){A_x, A_y}, (points){B_x, B_y}, (points){C_x, C_y});
    /**
     * Разбиваем треугольник на маленькие треугольники с шагом STEP_X
     * Число получееных узлов возвращаем в переменную 
     * +1 т.к. нумерация узлов в методе с нуля
     */
    uint_fast32_t n = triangle.GetGreed(STEP_X, debug) +1;
    /**
     * Число маленьких треугольников, полученых после триангуляции
     */
    uint_fast32_t k = triangle.triangles_array.size();
    
    /************************************************************************************/
    /*
     * Уравнение теплопроводности, МКЭ
     */
    /************************************************************************************/
    /**
     * Тепловой поток, заданный на границе
     */
    double_t heat_flow = 0;
    
    /**
     * Вектор F правых частей для элемента
     */
    double_t *F_elem = new double_t[DIMENSION];
    /**
     * Вектор F правой части для всех элементов
     */
    double_t *F = new double_t [n];
    /**
     * Вектор F правой части для всех элементов на предыдущем временном слое
     */
    double_t *F_old = new double_t[n];
    /**
     * 2/tau [C] - [K]
     */
    double_t **R = new double_t *[n];    
    /**
     * Вектор правых частей итоговой системы (F + R)
     */
    double_t *Result = new double_t[n];
    
    /**
     * Итоговая матрица теплопроводности
     */
    double_t **Main_Matrix= new double_t *[n];
    /**
     * Итоговая матрица теплопроводности, объединенная с итоговым вектором правой части
     */
    vector< vector<double_t> > Result_matrix(n, std::vector<double_t>(n+1));    
    
    /**
     * Распределение температуры для жлемента
     */
    double_t Temperature_elem[DIMENSION];    
    /**
     * Итоговое распределение температуры во всем треугольнике
     */
    double_t *Temperature = new double_t[n];
    
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
     * к которому нужно будет прибавить полученные правые части F_elem на элементе
     */
    uint_fast32_t index[DIMENSION];

    /**
     * Результат вычитания двух матриц. Вспомогательная переменная
     */
    double_t ** substracted_matrix = new double_t* [DIMENSION];
    /**
     * Получаем температуру на новом слое, решая уравнение ниже методом Гаусса
     * Main_Matrix * Temp_0 = Result = > Temp_1,
     * где Result_matrix = Main_Matrix | Result
     */
    Gauss equation;

    /**
     * Средняя температура на наконечнике 
     */
    double_t average_temp = 0;
    /* 
     * Время, которое уже успела пролететь ракета 
     */
    double_t flight_time = - TAU;
    

    /************************************************************************************/
    /*
     * Debug
     */
    /************************************************************************************/
    cout << "Number of dots = " << n << endl;
    cout << "Number of triangles = " << k << endl;
    cout << "Heat flow heat_flow = " << heat_flow << endl;
    cout << "Space step h = " << STEP_X << endl;
    cout << "Time step tau = " << TAU << endl;
    cout << endl;


    /************************************************************************************/
    /*
     * Различные файлы для визуализации в Wolfram Mathematica
     */
    /************************************************************************************/
    /**
     * Debug вывод распределения температуры на треугольнике во времени
     * для проверки полученного решения
     */
    ofstream triangle_temp("triangle_temp.dat", std::ofstream::out);
    uint_fast32_t dots_num = floor(1 / STEP_X) + 1;
    uint_fast32_t dots_num_next = dots_num;
    uint_fast32_t dots_shift = 0;
    /**
     * Координаты треугольников, полученных триангуляцией
     */
    ofstream triangles_coordinates("triangles_coordinates.dat", std::ofstream::out);
    /**
     * Изменение теплового потока во времени
     */
    ofstream q_difference("q_difference.dat", std::ofstream::out);
    /**
     * Изменение плотности воздуха во времени
     */
    ofstream rho_difference("rho_difference.dat", std::ofstream::out);
    /**
     * Коэффициенты функции формы для элементов
     */
    ofstream form_func_coeffs("form_func_coeffs.dat", std::ofstream::out);
    
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

    form_func_coeffs << k * 3 << endl;

    triangle_temp << "Number of dots = " << n << endl;
    triangle_temp << "Number of triangles = " << k << endl;
    triangle_temp << "Space step h = " << STEP_X << endl;
    triangle_temp << "Time step tau = " << TAU << endl;
    triangle_temp << "Thermal_Diffusivity = " << Thermal_Diffusivity << endl;
    triangle_temp << "Thermal_Conductivity = " << Thermal_Conductivity << endl;


    /************************************************************************************/
    /*
     * Ицициализация переменных
     */
    /************************************************************************************/

    for (auto i = 0; i < n; i++){
        F[i] = 0;
        R[i] = new double_t[n];
        Result[i] = 0;
        Main_Matrix[i] = new double_t[n];
        Temperature[i] = INITIAL_TEMPERATURE;
        for (auto j = 0; j < n; j++) {
            R[i][j] = 0;
            Main_Matrix[i][j] = 0;
        }
    }

    for (auto i = 0; i < DIMENSION; i++) {
        substracted_matrix[i] = new double_t[DIMENSION];
        K[i] = new double_t[DIMENSION];
        C[i] = new double_t[DIMENSION];
    }


    /************************************************************************************/
    /*
     * Ищем решения по временным слоям
     * Пока не начнется плавление наконечника
     */
    /************************************************************************************/
    while (average_temp < MELTING_TEMPERATURE) {

        average_temp = 0;
        flight_time += TAU;

        for(auto i = 0; i < n; i++)
            average_temp += Temperature[i] / n;

        heat_flow = -func_calculate_q(flight_time);

        rho_difference << flight_time << ' ' << func_calculate_rho(flight_time) << endl;
        q_difference << flight_time << ' '  << - heat_flow << endl;

        /**
         * Обнуляем матрицы коэффициентов и вектор правой части элементов
         * перед заходом на новый временной слой
         */
        for (auto i = 0; i < DIMENSION; i++) {
            F_elem[i] = 0;
            Temperature_elem[i] = 0;
            for (auto j = 0; j < DIMENSION; j++) {
                K[i][j] = 0;
                C[i][j] = 0;
                substracted_matrix[i][j] = 0;
            }
        }

        /**
         * Обнуление матрицы жесткости и вектора правой части итоговой системы
         * перед заходом на новый временной слой
         * Также запоминаем значение вектора правой части с предыдущего временного слоя
         */
        for (auto i = 0; i < n; i++) {
            F_old[i] = F[i];
            F[i] = 0;
            Result[i] = 0;
            for (int j = 0; j < n; j++) {
                R[i][j] = 0;
                Main_Matrix[i][j] = 0;
            }
        }

        /************************************************************************************/
        /*
         * Идем по всем элементам
         */
        /************************************************************************************/
        for (auto element = 0; element < k; element++) {
            /*
             * Вычисляем матрицы коэффициентов K, C и вектор правых частей F_elem для элемента
             */
            triangle.triangles_array[element].Matrix_K(K);
            triangle.triangles_array[element].Matrix_C(C);
            triangle.triangles_array[element].Column_F(F_elem, heat_flow);

            /**
             * Запоминаем индексы глобальной матрицы и глобального вектора,
             * с которыми потом будем объединять значения, полученные на элементе
             */
            index[0] = triangle.triangles_array[element].first_point.point_num;
            index[1] = triangle.triangles_array[element].second_point.point_num;
            index[2] = triangle.triangles_array[element].third_point.point_num;

            /**
             * Записываем температуры на текущем временном слое тех узлов,
             * которые входят в элемент
             */
            for (auto i = 0; i < DIMENSION; i++) {
                Temperature_elem[i] = Temperature[index[i]];
            }

            /**
             * Вычисляем разность матриц
             * substracted_matrix = 2.0/TAU [C] - [K]
             */
            func_substract_two_matrices(substracted_matrix, C, K, 3, 2.0/TAU);

            /*
             * Заполняем итоговую матрицу жесткости и вектор правой части  полученными на элементе значениями
             * Реализуем схему 11.23 стр. 206 Сегерлинд
             *
             * F = -2F*
             * Main_Matrix = [K] + 2/tau [C]
             * R = 2/tau [C] - [K]
             */
            for (auto i = 0; i < DIMENSION; i++) {
                /* по формуле 11.21 */
                F[index[i]] += - (F_elem[i]);
                for (auto j = 0; j < DIMENSION; j++) {
                    Result_matrix[index[i]][index[j]] += (C[i][j] * 2 / TAU) + K[i][j];
                    R[index[i]][index[j]] += substracted_matrix[i][j];
                }
            }
        } // Конец цикла по элементам

        /*
         * Заполняыем  итоговый вектор правой части уравнения 11.23
         * Он состоит из поэлементно накопленных значений R + F
         */
        for (auto i = 0; i < n; i++) {
            for (auto j = 0; j <n; j++)
                Result_matrix[i][n] += (R[i][j] * Temperature[j]) + (F[i] - F_old[i]);
        }

        /**
         * Вычисляем температуру на новом слое по времени
         */
        equation.Solve(Temperature, Result_matrix);

        /************************************************************************************/
        /*
         * Debug
         */
        /************************************************************************************/
        dots_num = floor(1 / STEP_X) + 1;
        dots_num_next = dots_num;
        dots_shift = 0;
        triangle_temp << endl;
        triangle_temp << "Flight time = " << flight_time << endl;
        triangle_temp << "Heat flow heat_flow = " << heat_flow << endl;
        triangle_temp << "Average temperature = " << average_temp << endl;
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
            triangle_temp << setw(10) << setprecision(8)<< Temperature[d];
        }

        /*
         * Вычисляем коэффициенты функции формы
         */
        /*double_t N[3];

        for (int j = 0; j<k; j++ ){

            index[0] = triangle.triangles_array[j].first_point.point_num;
            index[1] = triangle.triangles_array[j].second_point.point_num;
            index[2] = triangle.triangles_array[j].third_point.point_num;


            for (uint_fast32_t jj = 0; jj < DIMENSION; jj++)
                Temperature_elem[jj] = Temperature[index[jj]];

            triangle.triangles_array[j].GetN(N,Temperature_elem);
            for (auto pp = 0; pp < DIMENSION; pp++)
                form_func_coeffs << N[pp] << ' ' << endl;
        }*/

    } // Конец цикла по времени

    triangle_temp.close();
    triangles_coordinates.close();
    rho_difference.close();
    q_difference.close();
    form_func_coeffs.close();
    return 0;
}
