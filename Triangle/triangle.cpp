#include "triangle.h"

double_t IsoscelesTriangleGrid::LineFunction_ab(double_t x, double_t step)
{
    double_t y = RATIO_Y_TO_X * (x - step);
    return y;
}

double_t IsoscelesTriangleGrid::LineFunction_bc()
{
    return TRIANGLE_BASE;
}

uint_fast32_t IsoscelesTriangleGrid::GetGreed(double_t m_hx, points *cordinates, bool debug)
{
    /*
     * Массив координат треугольников
     * TODO изменить фиксированную размерность. Заранее вычислять число треугольников.
     */

    /* Координаты узлов сетки */
    double_t x = m_a.x;
    double_t y = 0;

    /*
     * t1 - разница между координатой y правой точки маленького треугольника и правой границей
     * основного треугольника. Если разница меньше eps, считаем что мы дошли до границы.
     * t2 - то же самое, для перевернутых маленьких треугольников.
     * Если выполнены условия для t1 и t2, значит этот уровень триангулирован и надо идти вверх по оси y
     * на следующий уровень.
     */
    double_t eps = EPS_T;
    double_t t1 = 0;
    double_t t2 = 0;

    /*
     * Счетчики
     * i - счетчик прямых, параллельных левой грани треугольника (идет через два шага, одинаково на всех уровнях)
     * обнуляется в конце уровня
     * j - счетчик треугольников в массиве треугольников
     * ii - счетчик сдвига по x при переходе на ряд выше ( счетчик иксов на уровне)
     * l - счетчик рядов по высоте (он идет до середины треугольника)
     * k - нумератор узлов
     * n - число узлов на уровне (на сколько узлов делится горизонтальная сторона треугольника)
     * length - длина уровня
     */
    uint_fast32_t i = 0;
    uint_fast32_t j = 0;
    uint_fast32_t ii = 0;
    uint_fast32_t l = 0;

    uint_fast32_t k = 0;
    uint_fast32_t n = 0;
    // m_b.x -1 потому что тут бага в координатах (см файл constants).
    double_t length = LineFunction_bc() + STEP_X;

    uint_fast32_t  coord_num = 0;

    /* Переход с ряда на ряд (вверх), пока не дойдем до середины треугольника */
    while (y < TRIANGLE_HEIGHT) {
        length -= STEP_X;
        n += (length / m_hx) + 1;
        /* Цикл по горизонтальному ряду */
        for (;;) {
            x = ii * m_hx;
            y = LineFunction_ab(x, i * m_hx);
            cordinates[coord_num].x = x;
            cordinates[coord_num].y = y;
            coord_num ++;
            triangle.first_point = { x, y, k };

            x = (ii + 1) * m_hx;
            y = LineFunction_ab(x, (i + 1) * m_hx);

            triangle.second_point = { x, y, k + 1 };

            t1 = abs(x - LineFunction_bc());

            x = (ii + 1) * m_hx;
            y = LineFunction_ab(x, i * m_hx);

            triangle.third_point = { x, y, k + ((length / m_hx) + 1) };

            triangles_array.push_back(triangle);

            if (debug == true) {
                cout << "First point of triange " << j << ": " << endl;
                cout << "x: " << triangles_array.at(j).first_point.x << endl;
                cout << "y: " << triangles_array.at(j).first_point.y << endl;
                cout << "k: " << triangles_array.at(j).first_point.point_num << endl;
                cout << "Second point of triange " << j << ": " << endl;
                cout << "x: " << triangles_array.at(j).second_point.x << endl;
                cout << "y: " << triangles_array.at(j).second_point.y << endl;
                cout << "k: " << triangles_array.at(j).second_point.point_num << endl;
                cout << "Third point of triange " << j << ": " << endl;
                cout << "x: " << triangles_array.at(j).third_point.x << endl;
                cout << "y: " << triangles_array.at(j).third_point.y << endl;
                cout << "k: " << triangles_array.at(j).third_point.point_num << endl;
            }

            t2 = abs(x - LineFunction_bc());

            if ((t2 < EPS_T) && (t1 < EPS_T)) {
                ii = l + 1;
                i = 0;
                j++;
                k = n;
                cordinates[coord_num].x = triangle.second_point.x;
                cordinates[coord_num].y = triangle.second_point.y;
                coord_num ++;
                if (debug == true) {
                    cout << "BREAK CONDITION." << endl;
                    cout << "n = " << n << endl;
                    cout << "k = " << k << endl;
                }
                break;
            }

            triangle.first_point.x = triangles_array.back().third_point.x;
            triangle.first_point.y = triangles_array.back().third_point.y;
            triangle.first_point.point_num = triangles_array.back().third_point.point_num;
            triangle.second_point.x = triangles_array.back().second_point.x;
            triangle.second_point.y = triangles_array.back().second_point.y;
            triangle.second_point.point_num = triangles_array.back().second_point.point_num;

            x = (ii + 2) * m_hx;
            y = LineFunction_ab(x, (i + 1) * m_hx);

            triangle.third_point = { x, y, k + 1 + ((length / m_hx) + 1) };
            triangles_array.push_back(triangle);

            if (debug == true) {
                cout << "Second point of triange " << j + 1 << ": " << endl;
                cout << "x: " << triangles_array.at(j + 1).first_point.x << endl;
                cout << "y: " << triangles_array.at(j + 1).first_point.y << endl;
                cout << "k: " << triangles_array.at(j + 1).first_point.point_num << endl;
                cout << "Second point of triange " << j + 1 << ": " << endl;
                cout << "x: " << triangles_array.at(j + 1).second_point.x << endl;
                cout << "y: " << triangles_array.at(j + 1).second_point.y << endl;
                cout << "k: " << triangles_array.at(j + 1).second_point.point_num << endl;
                cout << "Third point of triange " << j + 1 << ": " << endl;
                cout << "x: " << triangles_array.at(j + 1).third_point.x << endl;
                cout << "y: " << triangles_array.at(j + 1).third_point.y << endl;
                cout << "k: " << triangles_array.at(j + 1).third_point.point_num << endl;
            }
            i = i + 1;
            ii = ii + 1;
            j = j + 2;
            k++;
        }
        l++;
        if (debug == true) {
            cout << "New level l = " << l << endl;
        }
    }
    cordinates[coord_num].x = C_x;
    cordinates[coord_num].y = C_y;
    //cout << "coord_num = " << coord_num << endl;
    return triangles_array.back().third_point.point_num;
}
