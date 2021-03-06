//
// Created by aksol on 16.06.2017.
//

#ifndef DIPLOMSOLDATOV_STRUCTURES_H
#define DIPLOMSOLDATOV_STRUCTURES_H
#include <iostream>
#include "constants.h"

/*
 * Структура, описывающая координаты и номер узла, а также радиус-вектор узла.
 */
struct points {
    double_t x;
    double_t y;
    // TODO проверить, что радиус-вектор в цилиндрических координатах вычисляется так
    double_t rad_vector() { return x; }
    uint_fast32_t point_num;
};

/*
 * Структура, описывающая координаты точек и номер треугольника,
 * а также вычисляющая коэффициенты функции формы узлов, входящих в треугольник,
 * а также вычисляюая площадь треугольника и матрицы элементов [K], [C], [F].
 */
struct triangles {
    points first_point;
    points second_point;
    points third_point;

    void coef_a(double_t* a)
    {
        a[0] = second_point.x * third_point.y - third_point.x * second_point.y;
        a[1] = third_point.x * first_point.y - first_point.x * third_point.y;
        a[2] = first_point.x * second_point.y - second_point.x * first_point.y;
    }
    void coef_b(double_t* b)
    {
        b[0] = second_point.y - third_point.y;
        b[1] = third_point.y - first_point.y;
        b[2] = first_point.y - second_point.y;
    }
    void coef_c(double_t* c)
    {
        c[0] = third_point.x - second_point.x;
        c[1] = first_point.x - third_point.x;
        c[2] = second_point.x - first_point.x;
    }
    double_t GetSquareTriangleArea()
    {
        return STEP_X * (RATIO_Y_TO_X * STEP_X) * 0.5;
    }
    double_t GetMatrixADeterminant()
    {
        return 0.5 * (second_point.x * third_point.y
                      - third_point.x * second_point.y
                      - first_point.x * third_point.y
                      + first_point.x * second_point.y
                      + third_point.x * first_point.y
                      - second_point.x * first_point.y);
    }
    double_t GetR() {
        points p[3];

        p[0].x = first_point.x;
        p[0].y = first_point.y;
        p[1].x = second_point.x;
        p[1].y = second_point.y;
        p[2].x = third_point.x;
        p[2].y = third_point.y;

        double_t R = (0.0833333333333333333333333) * ((2 * p[0].rad_vector()
                                                       + p[1].rad_vector()
                                                       + p[2].rad_vector())
                                                      * p[0].rad_vector()
                                                      + (p[0].rad_vector()
                                                         + 2 * p[1].rad_vector()
                                                         + p[2].rad_vector())
                                                        * p[1].rad_vector()
                                                      + (p[0].rad_vector()
                                                         + p[1].rad_vector()
                                                         + 2 * p[2].rad_vector())
                                                        * p[2].rad_vector());
        return R;
    }

    /* Коэффициенты функции формы каждого симплекса */
    void GetN(double_t * N, double_t * Phi) {
        double_t a[3];
        double_t b[3];
        double_t c[3];
        coef_a(a);
        coef_b(b);
        coef_c(c);

        for (auto i = 0; i < DIMENSION; i++)
            N[i] = 0;

        for (auto i = 0; i < DIMENSION; i++) {
            N[0] += 1 / (2 * GetMatrixADeterminant()) * (a[i]) * Phi[i];
            N[1] += 1 / (2 * GetMatrixADeterminant()) * (b[i]) * Phi[i]; // при rho
            N[2] += 1 / (2 * GetMatrixADeterminant()) * (c[i]) * Phi[i]; // при z
        }
    }

    void Matrix_K(double_t** K)
    {
        double_t a[3];
        double_t b[3];
        double_t c[3];
        coef_a(a);
        coef_b(b);
        coef_c(c);

        double_t R = GetR();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {

                /* Для теста из учебника (стр. 95) подставить вместо Thermal_Diffusivity 40 */
                K[i][j] = (Thermal_Diffusivity * 2.0 * PI * R / (4.0 * this->GetMatrixADeterminant()))
                          * (b[i] * b[j] + c[i] * c[j]);
            }
        }
    }

    void Matrix_C(double_t** C)
    {
        points p[3];
        p[0].x = first_point.x;
        p[0].y = first_point.y;
        p[1].x = second_point.x;
        p[1].y = second_point.y;
        p[2].x = third_point.x;
        p[2].y = third_point.y;

        double_t D = 2 * PI * this->GetMatrixADeterminant() / 180;
        double_t R[3];
        for (int i = 0; i < 3; i++) {
            R[i] = p[i].rad_vector();
        }

        C[0][0] = D * (12 * R[0] * R[0]
                       + 2 * R[1] * R[1]
                       + 2 * R[2] * R[2]
                       + 6 * R[0] * R[1]
                       + 6 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[0][1] = D * (3 * R[0] * R[0]
                       + 3 * R[1] * R[1]
                       + R[2] * R[2]
                       + 4 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[0][2] = D * (3 * R[0] * R[0]
                       + 1 * R[1] * R[1]
                       + 3 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 4 * R[0] * R[2]
                       + 2 * R[1] * R[2]);
        C[1][0] = C[0][1];
        C[2][0] = C[0][2];
        C[1][1] = D * (2 * R[0] * R[0]
                       + 12 * R[1] * R[1]
                       + 2 * R[2] * R[2]
                       + 6 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 6 * R[1] * R[2]);
        C[1][2] = D * (1 * R[0] * R[0]
                       + 3 * R[1] * R[1]
                       + 3 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 2 * R[0] * R[2]
                       + 4 * R[1] * R[2]);
        C[2][2] = D * (2 * R[0] * R[0]
                       + 2 * R[1] * R[1]
                       + 12 * R[2] * R[2]
                       + 2 * R[0] * R[1]
                       + 6 * R[0] * R[2]
                       + 6 * R[1] * R[2]);
        C[2][1] = C[1][2];
    }

    void Column_F(double_t* F, double_t q)
    {
        double_t L = sqrt(STEP_X * STEP_X + (RATIO_Y_TO_X *STEP_X) * (RATIO_Y_TO_X * STEP_X));
        double_t k = L * q * 2.0 * PI / 6;

        F[0] = 0;
        F[1] = 0;
        F[2] = 0;
        /*
         * TODO Уточнить у Пугачева условие в левой нижней и верхней точках, где 
         * на точки действуют сразу два разных граничных условия
         * в этих местах температура сходит с ума (идет в минус)
         */
        if (fabs(RATIO_Y_TO_X * first_point.x - first_point.y) < EPS_T
            && fabs(RATIO_Y_TO_X * third_point.x - third_point.y) < EPS_T) {
            F[0] = k * (2.0 * first_point.rad_vector() + third_point.rad_vector());
            F[2] = k * (first_point.rad_vector() + 2.0 * third_point.rad_vector());

        }
    }

};
#endif //DIPLOMSOLDATOV_STRUCTURES_H
