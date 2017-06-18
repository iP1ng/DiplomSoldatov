//
// Created by aksol on 16.06.2017.
//

#ifndef DIPLOMSOLDATOV_CONSTANTS_H
#define DIPLOMSOLDATOV_CONSTANTS_H

/**
 * Шаг по пространству
 */
const double STEP_X = 0.05;
/**
 * Шаг по времени
 */
const double TAU = 0.001;


/**
 *         B
 *        /|
 *       / | TRIANGLE_HEIGHT
 *      /__|
 *     A    C
 *  TRIANGLE_BASE
 */
/**
 * Координаты треугольника
 */
const double_t A_x = 0;
const double_t A_y = 0;
const double_t B_x = 1;
const double_t B_y = 0;
const double_t C_x = 1;
const double_t C_y = 2;
/**
 * Высота прямоугольного треугольника
 */
const double TRIANGLE_HEIGHT = 2.0;
/**
 * Основание прямоугольного треугольника
 */
const double TRIANGLE_BASE = 1.0;
/**
 * Точность пересечения границы прямоугольного треугольника
 * (используется для контроля выхода за границу области при триангуляции)
 */
const double EPS_T = 0.0001;
/**
 * Число узлов в выбранном симплексе
 * (В нашем случае это треугольник)
 */
const uint_fast32_t DIMENSION = 3;


/**
 * Начальное распределение температур
 */
const double INITIAL_TEMPERATURE = 300;
/**
 * Коэффициент температуропроводности титана 9300000
 */
const double Thermal_Diffusivity = 9.3;
/**
 * Коэффициент теплопроводности титана
 */
const double Thermal_Conductivity = 22.3;
/* Коэффициент сопротивления формы конуса 2:1*/
const double_t Coeff_S = 0.5;

/**
 * Ускорение ракеты
 */
const double_t ACCELERATION = 90;

/**
 * начальная плотность
 */
const double_t RHO_0 = 1.2;

/**
 * Коэффициент плотности k
 */
const double_t Coeff_K = 0.00013;
/**
 * Время, когда кончится топливо
 */
const double TIME_WHEN_FUEL_ENDS = 39.134;

const double PI = 3.1415926535897932384626433832795;

/*Debug */
const uint_fast32_t DEBUG_TIME_STEPS = 1;

/* End debug */
#endif //DIPLOMSOLDATOV_CONSTANTS_H
