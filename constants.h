//
// Created by aksol on 16.06.2017.
//

#ifndef DIPLOMSOLDATOV_CONSTANTS_H
#define DIPLOMSOLDATOV_CONSTANTS_H

/**
 * Шаг по пространству
 */
const double_t STEP_X = 0.05;
/**
 * Коэффициент вытянутости конуса по z (в работе 2:1)
 * Задает шаг по y и форму прямоугольних треугольников (и углы)
 */
const double_t RATIO_Y_TO_X = 2;
/**
 * Шаг по времени
 */
const double_t TAU = 0.001;


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
const double_t TRIANGLE_HEIGHT = 2.0;
/**
 * Основание прямоугольного треугольника
 */
const double_t TRIANGLE_BASE = 1.0;
/**
 * Точность пересечения границы прямоугольного треугольника
 * (используется для контроля выхода за границу области при триангуляции)
 */
const double_t EPS_T = STEP_X * 0.001;
/**
 * Число узлов в выбранном симплексе
 * (В нашем случае это треугольник)
 */
const uint_fast32_t DIMENSION = 3;


/**
 * Начальное распределение температур
 */
const double_t INITIAL_TEMPERATURE = 300;
/**
 * Коэффициент температуропроводности титана (от 6,2 до 9,3 - проверить)
 */
const double_t Thermal_Diffusivity = 0.0000093;
/**
 * Коэффициент теплопроводности титана
 */
const double_t Thermal_Conductivity = 22.3;
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
const double_t TIME_WHEN_FUEL_ENDS = 39.134;

/**
 * Температура плавления (в задаче температура плавления титана)
 */
const double_t MELTING_TEMPERATURE = 1941.15;

const double_t PI = 3.1415926535897932384626433832795;

/*Debug */
const uint_fast32_t DEBUG_TIME_STEPS = 100;

/* End debug */
#endif //DIPLOMSOLDATOV_CONSTANTS_H
