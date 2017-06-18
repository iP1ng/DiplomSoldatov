
#ifndef DIPLOMSOLDATOV_CONSTANTS_H
#define DIPLOMSOLDATOV_CONSTANTS_H

/************************************************************************************/
/*
 * Сетка
 */
/************************************************************************************/
/**
 * Шаг по пространству
 */
const double_t STEP_X = 0.05;
/**
 * Шаг по времени
 */
const double_t TAU = 0.01;


/************************************************************************************/
/*
 * Треугольник
 *
 *         B
 *        /|
 *       / | TRIANGLE_HEIGHT
 *      /__|
 *     A    C
 *  TRIANGLE_BASE
 */
/************************************************************************************/
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
 * Высота треугольника
 */
const double_t TRIANGLE_HEIGHT = 2.0;
/**
 * Основание треугольника
 */
const double_t TRIANGLE_BASE = 1.0;
/**
 * Коэффициент вытянутости треугольника по высоте (в работе 2:1)
 * Задает шаг по y и форму прямоугольних треугольников (и углы)
 */
const double_t RATIO_Y_TO_X = 2;
/**
 * Точность пересечения границы прямоугольного треугольника
 * (используется для контроля выхода за границу области при триангуляции)
 */
const double_t EPS_T = STEP_X * 0.001;
/**
 * Число сторон треугольника
 */
const uint_fast32_t DIMENSION = 3;


/************************************************************************************/
/*
 * Коэффициенты уравнения теплопроводности
 */
/************************************************************************************/
/**
 * Начальное распределение температур (К)
 */
const double_t INITIAL_TEMPERATURE = 300;
/**
 * Коэффициент температуропроводности (м^2/c)
 * (в работе это титан)
 * Взят с http://ispu.ru/files/u2/SP._bez_nomera_-_Spravochn._materialy_dlya_resheniya_zadach_po_kursu_Teplomassoobmen..pdf
 */
const double_t Thermal_Diffusivity = 0.00000622;
/**
 * Коэффициент теплопроводности (Вт/(м*К))
 * (в работе это титан)
 * Взят с http://ispu.ru/files/u2/SP._bez_nomera_-_Spravochn._materialy_dlya_resheniya_zadach_po_kursu_Teplomassoobmen..pdf
 */
const double_t Thermal_Conductivity = 15.1;
/*
 * Коэффициент сопротивления формы конуса 2:1 (острием к потоку)
 * (безразмерная величина, определяющая реакцию среды на движение в ней тела)
 * Взят с https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D1%8D%D1%84%D1%84%D0%B8%D1%86%D0%B8%D0%B5%D0%BD%D1%82_%D1%81%D0%BE%D0%BF%D1%80%D0%BE%D1%82%D0%B8%D0%B2%D0%BB%D0%B5%D0%BD%D0%B8%D1%8F_%D1%84%D0%BE%D1%80%D0%BC%D1%8B
 */
const double_t Coeff_S = 0.5;
/**
 * Начальная плотность воздуха (кг/м^3)
 */
const double_t RHO_0 = 1.2;
/**
 * Коэффициент плотности воздуха (м^-1)
 */
const double_t Coeff_K = 0.00013;


/************************************************************************************/
/*
 * Коэффициенты полета ракеты
 */
/************************************************************************************/
/**
 * Ускорение ракеты (м/с^2)
 */
const double_t ACCELERATION = 90;


/************************************************************************************/
/*
 * Коэффициенты
 */
/************************************************************************************/
/**
 * Температура плавления
 * (в задаче температура плавления титана)
 * Взято с https://ru.wikipedia.org/wiki/%D0%A2%D0%B8%D1%82%D0%B0%D0%BD_(%D1%8D%D0%BB%D0%B5%D0%BC%D0%B5%D0%BD%D1%82)
 * Ракета летит до того момента, пока средняя температура на наконечнике ниже температуры плавления (по постановке)
 */
const double_t MELTING_TEMPERATURE = 1941.15;
/**
 * Число Пи
 */
const double_t PI = 3.1415926535897932384626433832795;
#endif //DIPLOMSOLDATOV_CONSTANTS_H
