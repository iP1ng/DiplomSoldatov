
#ifndef TRIANGLE_H
#define TRIANGLE_H

/* Подключаем быстрые int фиксированного размера */
#include <cstdint>
#include <vector>

/* Подключаем double_t тип для оптимальной работы на любой архитектуре компьютера. */
#include <math.h>

/* Мои хедеры */
#include "../constants.h"
#include "../structures.h"

/* Debug */
#include <iostream>
#include <fstream>

using std::vector;
using std::cout;
using std::endl;

/*
 * Построение треугольной сетки на равнобедренном треугольнике.
 * Левая координата треугольника должна совпадать с началом координат.
 */
class IsoscelesTriangleGrid {

private:
    /*
     * m_a - левая точка треугольника
     * m_b - правая точка треугольника
     * m_c - центральная точка треугольника
     * m_hx - шаг по оси x
     */
    points m_a;
    points m_b;
    points m_c;
    double_t m_hx;
    triangles triangle;

public:
    vector<triangles> triangles_array;

    /*
     * Конструктор класса IsoscelesTriangleGrid
     */
    IsoscelesTriangleGrid(points a, points b, points c, double_t hx)
    {
        m_a = a;
        m_b = b;
        m_c = c;
        m_hx = hx;
    }

    double_t GetStep() { return m_hx; }

    /*
     * Метод построения сетки
     * TODO массив задавать не хардкодом, а вычислять число треугольников заранее
     * На выход число узлов
     */
    uint_fast32_t GetGreed(bool debug);

    /*
     * Функция, вычисляющая значение координаты y по заданной координате x на прямой AB треугольника
     * или на прямой, сдвинутой относительно AB на step
     * (структура треугольника описана в constants.h)
     */
    double_t LineFunction_ab(double_t x, double_t step);

    /*
     * Функция, вычисляющая значение координаты y на прямой BC треугольника
     * (структура треугольника описана в constants.h)
     */
    double_t LineFunction_bc();
};

#endif //TRIANGLE_H
