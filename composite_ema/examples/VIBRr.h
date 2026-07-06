/*
  МИНИМАЛЬНЫЙ shim ТОЛЬКО ДЛЯ АВТОНОМНОЙ СБОРКИ ПРИМЕРА ema_models_demo.c.

  Файлы ../source_c/COMPOSIT.c и ../source_c/COMPOSITE_MODELS.c начинаются с
  #include "VIBRr.h" — в проекте Лерера это его собственный заголовок. Чтобы
  собрать демо отдельно (gcc), здесь объявлены только те имена, которые эти
  файлы используют: комплексные помощники divide_/mult_, диспетчер материалов
  libmat и глобальная длина волны f. Их определения даёт сам ema_models_demo.c.

  В реальный проект Лерера этот файл включать НЕ нужно.
*/
#ifndef VIBRR_H_DEMO_SHIM
#define VIBRR_H_DEMO_SHIM

#include <stdio.h>
#include <math.h>

extern double f;   /* длина волны, нм (используется обёртками COMPOSITE/...) */

void divide_(double x, double y, double x1, double y1, double *u, double *v);
void mult_  (double x, double y, double x1, double y1, double *u, double *v);
void libmat (int n_mat, double lam, double *Re_eps, double *Im_eps);

#endif
