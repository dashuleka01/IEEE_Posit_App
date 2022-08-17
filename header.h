#ifndef _header_h
#define _header_h

#include "softposit.h"
#include "internals.h"
#include <math.h>

void CalcIntegralFloat(float a, float b, int n, int step);
void CalcIntegralPosit(float a, float b, int n, int step);
void CalcIntegralQuire(float a, float b, int n, int step);

void JacobiCalcFloat(float matrix[100][101], int n, float e);
void JacobiCalcPosit(float matrix[100][101], int n, float e);
void JacobiCalcQuire(float matrix[100][101], int n, float e);

void RungeKutta2Float(float x0, float y0, float y_der, float x1, int n);
void RungeKutta2Posit(float x0, float y0, float y_der, float x1, int n);
void RungeKutta2Quire(float x0, float y0, float y_der, float x1, int n);

#endif