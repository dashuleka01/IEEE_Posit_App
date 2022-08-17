#include "header.h"

float InFunctionFloat(float x) {
	return 8 + 2 * x - x * x;
}

float InFunction1Float(float x) {
	return 1 / x;
}

void CalcIntegralFloat(float a, float b, int n, int step) {
	for(int i = step; i <= n; i += step) {
		float result = 0, h = (b - a) / i;

		for(int j = 0; j < i; j++) {
			result += InFunctionFloat(a + h * (j + 0.5));
		}
		result *= h;
		printf("Middle Rectange Method with IEEE-754: ");
		printf("%i elementary segments, result = %.20f\n", i, result);
	}
}


posit32_t InFunctionPosit(posit32_t x) {
	posit32_t const8 = convertDoubleToP32(8);
	posit32_t const2 = convertDoubleToP32(2);

	return p32_sub(p32_add(const8, p32_mul(const2, x)), p32_mul(x, x));
}

posit32_t InFunction1Posit(posit32_t x) {
	return p32_div(convertDoubleToP32(1), x);
}

void CalcIntegralPosit(float a, float b, int n, int step) {
	for(int i = step; i <= n; i += step) { 
		posit32_t result = convertDoubleToP32(0);
		posit32_t a_p = convertDoubleToP32(a);
		posit32_t b_p = convertDoubleToP32(b);
		posit32_t h = p32_div(p32_sub(b_p, a_p), convertDoubleToP32(i));
  
		for(int j = 0; j < i; j++) {
			posit32_t x = p32_add(a_p, 
				p32_mul(h, p32_add(convertDoubleToP32(j), convertDoubleToP32(0.5))));
			result = p32_add(result, InFunctionPosit(x));
		}
  
		result = p32_mul(result, h);
		printf("Middle Rectange Method with Posit: ");
		printf("%i elementary segments, result = %.20f\n", i, convertP32ToDouble(result));
	}
}


posit32_t InFunctionQuire(posit32_t x) {
	posit32_t const8 = convertDoubleToP32(8);
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const1 = convertDoubleToP32(1);
	
	quire32_t result = q32_clr(result);
	result = q32_fdp_add(result, const8, const1);
	result = q32_fdp_add(result, const2, x);
	result = q32_fdp_sub(result, x, x);

	return q32_to_p32(result);
}

posit32_t InFunction1Quire(posit32_t x) {
	return p32_div(convertDoubleToP32(1), x);
}

void CalcIntegralQuire(float a, float b, int n, int step) {
	for(int i = step; i <= n; i += step) { 
		quire32_t result = q32_clr(result);
		posit32_t a_p = convertDoubleToP32(a);
		posit32_t b_p = convertDoubleToP32(b);
		posit32_t h = p32_div(p32_sub(b_p, a_p), convertDoubleToP32(i));
  
		for(int j = 0; j < i; j++) {
			quire32_t x = q32_clr(x);
			x = q32_fdp_add(x, a_p, convertDoubleToP32(1));
			x = q32_fdp_add(x, h, convertDoubleToP32(j));
			x = q32_fdp_add(x, h, convertDoubleToP32(0.5));
			result = q32_fdp_add(result, 
				InFunctionQuire(q32_to_p32(x)), convertDoubleToP32(1));
		}
		printf("Middle Rectange Method with Posit + Quire: ");
		printf("%i elementary segments, result = %.20f\n", i, 
			convertP32ToDouble(p32_mul(q32_to_p32(result), h)));
	}
}
