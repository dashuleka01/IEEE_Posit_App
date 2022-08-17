#include "header.h"

float GetXFloat (float x0, int i, float h) {
	return x0 + i * h;
}

void FunctionFloat (float x, float *y, float *z) {
	z[0] = y[1];
	z[1] = (2 * x * z[0] - 2 * y[0]) / (x * x - 1);
}

void Function1Float (float x, float *y, float *z) {
	z[0] = y[1];
	z[1] = (x + 1) / x * z[0] + 2 * (x - 1) / x * y[0];
}

void RungeKutta2Float(float x0, float y0, float y_der, float x1, int n) {
	float h = (x1 - x0) / n;
	float k1[2], k2[2], k3[2], k4[2], cur_y[2];
	float x[n + 1], y[n + 1];

	x[0] = x0;
	y[0] = y0;
	cur_y[0] = y0;
	cur_y[1] = y_der;
	
	for(int i = 0; i < n; i++) {
		Function1Float(x[i], cur_y, k1);
		
		float cur_y_k2[] = {cur_y[0] + h * k1[0] / 2, cur_y[1] + h * k1[1] / 2};
		Function1Float(x[i] + h / 2, cur_y_k2, k2);
		
		float cur_y_k3[] = {cur_y[0] + h * k2[0] / 2, cur_y[1] + h * k2[1] / 2};
		Function1Float(x[i] + h / 2, cur_y_k3, k3);
		
		float cur_y_k4[] = {cur_y[0] + h * k3[0], cur_y[1] + h * k3[1]};
		Function1Float(x[i] + h, cur_y_k4, k4);
		
		cur_y[0] += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
		cur_y[1] += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
		
		y[i + 1] = cur_y[0];
		x[i + 1] = GetXFloat(x0, i + 1, h);
	}
	
	printf("Runge-Kutta Method with IEEE-754:\n");
	for(int i = 0; i <= n; i++) {
		printf("%i: x = %.20f, y = %.20f\n", i, x[i], y[i]);
	}
}



void FunctionPosit(posit32_t x, posit32_t *y, posit32_t *z) {
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const1 = convertDoubleToP32(1);
	z[0] = y[1];
	z[1] = p32_div(p32_sub(p32_mul(const2, p32_mul(x, z[0])), 
		p32_mul(const2, y[0])), p32_sub(p32_mul(x, x), const1));
}

void Function1Posit(posit32_t x, posit32_t *y, posit32_t *z) {
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const1 = convertDoubleToP32(1);
	z[0] = y[1];
	
	z[1] = p32_add(p32_mul(p32_div(p32_add(x, const1), x), z[0]),
		p32_mul(p32_div(p32_sub(x, const1), x), p32_mul(const2, y[0])));
}


posit32_t GetXPosit(posit32_t x0, int i, posit32_t h) {
	return p32_add(x0, p32_mul(i32_to_p32(i), h));
}

void RungeKutta2Posit(float x0, float y0, float y_der, float x1, int n) {
	posit32_t x0_p = convertDoubleToP32(x0);
	posit32_t x1_p = convertDoubleToP32(x1);
	posit32_t h = p32_div(p32_sub(x1_p, x0_p), convertDoubleToP32(n));
	
	posit32_t k1[2], k2[2], k3[2], k4[2], cur_y[2];
	posit32_t x[n + 1], y[n + 1];

	x[0] = x0_p;
	y[0] = convertDoubleToP32(y0);
	cur_y[0] = y[0];
	cur_y[1] = convertDoubleToP32(y_der);
	
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const6 = convertDoubleToP32(6);
	
	for(int i = 0; i < n; i++) {
		Function1Posit(x[i], cur_y, k1);
		
		posit32_t cur_y_k2[] = {
			p32_add(cur_y[0], p32_div(p32_mul(k1[0], h), const2)),
			p32_add(cur_y[1], p32_div(p32_mul(k1[1], h), const2))
		};
		Function1Posit(p32_add(x[i], p32_div(h, const2)), cur_y_k2, k2);
		
		
		posit32_t cur_y_k3[] = {
			p32_add(cur_y[0], p32_div(p32_mul(k2[0], h), const2)),
			p32_add(cur_y[1], p32_div(p32_mul(k2[1], h), const2))
		};
		Function1Posit(p32_add(x[i], p32_div(h, const2)), cur_y_k3, k3);
		
		posit32_t cur_y_k4[] = {
			p32_add(cur_y[0], p32_mul(k3[0], h)),
			p32_add(cur_y[1], p32_mul(k3[1], h))
		};
		Function1Posit(p32_add(x[i], h), cur_y_k4, k4);
		
		cur_y[0] = p32_add(cur_y[0], p32_mul(p32_div(h, const6),
			p32_add(p32_add(k1[0], p32_mul(const2, k2[0])),
			p32_add(p32_mul(const2, k3[0]), k4[0]))));
		
		cur_y[1] = p32_add(cur_y[1], p32_mul(p32_div(h, const6),
			p32_add(p32_add(k1[1], p32_mul(const2, k2[1])),
			p32_add(p32_mul(const2, k3[1]), k4[1]))));
		
		y[i + 1] = cur_y[0];
		x[i + 1] = GetXPosit(x0_p, i + 1, h);
	}
	
	printf("Runge-Kutta Method with Posit:\n");
	for(int i = 0; i <= n; i++) {
		printf("%i: x = %.20f, y = %.20f\n", i, 
			convertP32ToDouble(x[i]), convertP32ToDouble(y[i]));
	}
}



void FunctionQuire(posit32_t x, quire32_t *y, quire32_t *z) {
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const1 = convertDoubleToP32(1);
	
	z[0] = y[1];
	
	posit32_t temp = p32_div(const1, p32_sub(p32_mul(x, x), const1));
	quire32_t zero = q32_clr(zero);
	z[1] = q32_clr(z[1]);
	z[1] = q32_fdp_add(z[1], p32_mul(const2, x), q32_to_p32(z[0]));
	z[1] = q32_fdp_sub(z[1], const2, q32_to_p32(y[0]));
	z[1] = q32_fdp_add(zero, q32_to_p32(z[1]), temp);
}

void Function1Quire(posit32_t x, quire32_t *y, quire32_t *z) {
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const1 = convertDoubleToP32(1);
	
	z[0] = y[1];
	
	quire32_t zero = q32_clr(zero);
	z[1] = q32_clr(z[1]);
	z[1] = q32_fdp_add(z[1], p32_div(p32_add(x, const1), x), q32_to_p32(z[0]));
	z[1] = q32_fdp_add(z[1], p32_div(p32_sub(x, const1), x), 
		p32_mul(const2, q32_to_p32(y[0])));
}

posit32_t GetXQuire(posit32_t x0, int i, posit32_t h) {
	quire32_t result = q32_clr(result);
	result = q32_fdp_add(result, x0, i32_to_p32(1));
	result = q32_fdp_add(result, i32_to_p32(i), h);
	return q32_to_p32(result);
}

void RungeKutta2Quire(float x0, float y0, float y_der, float x1, int n) {
	posit32_t const1 = convertDoubleToP32(1);
	posit32_t const2 = convertDoubleToP32(2);
	posit32_t const6 = convertDoubleToP32(6);
	posit32_t const3 = convertDoubleToP32(3);
	
	posit32_t x0_p = convertDoubleToP32(x0);
	posit32_t x1_p = convertDoubleToP32(x1);
	posit32_t h = p32_div(p32_sub(x1_p, x0_p), convertDoubleToP32(n));
	
	quire32_t k1[2], k2[2], k3[2], k4[2], cur_y[2];
	posit32_t x[n + 1], y[n + 1];

	x[0] = x0_p;
	y[0] = convertDoubleToP32(y0);
	cur_y[0]= q32_clr(cur_y[0]);
	cur_y[1]= q32_clr(cur_y[1]);
	cur_y[0] = q32_fdp_add(cur_y[0], const1, y[0]);
	cur_y[1] = q32_fdp_add(cur_y[1], const1, convertDoubleToP32(y_der));
	
	for(int i = 0; i < n; i++) {
		Function1Quire(x[i], cur_y, k1);
		
		quire32_t cur_y_k2[] = {
			q32_fdp_add(cur_y[0], q32_to_p32(k1[0]), p32_div(h, const2)),
			q32_fdp_add(cur_y[1],  q32_to_p32(k1[1]), p32_div(h, const2))
		};
		Function1Quire(p32_add(x[i], p32_div(h, const2)), cur_y_k2, k2);
		
		quire32_t cur_y_k3[] = {
			q32_fdp_add(cur_y[0], q32_to_p32(k2[0]), p32_div(h, const2)),
			q32_fdp_add(cur_y[1], q32_to_p32(k2[1]), p32_div(h, const2))
		};
		Function1Quire(p32_add(x[i], p32_div(h, const2)), cur_y_k3, k3);
		
		quire32_t cur_y_k4[] = {
			q32_fdp_add(cur_y[0], q32_to_p32(k3[0]), h),
			q32_fdp_add(cur_y[1], q32_to_p32(k3[1]), h)
		};
		Function1Quire(p32_add(x[i], h), cur_y_k4, k4);
		
		cur_y[0] = q32_fdp_add(cur_y[0], p32_div(h, const6), q32_to_p32(k1[0]));
		cur_y[0] = q32_fdp_add(cur_y[0], p32_div(h, const3), q32_to_p32(k2[0]));
		cur_y[0] = q32_fdp_add(cur_y[0], p32_div(h, const3), q32_to_p32(k3[0]));
		cur_y[0] = q32_fdp_add(cur_y[0], p32_div(h, const6), q32_to_p32(k4[0]));
		
		cur_y[1] = q32_fdp_add(cur_y[1], p32_div(h, const6), q32_to_p32(k1[1]));
		cur_y[1] = q32_fdp_add(cur_y[1], p32_div(h, const3), q32_to_p32(k2[1]));
		cur_y[1] = q32_fdp_add(cur_y[1], p32_div(h, const3), q32_to_p32(k3[1]));
		cur_y[1] = q32_fdp_add(cur_y[1], p32_div(h, const6), q32_to_p32(k4[1]));
		
		y[i + 1] = q32_to_p32(cur_y[0]);
		x[i + 1] = GetXQuire(x0_p, i + 1, h);
	}
	
	printf("Runge-Kutta Method with Posit + Quire:\n");
	for(int i = 0; i <= n; i++) {
		printf("%i: x = %.20f, y = %.20f\n", i, 
			convertP32ToDouble(x[i]), convertP32ToDouble(y[i]));
	}
}