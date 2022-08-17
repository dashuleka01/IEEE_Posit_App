#include "header.h"

void JacobiCalcFloat(float matrix[100][101], int n, float e) {
	float new_matrix[n][n + 1], result[n], prev_result[n];
	float max_delta = 0;
	int iteration = 0;
	
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n + 1; j++) {
			new_matrix[i][j] = matrix[i][j] / matrix[i][i];
		}
		result[i] = new_matrix[i][n];
	}
	
	do {
		iteration++;
		for(int i = 0; i < n; i++) {
			prev_result[i] = result[i];
		}
		
		for(int i = 0; i < n; i++) {
			result[i] = new_matrix[i][n];
			for(int j = 0; j < n; j++) {
				if(i != j) {
					result[i] -= new_matrix[i][j] * prev_result[j];
				}
			}
		}
		
		max_delta = 0;
		for(int i = 0; i < n; i++) { 
			max_delta = fabs(prev_result[i] - result[i]) > max_delta 
				? fabs(prev_result[i] - result[i]) : max_delta;
		}
	} while(e <= max_delta);
	
	printf("Jacobi Method IEEE-754 with error %f:\n%i steps, ", e, iteration);
	for(int i = 0; i < n; i++) {
		printf("x%i = %.20f ", i, result[i]);
	}
}

posit32_t fabsPosit(posit32_t A) {
	union ui32_p32 uA;
	uA.p = A;
	if(uA.ui >> 31) {
		return p32_mul(A, convertDoubleToP32(-1));
	}
	return A;
}


void JacobiCalcPosit(float matrix[100][101], int n, float e) {
	posit32_t new_matrix[n][n + 1], result[n], prev_result[n];
	posit32_t max_delta = convertDoubleToP32(0);
	int iteration = 0;
	
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n + 1; j++) {
			new_matrix[i][j] = p32_div(
				convertDoubleToP32(matrix[i][j]), convertDoubleToP32(matrix[i][i]));
		}
		result[i] = new_matrix[i][n];
	}
	
	do {
		iteration++;
		for(int i = 0; i < n; i++) {
			prev_result[i] = result[i];
		}
		
		for(int i = 0; i < n; i++) {
			result[i] = new_matrix[i][n];
			for(int j = 0; j < n; j++) {
				if(i != j) {
					result[i] = p32_sub(result[i], 
						p32_mul(new_matrix[i][j], prev_result[j]));
				}
			}
		}
		
		max_delta = convertDoubleToP32(0);
		for(int i = 0; i < n; i++) { 
			posit32_t delta = fabsPosit(p32_sub(prev_result[i], result[i]));
			max_delta = p32_lt(max_delta, delta) ?
				delta : max_delta;
		}
	} while(p32_lt(convertDoubleToP32(e), max_delta));
	
	printf("Jacobi Method Posit with error %f:\n%i steps, ", e, iteration);
	for(int i = 0; i < n; i++) {
		printf("x%i = %.20f ", i, convertP32ToDouble(result[i]));
	}
}

void JacobiCalcQuire(float matrix[100][101], int n, float e) {
	quire32_t new_matrix[n][n + 1], result[n], prev_result[n];
	posit32_t max_delta = convertDoubleToP32(0);
	int iteration = 0;
	
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n + 1; j++) {
			new_matrix[i][j] = q32_clr(new_matrix[i][j]);
			
			new_matrix[i][j] = q32_fdp_add(new_matrix[i][j], 
				convertDoubleToP32(1), 
				p32_div(convertDoubleToP32(matrix[i][j]), 
					convertDoubleToP32(matrix[i][i]))); 
		}
		result[i] = new_matrix[i][n];
	}
	
	do {
		iteration++;
		for(int i = 0; i < n; i++) {
			prev_result[i] = result[i];
		}
		
		for(int i = 0; i < n; i++) {
			result[i] = new_matrix[i][n];
			for(int j = 0; j < n; j++) {
				if(i != j) {
					result[i] = q32_fdp_sub(result[i], 
						q32_to_p32(new_matrix[i][j]), q32_to_p32(prev_result[j]));
				}
			}
		}
		
		max_delta = convertDoubleToP32(0);
		for(int i = 0; i < n; i++) { 
			posit32_t delta = fabsPosit(p32_sub(
				q32_to_p32(prev_result[i]), q32_to_p32(result[i])));
			max_delta = p32_lt(max_delta, delta) ?
				delta : max_delta;
		}
	} while(p32_lt(convertDoubleToP32(e), max_delta));
	
	printf("Jacobi Method Posit + Quire with error %f:\n%i steps, ", e, iteration);
	for(int i = 0; i < n; i++) {
		printf("x%i = %.20f ", i, convertP32ToDouble(q32_to_p32(result[i])));
	}
}
