#include "header.h"

int main (int argc, char *argv[]){
	// Middle Rectangle Method
	CalcIntegralFloat(-2, 4, 3000, 1);
	CalcIntegralPosit(-2, 4, 3000, 1);
	CalcIntegralQuire(-2, 4, 3000, 1);
	
	// Ringe-Kutta Method
	RungeKutta2Float(1, exp(2.0), 2 * exp(2.0), 2, 100);
	printf("\n\n\n");
	RungeKutta2Posit(1, exp(2.0), 2 * exp(2.0), 2, 100);
	printf("\n\n\n");
	RungeKutta2Quire(1, exp(2.0), 2 * exp(2.0), 2, 100);
	
	// Jacobi Method
	float matrix[100][101];
	matrix[0][0] = 6;
	matrix[0][1] = 5;
	matrix[0][2] = 10;
	matrix[1][0] = -5;
	matrix[1][1] = 6;
	matrix[1][2] = 12;
	
	for(int i = 1000; i >= 1; i--){
		JacobiCalcFloat(matrix, 2, powf(10, -5) * i);
		printf("\n");
	}
	
	for(int i = 1000; i >= 1; i--){
		JacobiCalcPosit(matrix, 2, powf(10, -5) * i);
		printf("\n");
	}
	
	for(int i = 1000; i >= 1; i--){
		JacobiCalcQuire(matrix, 2, powf(10, -5) * i);
		printf("\n");
	}

	getchar();    
    return 0;
}

// gcc -lm -o main main.c RungeKutta2.c MiddleRectange.c Jacobi.c ../SoftPosit/build/Linux-x86_64-GCC/softposit.a  -I ../SoftPosit/source/include -O2 
