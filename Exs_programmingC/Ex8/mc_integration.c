#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 2000

float func(float x)
{
  return 4*x-2*x*x;
}

int main(int arg, const int *argv[])
{
  int i,ok=0;
  float li=0.,lf=2.,func_max=0.,rand_number,rand_func;
  float x[N];
  printf("Integration of the parabola 4x-2xx between xi=0 and xf=2 via MC:\n");
  for(i=0; i<N; i++){
    if( func(i/1000.) > func_max ) func_max = func(i/1000.);
  }
  //printf("func_max: %f",func_max);
  srand(time(NULL));
  for(i=0; i<N; i++){
     rand_number = (rand() % 1000+1)/1000.;
     rand_func = (rand() % 1000+1)/1000.;
     //printf("rand_number: %f,  rand_func: %f\n",rand_number,rand_func);
     if( rand_func*func_max < func(rand_number) ) ok += 1;
  }
  printf("Integral: %.3f\n\n",(ok/(float)N)*(lf-li)*func_max);
  return 0;
}