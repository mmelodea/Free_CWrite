#include <stdio.h>
#include <math.h>
#define N 200

float func(float a){return 4*a-2*a*a;}

float intg(float x[], int n, float li, float lf)
{
  int id, res=0;
  for(id=1; id<n; id++) res += 2*func(x[id]);
  res += func(x[0])+func(x[n]);
  return 0.5*((lf-li)/(float)n)*res;
}

int main(int arg, const int *argv[])
{
  int i;
  float li=0,lf=2;
  float x[N+1];
  printf("Integration of the parabola 4x-2xx between xi=0 and xf=2:\n");
  for(i=0; i<=N; i++) x[i] = i/100.;
  //for(i=0; i<=N; i++) printf("x[%d]: %f\n",i,x[i]);
  printf("Integral: %.2f\n\n",intg(x,N,li,lf));
  return 0;
}