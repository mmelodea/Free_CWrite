#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int prime_test(int a)
{
  //printf("\nTestes: %i, %i, %i, %i",a%2,a%3,a%5,a%7);
  if( a==2 || a==3 || a==5 || a==7 ) return 0;
  if( (a%2)==0 ) return 1;
  if( (a%3)==0 ) return 1;
  if( (a%5)==0 ) return 1;
  if( (a%7)==0 ) return 1;
  else return 0;
}

int main(int arg, const int *argv[])
{
  int i=0, t=0, primes[100];
  for(i=0; i<100; i++) primes[i] = -1;
  for(i=2; i<=100; i++){
    //printf("\nPrime_test(%d): %i",i,prime_test(i));
    if( prime_test(i) == 0 ){ 
      //printf("\nprime_test(%d): %f",i,prime_test(i));
      primes[t] = i;
      t += 1;
    }
  }
  printf("\n\t\t The prime numbers between 1-100:\n");
  for(i=0; i<t; i++) printf("  %d",primes[i]);
  printf("\n\n");
  
  return 0;
}