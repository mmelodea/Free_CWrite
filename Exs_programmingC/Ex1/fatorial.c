#include <stdio.h>

int recursive_factorial(int num) {
  if(num==1) return 1;
  return num * recursive_factorial(num-1);
}

int iterative_factorial(int num) {
  int i, result=1;
    for (i=1; i<num+1; i++) result *= i;
    return result;
}

int main(int argc, const char * argv[]) {
  printf("Iterative 6! :%d\n", iterative_factorial(6));
  printf("Recursive 6! :%d\n", recursive_factorial(6));
  return 0;
}

