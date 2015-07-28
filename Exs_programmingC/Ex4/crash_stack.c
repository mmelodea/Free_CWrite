#include <stdio.h>

int recursive_factorial(int num) {
  printf("\nnum: %d",num);
  return num * recursive_factorial(num+1);
}

int main(int argc, const char * argv[]) {
  printf("Recursive infinity factorial",recursive_factorial(0));
  return 0;
}

