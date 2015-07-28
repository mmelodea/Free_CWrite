#include <stdio.h>

int recursive_fibonacci(const int n) {
   if ( n == 0 )
      return 0;
   else if ( n == 1 )
      return 1;
   else
      return ( recursive_fibonacci(n-1) + recursive_fibonacci(n-2) );
}

int iterative_fibonacci(const int n) {
  int i, Fnv[n+1]; Fnv[0] = 0; Fnv[1] = 1;
  for(i=2; i<=n; i++) Fnv[i] = Fnv[i-1] + Fnv[i-2];
  return Fnv[n];
}

int main(int argc, const char * argv[]) {
  int n;
  printf("What Fibonacci serie element? ");
  scanf("%d",&n);
  printf("Iterative Fibonacci of F(n=%d): %d\n", n, iterative_fibonacci(n));
  printf("Recursive Fibonacci of F(n=%d): %d\n", n, recursive_fibonacci(n));
  return 0;
}