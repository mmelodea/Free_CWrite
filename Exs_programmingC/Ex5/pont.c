#include <stdio.h>

void func(int *sa, int *sb)
{
  if(*sa<0 || *sb<0) printf("Error! at least one argument negative");
  printf("Sum: %d\n",*sa+*sb);
  if(*sa > *sb) printf("Difference: %d\n",*sa-*sb);
  if(*sa < *sb) printf("Difference: %d\n",*sb-*sa);
  printf("Mean: %.2f\n",(*sa+*sb)/2.);
  //printf("Val: %d",sa); acessar o endereço da variavel na memoria
}

int main(int arg, const char *argv[])
{
  int a, b;
  printf("Put two numbers:\n");
  scanf("%d",&a);
  scanf("%d",&b);
  func(&a,&b);
  return 0;
}