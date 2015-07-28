#include <stdio.h>
int main ()
{
  unsigned int x = 0x76543210;
  char *c = (char*) &x;

  printf("Testing x=0x%x\n",x);
  printf("*c is: 0x%x\n", *c);
  if(*c == 0x10) printf("Underlying architecture is little endian. \n");
  else if (*c == 0x76) printf("Underlying architecture is big endian. \n");
 
  return 0;
}