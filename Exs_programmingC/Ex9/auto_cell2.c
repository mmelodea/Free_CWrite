#include <stdio.h>
#include <stdlib.h>

#define T 100
#define N 121

/*int rules(int a, int b, int c)
{
   int ruleset[] = {0,1,0,1,1,0,1,0};
   if      (a == 1 && b == 1 && c == 1) return ruleset[0];
   else if (a == 1 && b == 1 && c == 0) return ruleset[1];
   else if (a == 1 && b == 0 && c == 1) return ruleset[2];
   else if (a == 1 && b == 0 && c == 0) return ruleset[3];
   else if (a == 0 && b == 1 && c == 1) return ruleset[4];
   else if (a == 0 && b == 1 && c == 0) return ruleset[5];
   else if (a == 0 && b == 0 && c == 1) return ruleset[6];
   else if (a == 0 && b == 0 && c == 0) return ruleset[7];
   return 0;
}
*/
int rules(int a, int b, int c)
{
  int ruleset[] = {0,1};
  float weight = (a+b+c)/2.;
  if(weight == 0.) return ruleset[1];
  if(weight == 0.5) return ruleset[0];
  if(weight == 1.) return ruleset[0];
  if(weight == 1.5) return ruleset[1];
  return 0;
}

int main()
{
  int g,h,i,left, middle, right,cell;
  int cells[N], next_generation[N];
  
  srand(time(NULL));
  for(h=0; h<N; h++)
  {
    cells[h] = rand() % 2;
  }

/*  for(h=0; h<N; h++)
  {
    if(h==60) cells[h] = 1;
    else cells[h] = 0;
  }*/

  printf("\t******** 1D Cellular Automata ********\n\t******** Ploting generations\n\n");
  for(g=0; g<T; g++)//Time
  {
    for(i=0; i<N; i++)
    {
      //printf("%d",cells[i]);
      if(i==0)
      {
	left = cells[N-1];
	middle = cells[0];
	right = cells[1];
	//next_generation[0] = cells[0];
      }
      if(i==N-1)
      {
	left = cells[N-2];
	middle = cells[N-1];
	right = cells[0];
	//next_generation[N-1] = cells[N-1];
      }
      else
      {
	left = cells[i-1];
	middle = cells[i];
	right = cells[i+1];
	//next_generation[i] = rules(left,middle,right);
      }
      next_generation[i] = rules(left,middle,right);
    }
    for(i=0; i<N; i++)
    {
      if(cells[i]==0) printf(" ");
      if(cells[i]==1) printf("H");
      cells[i] = next_generation[i];  
    }
    printf("      Time g= %d\n",g);
  }
  
  return 0;
}