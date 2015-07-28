#include <stdio.h>
void swap(int* as, int* bs) {
int tmp = *as;
*as = *bs;
*bs = tmp;
printf("swap: a -> %d - b -> %d\n", *as, *bs);
}
int main(int argc, const char * argv[]) {
int a = 42;
int b = 21;
printf("a -> %d - b -> %d\n", a, b);
swap(&a, &b);
printf("a -> %d - b -> %d\n", a, b);
return 0;
}
