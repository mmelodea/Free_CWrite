#include <stdio.h>
#include "weight.h"

float pound_kilogram(float w)
{
  return printf("\nConversion: %.2f lb = %.2f kg",w,w/2.205);
}

float kilogram_pound(float w)
{
  return printf("\nConversion: %.2f kg = %.2f lb",w,w*2.205);
}

float gram_ounce(float w)
{
  return printf("\nConversion: %.2f g = %.2f oz",w,w*0.035);
}

float ounce_gram(float w)
{
  return printf("\nConversion: %.2f oz = %.2f g",w,w/0.035);
}