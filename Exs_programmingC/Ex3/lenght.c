#include <stdio.h>
#include "lenght.h"

float centimeter_meter(float l)
{
  return l/100.;
}

float centimeter_inche(float l)
{
  return printf("\nConversion: %.2f cm = %.2f in",l,l*0.394);
}

float inche_centimeter(float l)
{
  return printf("\nConversion: %.2f in = %.2f cm",l,l/0.394);
}

float meter_feet(float l)
{
  return printf("\nConversion: %.2f m = %.2f ft",l,l*3.281);
}

float feet_meter(float l)
{
  return printf("\nConversion: %.2f ft = %.2f m",l,l/3.281);
}

float kilometer_mile(float l)
{
  return printf("\nConversion: %.2f km = %.2f mi",l,l*0.621);
}

float mile_kilometer(float l)
{
  return printf("\nConversion: %.2f mi = %.2f km",l,l/0.621);
}