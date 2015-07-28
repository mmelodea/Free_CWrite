#include <stdio.h>
#include "temperature.h"

float celsius_fahrenheit(float t)
{
  return printf("\nConversion: %.1f°C = %.1f°F",t,180*t/100.+32);
}

float fahrenheit_celsius(float t)
{
  return printf("\nConversion: %.1f°F = %.1f°C",t,10/18.*(t-32));
}

float celsius_kelvin(float t)
{
  return printf("\nConversion: %.1f°C = %.1fK",t,t+273);
}

float kelvin_celsius(float t)
{
  return printf("\nConversion: %.1fK = %.1f°C",t,t-273); 
}

float kelvin_fahrenheit(float t)
{
  return printf("\nConversion: %.1fK = %.1f°F",t,(9/5.)*(t-273)+32);
}

float fahrenheit_kelvin(float t)
{
  return printf("\nConversion: %.1f°F = %.1fK",t,(5/9.)*(t-32)+273);
}