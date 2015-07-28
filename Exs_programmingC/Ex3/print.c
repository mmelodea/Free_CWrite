#include <stdio.h>
#include "lenght.h"
#include "weight.h"
#include "temperature.h"

int main(int argc, const char * argv[])
{
  printf("\n**** Lenght Conversions ****");
  centimeter_inche(5);
  inche_centimeter(8);
  meter_feet(35);
  feet_meter(56);
  kilometer_mile(4);
  mile_kilometer(7);
  
  printf("\n\n**** Weight Conversions ****");
  pound_kilogram(10);
  kilogram_pound(6);
  gram_ounce(8);
  ounce_gram(1);
  
  printf("\n\n**** Temperature Conversions ****");
  celsius_fahrenheit(37);
  fahrenheit_celsius(32);
  celsius_kelvin(100);
  kelvin_celsius(273);
  kelvin_fahrenheit(273);
  fahrenheit_kelvin(212);
  
  printf("\n\n");
  return 0;
}
