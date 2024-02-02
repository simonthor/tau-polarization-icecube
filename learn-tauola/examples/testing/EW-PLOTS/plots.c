/**
 * EW plots generation.
 * Draw plots using 'draw.C' root script
 *
 * @author Tomasz Przedzinski
 * @date 28 January 2011
 */

#include "Tauola/Plots.h"
#include "Tauola/Tauola.h"

using namespace Tauolapp;

int main(int argc,char **argv){

  Tauola::initialize();

  Plots p;
  // Set incoming PDG ID and cosTheta for plots 1 and 2
  p.setSancVariables(2,-0.2);

  p.SANCtest1();
  p.SANCtest2();
  p.SANCtest3(); //longer
  p.SANCtest4(); //medium
}

