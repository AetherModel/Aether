// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

/*
 This is an example code that highlights how to use the armadillo code.
 This allows you to use loop-less math, so you can add and do math
 with the arrays.  It makes everything much more streamlined.

 We use 3D arrays alot in this code, and this library makes it so it
 is each to use these 3d arrys.

 Here is how I compiled the code (clearly, I installed it in /usr/local):
 g++ arm.cpp -o arm.exe -std=c++11 -O2 -larmadillo -I/usr/local/include -L/usr/local/lib

*/

#include <iostream>
#include <armadillo>
#include <cmath>

#include "../../include/sizes.h"
#include "../../include/constants.h"

#include "arm_vars.h"

// -----------------------------------------------------
// These are defined for the class
// -----------------------------------------------------

Grid::Grid(int nX, int nY, int nZ) {

  radius3d.set_size(nX, nY, nZ);
  radius3d.zeros();

  fcube tmp(nX, nY, nZ);
  tmp.zeros();
  
  // Create a 3-element vector of 3D variables:
  gravity3dv.push_back(tmp);
  gravity3dv.push_back(tmp);
  gravity3dv.push_back(tmp);
  
  return;
  
}

void Grid::set_radius(float planet_radius, fcube alts3d) {

  radius3d = planet_radius + alts3d;
  return;

}

fcube Grid::get_radius() {

  return radius3d;

}

// -----------------------------------------------------
// main code to demonstrate how armadillo works.
// -----------------------------------------------------


int main(int argc, char** argv) {

  // cubes are defined as row, column, slice
  // so, rows = lons
  //     cols = lats
  //     slices = const alt

  // fcube is float cube
  // fixed is not changable, and the size has to be constants
  //   - this makes it difficult to make it generalizable.

  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> lon3d;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> lat3d;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> alt3d;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> local_time;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> sza;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> cos_sza;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> temperature;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> radius;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> gravity;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> scale_height;
  //fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> density;

  long nLons = nGeoLonsG;
  long nLats = nGeoLatsG;
  long nAlts = nGeoAltsG;

  Grid grid(nLons, nLats, nAlts);  

  // Ok, it doesn't seem like it makes a difference whether we use a
  // static variable (nGeoAltsG) vs a dynamic variable (nAlts)
  // Timing is very similar!
  
  fcube lon3d(nLons, nLats, nAlts);
  fcube lat3d(nLons, nLats, nAlts);
  fcube alt3d(nLons, nLats, nAlts);
  fcube local_time(nLons, nLats, nAlts);
  fcube sza(nLons, nLats, nAlts);
  fcube cos_sza(nLons, nLats, nAlts);
  fcube radius(nLons, nLats, nAlts);
  fcube gravity(nLons, nLats, nAlts);

  float earth_radius = 6372000.0;
  float declination = 0.0;
  float sin_dec = sin(declination);
  float cos_dec = cos(declination);

  float ut = 0.0;
  
  //fvec::fixed<nLons> lon1d;
  //fvec::fixed<nLats> lat1d;
  //fvec::fixed<nAlts> alt1d;
  //fvec::fixed<nAlts> temp1d;

  fvec lon1d(nLons);
  fvec lat1d(nLats);
  fvec alt1d(nAlts);
  fvec temp1d(nAlts);

  std::cout << "Lons :\n";
  float dLon = twopi/nGeoLons;
  for (int i=0; i < nLons; i++)
    lon1d(i) = (i-nGeoGhosts+0.5) * dLon;
  (lon1d*rtod).print();

  for (int i=0; i < nLats; i++) {
    for (int j=0; j < nAlts; j++) {
      lon3d.subcube(0,i,j,nLons-1,i,j) = lon1d;
    }
  }
  
  std::cout << "Lats :\n";
  float dLat = pi/nGeoLats;
  for (int i=0; i < nLats; i++)
    lat1d(i) = -pi/2 + (i-nGeoGhosts+0.5) * dLat;
  (lat1d*rtod).print();
  
  for (int i=0; i < nLons; i++) {
    for (int j=0; j < nAlts; j++) {
      lat3d.subcube(i,0,j,i,nLats-1,j) = lat1d;
    }
  }
  
  std::cout << "Alts :\n";
  float dAlt = 5000.0;
  for (int i=0; i < nAlts; i++)
    alt1d(i) = 100000.0 + (i-nGeoGhosts) * dAlt;
  (alt1d/1000.0).print();

  for (int i=0; i < nLons; i++) {
    for (int j=0; j < nLats; j++) {
      alt3d.tube(i,j) = alt1d;
    }
  }

  // Let's see if we can do some manipulation, just to see what happens:

  // If we want to time this, just do it over and over again.
  // Commented it out for now:
  // for (long iTime = 0; iTime < 1000; iTime++) {

  // This just demonstrates a set and get functional pair:
  grid.set_radius(earth_radius, alt3d);
  radius = grid.get_radius();

  // radius = earth_radius + alt3d;

  gravity = 9.8 * (earth_radius / radius) % (earth_radius / radius);

  fvec tmp1d(nAlts);
  tmp1d = gravity.tube(2,2);
  tmp1d.print("Gravity as a function of altitude :");

  local_time = (lon3d*rtod/15.0 + ut)*pi/12.0;
  
  cos_sza =
    sin_dec * sin(lat3d) +
    cos_dec * cos(lat3d) % cos(local_time - pi);
  sza = acos(cos_sza);

  // --------------------------------------------------------------
  // Now, let's make some species.  
  // --------------------------------------------------------------
  
  struct species_chars {
    float mass;
    float density_bc;
    // It seems like you have to declare cubes like this for
    // a structure.
    fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> scale_height;
    fcube::fixed<nGeoLonsG, nGeoLatsG, nGeoAltsG> density;
    //fcube scale_height(nLons, nLats, nAlts);
    //fcube density(nLons, nLats, nAlts);
  };

  std::vector<species_chars> neutrals;
  species_chars tmp_species;
  int nSpecies = 0;
  
  fcube temperature(nLons, nLats, nAlts);
  //density.ones();

  // make up some species:
  // O:
  tmp_species.mass = 16.0 * amu;
  tmp_species.density_bc = 1.0e18;
  tmp_species.density.ones();
  tmp_species.scale_height.zeros();

  // Define O as #0:
  neutrals.push_back(tmp_species);
  nSpecies++;
  
  // O2:
  tmp_species.mass = 2 * 16.0 * amu;
  tmp_species.density_bc = 2.0e18;
  // Define O2 as #1:
  neutrals.push_back(tmp_species);
  nSpecies++;

  // N2:
  tmp_species.mass = 2 * 14.0 * amu;
  tmp_species.density_bc = 1.0e19;
  // Define N2 as #2:
  neutrals.push_back(tmp_species);
  nSpecies++;
  
  // Here we make some temperature dependence so that the dayside
  // is hotter than the nightside:
  float alt_min = alt3d.min();
  temperature = 200.0 +
    (700.0 + 100.0*cos_sza) % (1.0-exp(-(alt3d-alt_min)/(dAlt*15.0)));

  // Print out a strip near the equator:
  //int i = float(nLats)/2.0;
  //int j = nAlts-3;
  //std::cout << "Latitudes around the globe near equator : \n";
  //(lat3d.subcube(0, i, j, nLons-1, i, j)*rtod).print();
  //std::cout << "local times around the globe near equator : \n";
  //local_time.subcube(0, i, j, nLons-1, i, j).print();
  //std::cout << "cos(sza) around the globe near equator : \n";
  //cos_sza.subcube(0, i, j, nLons-1, i, j).print();
  //std::cout << "temperatures around the globe near equator : \n";
  //temperature.subcube(0, i, j, nLons-1, i, j).print();

  for (int iSpecies=0; iSpecies<nSpecies; iSpecies++) {

    neutrals[iSpecies].scale_height =
      boltzmanns_constant * temperature /
      (neutrals[iSpecies].mass * gravity);

    tmp1d = neutrals[iSpecies].scale_height.tube(2,2);
    (tmp1d/1000.0).print("scale height as a function of altitude (km) :");
  
    neutrals[iSpecies].density.slice(0).fill(neutrals[iSpecies].density_bc);
    for (int i=1; i < nAlts; i++) {
      neutrals[iSpecies].density.slice(i) = 
	neutrals[iSpecies].density.slice(i-1) %
	exp(-dAlt / neutrals[iSpecies].scale_height.slice(i-1));
    }
    tmp1d = neutrals[iSpecies].density.tube(2,2);
    tmp1d.print("density as a function of altitude :");
  }

  // This is for the timing loop...
  //}
  
  return 0;

}
