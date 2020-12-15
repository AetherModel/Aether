// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include <iostream>

#include "../include/constants.h"
#include "../include/neutrals.h"
#include "../include/earth.h"
#include "../include/report.h"
#include "../include/time.h"
#include "../include/solvers.h"

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

void Neutrals::calc_mass_density(Report &report) {

  long iLon, iLat, iAlt, index, iSpecies;

  std::string function="Neutrals::calc_mass_density";
  report.enter(function);

  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	density_s3gc[index] = 0.0;
	rho_s3gc[index] = 0.0;

	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  density_s3gc[index] = density_s3gc[index] +
	    neutrals[iSpecies].density_s3gc[index];
	  rho_s3gc[index] = rho_s3gc[index] +
	    neutrals[iSpecies].mass *
	    neutrals[iSpecies].density_s3gc[index];
	}

	mean_major_mass_s3gc[index] = rho_s3gc[index] / density_s3gc[index];
	pressure_s3gc[index] =
	  density_s3gc[index] *
	  boltzmanns_constant *
	  temperature_s3gc[index];
      }
    }
  }
    
  report.exit(function);
  return;
  
}
  
//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

void Neutrals::calc_specific_heat(Report &report) {

  long iLon, iLat, iAlt, index, iSpecies;

  double t, p, r;
  
  std::string function="Neutrals::calc_specific_heat";
  report.enter(function);

  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	Cv_s3gc[index] = 0.0;
	gamma_s3gc[index] = 0.0;
	kappa_s3gc[index] = 0.0;

	t = temperature_s3gc[index];
		
	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  p = neutrals[iSpecies].thermal_exp;
	  r = pow(t,p); // temp ^ exp; <--- This line slows the code way down.
	  Cv_s3gc[index] = Cv_s3gc[index] +
	    (neutrals[iSpecies].vibe - 2) *
	    neutrals[iSpecies].density_s3gc[index] *
	    boltzmanns_constant / neutrals[iSpecies].mass;
	  gamma_s3gc[index] = gamma_s3gc[index] +
	    neutrals[iSpecies].density_s3gc[index] / (neutrals[iSpecies].vibe-2);
	  kappa_s3gc[index] = kappa_s3gc[index] +
	    neutrals[iSpecies].density_s3gc[index] *
	    neutrals[iSpecies].thermal_cond *
	    r;
	}

	Cv_s3gc[index] = Cv_s3gc[index] / (2*density_s3gc[index]);
	gamma_s3gc[index] = gamma_s3gc[index] * 2.0 / density_s3gc[index] + 1.0;
	kappa_s3gc[index] = kappa_s3gc[index] / density_s3gc[index];

	r = gamma_s3gc[index] *
	  boltzmanns_constant *
	  temperature_s3gc[index] /
	  mean_major_mass_s3gc[index];
	
	sound_s3gc[index] = sqrt(r); // <---- Same with this line

      }
    }
  }

  report.exit(function);
  return;
  
}

//----------------------------------------------------------------------
// Calculate the altitude integral of the different species for EUV
// calculations.  Scale these to be the slant path through the
// atmosphere.  This gets complicated behind the terminator.  All of
// this is taken from Smith and Smith, JGR 1972, vol. 77, page 3592
// ----------------------------------------------------------------------

void Neutrals::calc_chapman(Grid grid, Report &report) {

  long iAlt, iLon, iLat, index, indexp;
  float H;  // scale height

  // This is all from Smith and Smith, JGR 1972, vol. 77, page 3592
  // "Numerical evaluation of chapman's grazing incidence integral ch(X,x)"
  // Xp is supposed to be R/H
  // JMB Update: 05/2017.  Corrected a small error in the y-calc for 
  // erfc(y)
  //
  // Also Updated the Grazing Integral for SZA > 90.0
  // We now do log-linear interpolation for smoother transitions
  
  double a = 1.06069630;
  double b = 0.55643831;
  double c = 1.06198960;
  double d = 1.72456090;
  double f = 0.56498823;
  double g = 0.06651874;

  double integral[nGeoAltsG], xp[nGeoAltsG], erfcy[nGeoAltsG];
  double log_int[nGeoAltsG];
  double y, dy;

  float Hp_up, Hp_dn, grad_hs, grad_xp, grad_in, Hg, Xg, in, int_g, int_p;
  long index_bottom, iindex, iindexp, iiAlt;
    
  std::string function="Neutrals::calc_chapman";
  report.enter(function);

  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {

    for (iLon = 0; iLon < nGeoLonsG; iLon++) {
      for (iLat = 0; iLat < nGeoLatsG; iLat++) {

	iAlt = nGeoAltsG-1;

	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	
	// The integral from the top to infinity is the density * H (scale height)
	// So calculate the scale height:

	H = calc_scale_height(iSpecies, index, grid);
    
	integral[iAlt] = neutrals[iSpecies].density_s3gc[index] * H;
	log_int[iAlt] = log(integral[iAlt]);

	// if we wanted to do this properly, then the integral should be
	// with cell edge values and the distances from cell centers to
	// cell centers. But, we are approximating here:

	for (int iAlt=nGeoAltsG-1; iAlt>=0; iAlt--) {

	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	  indexp = ijk_geo_s3gc(iLon,iLat,iAlt+1);

	  if (iAlt < nGeoAltsG-1) {
	    integral[iAlt] = integral[iAlt+1] +
	      neutrals[iSpecies].density_s3gc[index] * grid.dalt_lower_s3gc[indexp];
	    log_int[iAlt] = log(integral[iAlt]);
	  }

	  H = calc_scale_height(iSpecies, index, grid);
      
	  xp[iAlt] = grid.radius_s3gc[index] / H;

	  // Eqn (10) Smith & Smith
	  y = sqrt(0.5 * xp[iAlt]) * fabs(grid.cos_sza_s3gc[index]);

	  // Eqn (12) Smith and Smith
	  if (y < 8) erfcy[iAlt] = (a + b*y) / (c + d*y + y*y);
	  else erfcy[iAlt] = f / (g + y);

	}
    
	// Don't need chapman integrals in the lower ghostcells:

	for (int iAlt=0; iAlt < nGeoGhosts; iAlt++) { 
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	  neutrals[iSpecies].chapman_s3gc[index] = max_chapman;
	}
	
	// Rest of domain:

	for (int iAlt=nGeoGhosts; iAlt < nGeoAltsG; iAlt++) {

	  index = ijk_geo_s3gc(iLon,iLat,iAlt);

	  // This is on the dayside:
	  if (grid.sza_s3gc[index] < pi/2 || grid.sza_s3gc[index] > 3*pi/2) {

	    neutrals[iSpecies].chapman_s3gc[index] =
	      integral[iAlt] * sqrt(0.5 * pi * xp[iAlt]) * erfcy[iAlt];

	  } else {

	    // This is on the nghtside of the terminator:

	    y = grid.radius_s3gc[index] * abs(cos(grid.sza_s3gc[index]-pi/2));

	    // This sort of assumes that nGeoGhosts >= 2:
	    index_bottom = ijk_geo_s3gc(iLon,iLat,nGeoGhosts);
	    if (y > grid.radius_s3gc[index_bottom]) {

	      iiAlt = iAlt;
	      iindex = ijk_geo_s3gc(iLon,iLat,iiAlt-1);
	      while (grid.radius_s3gc[iindex]>y) {
		iiAlt--;
		iindex = ijk_geo_s3gc(iLon,iLat,iiAlt-1);
	      }
	      iiAlt--;

	      iindexp = ijk_geo_s3gc(iLon,iLat,iiAlt+1);
	      iindex = ijk_geo_s3gc(iLon,iLat,iiAlt);

	      Hp_up = calc_scale_height(iSpecies, iindexp, grid);
	      Hp_dn = calc_scale_height(iSpecies, iindex, grid);

	      // make sure to use the proper cell spacing (iiAlt+1 & lower):
	      grad_hs = (Hp_up - Hp_dn) / grid.dalt_lower_s3gc[iindexp];
	      grad_xp = (xp[iiAlt+1]-xp[iiAlt]) / grid.dalt_lower_s3gc[iindexp];
	      grad_in = (log_int[iiAlt+1] - log_int[iiAlt]) / grid.dalt_lower_s3gc[iindexp];
	  
	      // Linearly interpolate H and X:
	      dy = y - grid.radius_s3gc[iindex];
	      Hg = Hp_dn + grad_hs * dy;
	      Xg = xp[iiAlt] + grad_xp * dy;
	      in = log_int[iiAlt] + grad_in * dy;

	      int_g = exp(in);
	      int_p = integral[iAlt];
	      // Equation (19) Smith & Smith
	      neutrals[iSpecies].chapman_s3gc[index] =
		sqrt(0.5 * pi * Xg) * (2.0 * int_g - int_p * erfcy[iAlt]);

	      if (neutrals[iSpecies].chapman_s3gc[index] > max_chapman)
		neutrals[iSpecies].chapman_s3gc[index] = max_chapman;
	  
	    } else {

	      // This says that we are in the shadow of the planet:

	      neutrals[iSpecies].chapman_s3gc[index] = max_chapman;

	    }

	  }

	  if (report.test_verbose(10))
	    std::cout << "iSpecies, iAlt, chap : " << iSpecies << " " << iAlt << " " <<
	      grid.sza_s3gc[index]*rtod << " " << 
	      xp[iAlt] << " " << 
	      erfcy[iAlt] << " " << 
	      neutrals[iSpecies].chapman_s3gc[index] << " " << integral[iAlt] << "\n";

	} // iAlt

      } // iLat
    } // iLon
    
  } // iSpecies
    
  report.exit(function);
  return;
  
}

// -----------------------------------------------------------------------------
// Calculate thermal conduction
// -----------------------------------------------------------------------------

void Neutrals::calc_conduction(Grid grid, Times time, Report &report) {

  float prandtl[nGeoAltsG];
  float rhocv[nGeoAltsG];
  float lambda[nGeoAltsG];
  float conduction[nGeoAltsG];
  float temp[nGeoAltsG];
  float dalt_lower[nGeoAltsG];
  float dt;
  
  long iLon, iLat, iAlt, index;

  std::string function="Neutrals::calc_conduction";
  report.enter(function);
  
  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	conduction_s3gc[index] = 0.0;
      }
    }
  }
	
  for (iLon = iGeoLonStart_; iLon < iGeoLonEnd_; iLon++) {
    for (iLat = iGeoLatStart_; iLat < iGeoLatEnd_; iLat++) {
      
      // Treat each altitude slice individually:
      
      for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	rhocv[iAlt] = rho_s3gc[index] * Cv_s3gc[index];
	// rhocv needs to be scaled by radius squared:
	rhocv[iAlt] = rhocv[iAlt] * grid.radius_sq_s3gc[index];
    
	// Need to make this eddy * rho * cv:
	prandtl[iAlt] = 0.0;

	lambda[iAlt] = kappa_s3gc[index] + prandtl[iAlt];
	// lambda needs to be scaled by radius squared:
	lambda[iAlt] = lambda[iAlt] * grid.radius_sq_s3gc[index];
    
	conduction[iAlt] = 0.0;

	temp[iAlt] = temperature_s3gc[index];

	dalt_lower[iAlt] = grid.dalt_lower_s3gc[index];
	
      }

      dt = time.get_dt();

      // if (iLon == nGeoLonsG/2 && iLat == nGeoLatsG/2) {
      //   for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
      // 	  std::cout << "before conduction : "
      // 		    << iLon << " " << iLat << " " << iAlt << " "
      // 		    << temp[iAlt] << " "
      // 		    << lambda[iAlt] << " " 
      // 		    << rhocv[iAlt] << " " 
      // 		    << dt << " " 
      // 		    << dalt_lower[iAlt] << "\n";
      // 	}
      // }
	    
      solver_conduction(temp, lambda, rhocv, dt, dalt_lower, conduction);

      // if (iLon == nGeoLonsG/2 && iLat == nGeoLatsG/2) {
      //   for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
      // 	  std::cout << "after conduction : "
      // 		    << iLon << " " << iLat << " " << iAlt << " "
      // 		    << conduction[iAlt] << "\n";
      // 	}
      // }
      // We want the sources to be in terms of dT/dt, while the
      // conduction actually solves for Tnew-Told, so divide by dt

      for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
	conduction_s3gc[index] = conduction[iAlt]/dt;

	if (report.test_verbose(10))
	  std::cout << "conduction : " << index
		    << " " << conduction_s3gc[index]*seconds_per_day << " deg/day\n";
	
      }

      // End of calculation for each altitude

    } // lat

  } // lon

  report.exit(function);

}

  

// -----------------------------------------------------------------------------
// Calculate EUV driven ionization and heating rates
// -----------------------------------------------------------------------------

void Neutrals::calc_ionization_heating(Euv euv, Ions &ions, Report &report) {

  long iAlt, iLon, iLat, iWave, iSpecies, index, indexp;
  int i_, idion_, ideuv_, nIonizations, iIon, iIonization;
  float tau, intensity, photoion;

  float ionization;

  std::string function="calc_ionization_heating";
  report.enter(function);

  // Zero out all source terms:
  
  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {

	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	heating_euv_s3gc[index] = 0.0;
	for (iSpecies=0; iSpecies < nSpecies; iSpecies++)
	  neutrals[iSpecies].ionization_s3gc[index] = 0.0;

      }
    }
  }
	
  for (iLon = iGeoLonStart_; iLon < iGeoLonEnd_; iLon++) {
    for (iLat = iGeoLatStart_; iLat < iGeoLatEnd_; iLat++) {
      for (iAlt = iGeoAltStart_; iAlt < iGeoAltEnd_; iAlt++) {

	index = ijk_geo_s3gc(iLon,iLat,iAlt);
      
	heating_euv_s3gc[index] = 0.0;
	for (iSpecies=0; iSpecies < nSpecies; iSpecies++)
	  neutrals[iSpecies].ionization_s3gc[index] = 0.0;
    
	for (iWave=0; iWave < euv.nWavelengths; iWave++) {

	  // Need to calculate the intensity at particular wavelength
	  // (iWave), for each species (iSpecies), which is dependent
	  // on the chapman integrals at that particular location (index)
	  // and the cross section (i_):
	
	  tau = 0.0;
	  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	    if (neutrals[iSpecies].iEuvAbsId_ > -1) {
	      i_ = neutrals[iSpecies].iEuvAbsId_;
	      tau = tau +
		euv.waveinfo[i_].values[iWave] *
		neutrals[iSpecies].chapman_s3gc[index];
	    }	  
	  }

	  intensity = euv.wavelengths_intensity_top[iWave] * exp(-1.0*tau);
	  
	  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {

	    // Calculate Photo-Absorbtion for each species and add them up:
	    i_ = neutrals[iSpecies].iEuvAbsId_; // index of photo abs cross section
	    if (i_ > -1) {
	      heating_euv_s3gc[index] = heating_euv_s3gc[index] +
		intensity *
		euv.wavelengths_energy[iWave] * 
		euv.waveinfo[i_].values[iWave] *  // cross section
		neutrals[iSpecies].density_s3gc[index];
	    }

	    for (iIonization = 0; iIonization < neutrals[iSpecies].nEuvIonSpecies; iIonization++) {

	      i_ = neutrals[iSpecies].iEuvIonId_[iIonization];
	      // std::cout << iSpecies << " " << iIonization << " " << i_ << "\n";
	      
	      ionization = 
		intensity *
		euv.waveinfo[i_].values[iWave] *  // cross section
		neutrals[iSpecies].density_s3gc[index];

	      neutrals[iSpecies].ionization_s3gc[index] =
		neutrals[iSpecies].ionization_s3gc[index] + ionization;

	      iIon = neutrals[iSpecies].iEuvIonSpecies_[iIonization];
	      ions.species[iIon].ionization_s3gc[index] = 
		ions.species[iIon].ionization_s3gc[index] + ionization;	

	    }
	     
	  } // Each species

	  
	} // Each wavelength

	// Scale heating with efficiency, and 
	// convert energy deposition to change in temperature:

	heating_euv_s3gc[index] =
	  heating_efficiency *
	  heating_euv_s3gc[index] / rho_s3gc[index] / Cv_s3gc[index];

	if (report.test_verbose(10))
	  std::cout << "heating : " << index
	       << " " << heating_euv_s3gc[index]*seconds_per_day << " deg/day\n";
	  
      } // Alts
    } // Lats
  } // Lons
	    
//	// We need to do things here:
//	// - Identify where the ionization cross section is stored
//	// - Identify which ion gets the ionization
//	// - Keep track of neutral losses
//
//	// How many ionizations do we have:
//	nIonizations = neutrals[iSpecies].iEuvIonId_.size();
//	for (int iIon=0; iIon < nIonizations; iIon++) {
//
//	  // Find cross section:
//	  ideuv_ = neutrals[iSpecies].iEuvIonId_[iIon];
//	  // Find ion to put sources:
//	  idion_ = neutrals[iSpecies].iEuvIonSpecies_[iIon];
//
//	  photoion = euv.waveinfo[ideuv_].values[iWave];
//
//	  // This is the ionization rate:
//	  ionization = intensity * photoion * 
//	    neutrals[iSpecies].density[iAlt];
//	    
//	  // Keep track of losses for neutrals:
//	  neutrals[iSpecies].ionization[iAlt] =
//	    neutrals[iSpecies].ionization[iAlt] +
//	    ionization;
//
//	  // Keep track of sources for ions:
//	  ions.species[idion_].ionization[iAlt] = 
//	    ions.species[idion_].ionization[iAlt] +
//	    ionization;
//
//	}
//	  

  report.exit(function);
  return;

}
