// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_mag_grid(Planets planet, Inputs input, Report &report) {
  
  std::string function = "Grid::init_mag_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
  // This is just an example:
  
  Inputs::grid_input_struct grid_input = input.get_mag_grid_inputs();
  
  int64_t iLon, iLat, iAlt;
  
  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lon1d(nLons);
  float dlon = (grid_input.lon_max - grid_input.lon_min) / (nLons-2*nGCs);
  for (iLon=0; iLon < nLons; iLon++)
    lon1d(iLon) = grid_input.lon_min + (iLon-nGCs+0.5) * dlon;
  
  //  for (iLat=0; iLat < nLats; iLat++) {
  //  for (iAlt=0; iAlt < nAlts; iAlt++) {
  //    geoLon_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;
  //  }
  
  
  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lat1d(nLats);
  float dlat = (grid_input.lat_max - grid_input.lat_min) / (nLats-2*nGCs);
  for (iLat=0; iLat < nLats; iLat++)
    lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;
  
  //for (iLon=0; iLon < nLons; iLon++) {
  //  for (iAlt=0; iAlt < nAlts; iAlt++) {
  //   geoLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats-1, iAlt) = lat1d;
  //  }
  //}
  
  fvec alt1d(nAlts);
  
  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) =
        grid_input.alt_min +
        (iAlt-nGeoGhosts) * grid_input.dalt;
  }
  //for (iLon = 0; iLon < nLons; iLon++) {
  //  for (iLat = 0; iLat < nLats; iLat++) {
  //    geoAlt_scgc.tube(iLon, iLat) = alt1d;
  //}
  //}

  IsGeoGrid = 0;
  
  // Calculate the radius, etc:
  
  //  fill_grid_radius(planet, report);
  
  //  fill_grid_bfield(planet, input, report);
  
  
  // We want to now set up our p, phi, s coordinates as defined by Huba et al 2000
  // this cooresponds to x, y, and z in our grid object
  

  float AltMin = grid_input.alt_min; 
  float qS;
  float qN;
  float Lshell;
  Lshell = 4.0;
  lshell_to_qn_qs(planet, Lshell, AltMin,qN, qS,report);
  
  report.exit(function);
}

// ----------------------------------------------------------------------
// Routine to find q_N and q_S for a given L 
// 
// ----------------------------------------------------------------------
void Grid::lshell_to_qn_qs(Planets planet, float Lshell, float AltMin, float qN, float qS,Report &report) {
  std::string function = "Grid::lshell_to_qn_qs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  float XyzDipoleLeft[3], XyzDipoleMid[3], XyzDipoleRight[3];
  float XyzGeoLeft[3], XyzGeoMid[3], XyzGeoRight[3];
  float rGeoLeft, rGeoMid, rGeoRight;
  float LlrDipoleLeft[3], LlrDipoleMid[3], LlrDipoleRight[3];
  float ThetaTilt, PhiTilt;
  float Lat, Radius, rMin;
  // Named dimension constants
  static int Lon_= 0, Lat_= 1, Radius_= 2;
 
  //bound vars for bisection search
  float ThetaRight, ThetaLeft, ThetaMid;
  float rDipoleLeft,rDipoleMid,rDipoleRight;
  
  //Stopping condition for bisection search
  float DeltaTheta;
  float Tolerance = 1e-4;
    
  // status vars for bisection search
  int iStatusLeft, iStatusRight, iStatusMid;
  // note we normalize Lshell by equatorial radius
  float RadiusEq = planet.get_radius(0.0);


  

  // loop for qN and qS
  for(int iQ = 0; iQ < 2; iQ++){

    if (iQ == 0){
      // set initial left, mid, right bounds for bisection search for qN
      ThetaRight = 0.5*cPI;
      ThetaLeft = 1.0*cDtoR;
      ThetaMid = 0.5*(ThetaRight+ThetaLeft);
    }else{
      // set initial left, mid, right bounds for bisection search for qS
      ThetaLeft = 0.5*cPI;
      ThetaRight = 179.0*cDtoR;
      ThetaMid = 0.5*(ThetaRight+ThetaLeft);
    }
      
    // Initial stopping condition stopping condition
    DeltaTheta = abs(ThetaLeft-ThetaRight);
    
    
    // start bisection search for qN
    while( DeltaTheta > Tolerance ) {
      
      // find rDipole that cooresponds to these Left,Mid,Right
      // ThetaDipole values 
      rDipoleLeft   = Lshell * pow(sin(ThetaLeft),2.0);
      rDipoleMid    = Lshell * pow(sin(ThetaMid),2.0);
      rDipoleRight  = Lshell * pow(sin(ThetaRight),2.0);
      
      // Compute XyzDipole for left, mid,right states
      LlrDipoleLeft[Lon_] = 0.0;
      LlrDipoleLeft[Lat_] = 0.5*cPI-ThetaLeft;
      LlrDipoleLeft[Radius_] = rDipoleLeft;
      transform_llr_to_xyz(LlrDipoleLeft, XyzDipoleLeft);
      
      LlrDipoleMid[Lon_] = 0.0;
      LlrDipoleMid[Lat_] = 0.5*cPI-ThetaMid;
      LlrDipoleMid[Radius_] = rDipoleMid;
      transform_llr_to_xyz(LlrDipoleMid, XyzDipoleMid);
      
      LlrDipoleRight[Lon_] = 0.0;
      LlrDipoleRight[Lat_] = 0.5*cPI-ThetaRight;
      LlrDipoleRight[Radius_] = rDipoleRight;
      transform_llr_to_xyz(LlrDipoleRight, XyzDipoleRight);
      
      // Transform to XyzGeo and unnormalize
      convert_dipole_geo_xyz(planet, XyzDipoleLeft, XyzGeoLeft);
      convert_dipole_geo_xyz(planet, XyzDipoleMid, XyzGeoMid);
      convert_dipole_geo_xyz(planet, XyzDipoleRight, XyzGeoRight);
      
      XyzGeoLeft[0]=XyzGeoLeft[0]*RadiusEq;
      XyzGeoLeft[1]=XyzGeoLeft[1]*RadiusEq;
      XyzGeoLeft[2]=XyzGeoLeft[2]*RadiusEq;
      
      XyzGeoMid[0]=XyzGeoMid[0]*RadiusEq;
      XyzGeoMid[1]=XyzGeoMid[1]*RadiusEq;
      XyzGeoMid[2]=XyzGeoMid[2]*RadiusEq;
      
      XyzGeoRight[0]=XyzGeoRight[0]*RadiusEq;
      XyzGeoRight[1]=XyzGeoRight[1]*RadiusEq;
      XyzGeoRight[2]=XyzGeoRight[2]*RadiusEq;
      
      // Compute radius in geo coordinate for comparison to rmin
      rGeoLeft = 
	sqrt(pow(XyzGeoLeft[0],2)+pow(XyzGeoLeft[1],2)+pow(XyzGeoLeft[2],2));
      rGeoMid = 
	sqrt(pow(XyzGeoMid[0],2)+pow(XyzGeoMid[1],2)+pow(XyzGeoMid[2],2));
      rGeoRight =
	sqrt(pow(XyzGeoRight[0],2)+pow(XyzGeoRight[1],2)+pow(XyzGeoRight[2],2));
      
      // get rmin for given latitude. Radius is lat dependent in general.
      // also find status in (0) or out (1) of rMin
      Lat    = 0.5*cPI-acos(XyzGeoLeft[2]/rGeoLeft);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoLeft < rMin){
	iStatusLeft = 0;  
      }else{
	iStatusLeft = 1;  
      }
      
      Lat    = 0.5*cPI-acos(XyzGeoMid[2]/rGeoMid);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoMid < rMin){
	iStatusMid = 0;  
      }else{
	iStatusMid = 1;  
      }
      
      Lat    = 0.5*cPI-acos(XyzGeoRight[2]/rGeoRight);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoRight < rMin){
	iStatusRight = 0;  
      }else{
	iStatusRight = 1;  
      }
      
      // Use status values to update left, right and mid values of theta
      if (iStatusMid == 0) {
	if (iStatusRight == 1){
	  // Mid becomes left and right stays right
	  ThetaLeft = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}else{
	  // Mid becomes right and left stays left
	  ThetaRight = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}
      }else{
	if (iStatusRight == 0){
	  // Mid becomes left and right stays right
	  ThetaLeft = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}else{
	  // Mid becomes right and left stays left
	  ThetaRight = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}
      }
      // Update stopping condition
      DeltaTheta = abs(ThetaLeft-ThetaRight);
    }
    
    //set the q value
    if (iQ ==0){
      qN = pow(Lshell,-2.0)*cos(ThetaMid);
      cout << "!!! qN = " << qN << endl;
    }else{
      qS = pow(Lshell,-2.0)*cos(ThetaMid);
      cout << "!!! qS = " << qS << endl;
    }
  }
  report.exit(function);
  return ;
  
}

// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//  
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, float XyzDipole[3], float XyzGeo[3]) {
  float XyzRemoveShift[3];
  float XyzRemoveTilt[3];
  float XyzRemoveRot[3];

  // get planetary parameters, use radius at equator for Lshell reference
  float magnetic_pole_tilt = planet.get_dipole_tilt();
  float magnetic_pole_rotation = planet.get_dipole_rotation();
  float radius = planet.get_radius(0.0);

  
  // get the dipole shift, but normalize it to equatorial radius 
  float dipole_center[3];
  std::vector<float> temp_dipole_center = planet.get_dipole_center();
  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  dipole_center[0]=dipole_center[0]/radius;
  dipole_center[1]=dipole_center[1]/radius;
  dipole_center[2]=dipole_center[2]/radius;

  // Remove Tilt
  transform_rot_y(XyzDipole, magnetic_pole_tilt, XyzRemoveTilt);

  // Remove Rot
  transform_rot_z(XyzRemoveTilt, magnetic_pole_rotation, XyzRemoveRot);

  // Remove Shift
  vector_add(XyzRemoveRot, dipole_center, XyzGeo);

//  cout << "XyzDipole[0]" << XyzDipole[0] << endl;
//  cout << "XyzDipole[1]" << XyzDipole[1] << endl;
//  cout << "XyzDipole[2]" << XyzDipole[2] << endl;
//
//  cout << "XyzGeo[0]" << XyzGeo[0] << endl;
//  cout << "XyzGeo[1]" << XyzGeo[1] << endl;
//  cout << "XyzGeo[2]" << XyzGeo[2] << endl;
  
}



//// ----------------------------------------------------------------------
//// A routine for converting spherical geographic to spherical tilted
//// using equations 81-83 from Huba et al 2000
//// ----------------------------------------------------------------------
//void Grid::sph_geo_to_sph_tilted(Planets planet, float ThetaGeo, float PhiGeo,Report &report) {
//  std::string function = "Grid::sph_geo_to_sph_tilted";
//  static int iFunction = -1;
//  report.enter(function, iFunction);
//  
//  float ThetaTilt, PhiTilt;
//
//  float DipoleTilt = planet.get_dipole_tilt();
//  float magnetic_pole_rotation = planet.get_dipole_rotation();
//  
//  ThetaTilt = asin(sin(ThetaGeo)*sin(ThetaN);
//  
//  
//  
//  
//  report.exit(function);
//  return ThetaTilt, PhiTilt;
//  
//}
