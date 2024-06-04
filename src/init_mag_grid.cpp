// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_mag_grid(Planets planet) {
  
  std::string function = "Grid::init_mag_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
// turn the switch on! 
  set_IsGeoGrid(false);

  // This is just an example:  
  Inputs::grid_input_struct grid_input = input.get_grid_inputs("MagGrid");
  
  int64_t iLon, iLat, iAlt;
  
  // SHOW(grid_input.dalt); exit(10);

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  arma_vec lon1d(nLons);
  precision_t dlon = (grid_input.lon_max - grid_input.lon_min) / (nLons-2*nGCs);

  for (iLon=0; iLon < nLons; iLon++)
    lon1d[iLon] = grid_input.lon_min + (iLon-nGCs+0.5) * dlon;
  
  for (iLat=0; iLat < nLats; iLat++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {

      magLon_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;
      
      // AD: does below make sense?
      // magPhi_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;

    }
  }
  

  // cout << "from"<< function<< "nLats=  "<< nLats<<"\n"<<endl;

  magPhi_scgc = magLon_scgc;

  
  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  arma_vec lat1d(nLats);
  
  //cout << "!!!!!!!!!!! HA1a "<<endl;


  arma_vec lshell(nLats);

  precision_t dlat = (grid_input.lat_max - grid_input.lat_min) / (nLats-2*nGCs);
  
  lshell(0) = 1/pow(cos(grid_input.lat_min),2.0);
  
  lshell(nLats-1) = 1/pow(cos(grid_input.lat_max),2.0);

  precision_t dlshell = (lshell(nLats-1)-lshell(0))/nLats;
  
  for (iLat=1; iLat < nLats; iLat++){

    lshell[iLat] = lshell[iLat-1]+dlshell;
        
    //<
    lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;
    //>
  }
  
  // SHOW(lshell);
  // SHOW(lat1d); exit(10);


  //<
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      this->magLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }
  //>


  for (iLon=0; iLon < nLons; iLon++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
    
      //magP_scgc.subcube(iLon, 0, iAlt, iLon, nLats-1, iAlt) = lshell;
    
      for (iLat=0; iLat<nLats; iLat++){
	      //cout << "L fill " << lshell[iLat] << endl;
	      
        magP_scgc(iLon, iLat, iAlt) = lshell[iLat];

	      //if (iLon==12 and iLat==10 and iAlt==1){
	      //cout << iLon<<" "<< iLat <<" "<< iAlt <<" P fill " << magP_scgc(iLon,iLat,iAlt) << endl;
	      //}
	      //cout << iLon<<" "<< iLat <<" "<< iAlt <<" P fill " << magP_scgc(12,10,1) << endl;
	      //cout << "magP_scgc[iLon, iLat, iAlt] = lshell[iLat];"<<" "<<magP_scgc[iLon, iLat, iAlt] <<" " <<lshell[iLat] << endl;
      }
    }
  }

  //cout << magP_scgc(12, 10, 1)<<endl;
  
  
  
  // fill along the field 
  precision_t qS;
  precision_t qN;
  precision_t Lshell;
  precision_t Lon;
  precision_t Gamma;
  double q[nZ];  

  // set the min alt and gamma factor for filling line
  precision_t AltMin = grid_input.alt_min; 
  Gamma = 2.0;
  
  
  for (iLon=0; iLon < nLons; iLon++) {
    for (iLat=0; iLat < nLats; iLat++) {

      Lshell  = magP_scgc(iLon,iLat,1);
      Lon     = magPhi_scgc(iLon,iLat,1);
      
      // get the q value for N and S hemisphere
      // cout << iLon << " " << iLat << " L= "<<Lshell<<" "<< magP_scgc(iLon,iLat,1)<<endl;
      
      auto Qvals = lshell_to_qn_qs(planet, Lshell, Lon, AltMin);

      qN = Qvals.first;
      qS = Qvals.second;

      //cout << qN << endl;
      //cout << qS << endl;

      //fill in the q array for this P and Phi
      fill_dipole_q_line(qN, qS, Gamma, nZ, Lshell, Lon, q);

      //copy this q array into
      for (iAlt=0;iAlt<nAlts; iAlt++){
	      magQ_scgc(iLon, iLat, iAlt) = q[iAlt];
      }
    }
  }
  

  // fill x y z values
  precision_t Llr[3], Xyz[3];
  int iX, iY, iZ;
  
  for (iX=0; iX<nX; iX++){
    for (iY=0; iY<nY; iY++){
      for (iZ=0; iZ<nZ; iZ++){
	    
	    
        // For given q and p we can now find cooresponding r and theta in
        // dipole coord. Starty by numerically solving for r (normalized) 
        // using equation 4 of Huba et al 2000 given by q^2r^4+1/p r -1 = 0
        auto rtheta =
          p_q_to_r_theta(magP_scgc(iX,iY,iZ), magQ_scgc(iX,iY,iZ));
        
        precision_t r = rtheta.first;
        precision_t theta = rtheta.second;
     
        //cout << "i, q " << i << "  " << q[i] << endl;
        //cout << "i, x " << i << "  " << x[i] << endl;
        //cout << "i, r, theta " << i << "  " << r[i]<<" "<<theta[i] << endl << endl;
        
        //cout << "iX+ magPhi: " << iX <<" "<< iY << " " << iZ << " " << magPhi_scgc(iX,iY,iZ) << endl;
        
        // Llr: lat, lon, rad

        Llr[0] = magPhi_scgc(iX,iY,iZ);
        Llr[1] = 0.5*cPI - theta;
        Llr[2] = r;

        //if (iZ==nZ/2){
        //  cout << magP_scgc(iX,iY,iZ) <<endl;
        //  cout << Llr[0]<<" "<<Llr[1]<<" "<<Llr[2] << endl;
        //}
        
        transform_llr_to_xyz(Llr, Xyz);

        magX_scgc(iX,iY,iZ)=Xyz[0];
        magY_scgc(iX,iY,iZ)=Xyz[1];
        magZ_scgc(iX,iY,iZ)=Xyz[2];


        precision_t radius0 = planet.get_radius(0.0);
        // SHOW(radius0)
        // exit(10);

        //<


        // if (iX == 5 & iY == 5)
        //   cout << "lon, lat, alt: " << magPhi_scgc(iX,iY,iZ) << " "
        //   << theta << " " << r 
        //  << " " << magP_scgc(iX,iY,iZ) 
        //  << " " << magQ_scgc(iX,iY,iZ) << "\n";

        this->geoLon_scgc(iX,iY,iZ) = magPhi_scgc(iX,iY,iZ);
        this->geoLat_scgc(iX,iY,iZ) = theta;
        this->geoAlt_scgc(iX,iY,iZ) = r;

        this->geoAlt_scgc(iX,iY,iZ) *= radius0;


        // - grid_input.alt_min ?

        //>
        
	    }
    }
  }





  // save 3D mag grid for examination
  std::fstream gridfile;
  gridfile.open ("grid3D.dat",ios::out);
  gridfile.precision(std::numeric_limits<long double>::digits10);
  
  // write header
  gridfile << "VARIABLES = \"X\", \"Y\", \"Z\" " << endl;
  gridfile << "Zone I = "<< nZ << ",J = " << nY << ",K = "<< nX
	   <<", DATAPACKING=POINT" << endl;
  
  // write grid data
  for (iX=0; iX<nX; iX++){
    for (iY=0; iY<nY; iY++){
      for (iZ=0; iZ<nZ; iZ++){
	
        gridfile << std::fixed << magX_scgc(iX,iY,iZ)
          <<" "<< std::fixed << magY_scgc(iX,iY,iZ)
          <<" "<< std::fixed << magZ_scgc(iX,iY,iZ)
          << endl;
      }
    }
  }
  gridfile.close();


  // save 3D mag grid slice for examination
  std::fstream gridfileslice;
  gridfileslice.open ("grid_slice.dat",ios::out);

  // write header
  gridfileslice << "VARIABLES = \"X\", \"Y\", \"Z\" " << endl;
  gridfileslice << "Zone I = " << nZ << ",J = "<< nY
	   <<", DATAPACKING=POINT" << endl;
  
  // write grid data
  for (iY=0; iY<nY; iY++){
    for (iZ=0; iZ<nZ; iZ++){
      gridfileslice << magX_scgc(1,iY,iZ)
		    <<" "<< magY_scgc(1,iY,iZ)
		    <<" "<< magZ_scgc(1,iY,iZ)
		    << endl;
    }
  }
  
  gridfileslice.close();

  IsGeoGrid = false;
  
  // Calculate the radius, etc:
  
  fill_grid_radius(planet);
  
  //  fill_grid_bfield(planet, input, report);
  
  
  // We want to now set up our p, phi, s coordinates as defined by Huba et al 2000
  // this cooresponds to x, y, and z in our grid object
  

//  float AltMin = grid_input.alt_min; 
//  float qS;
//  float qN;
//  float Lshell;
//  float Lon;
//  float Gamma;
//  
//  Lshell  = 4.0;
//  Lon     = 0.0;
//  Gamma = 2.0;
//  
//  auto Qvals = lshell_to_qn_qs(planet, Lshell, Lon, AltMin, report);
//
//  qN = Qvals.first;
//  qS = Qvals.second;
//
//  cout << qN << endl;
//  cout << qS << endl;
//
//  double q[nZ];
//  fill_dipole_q_line(qN, qS, Gamma, nZ, Lshell, Lon, q, report);
//  
  report.exit(function);
}

// ----------------------------------------------------------------------
// Routine to find q_N and q_S for a given L 
// 
// ----------------------------------------------------------------------
std::pair<precision_t,precision_t> Grid::lshell_to_qn_qs(Planets planet, precision_t Lshell, precision_t Lon, precision_t AltMin) {
  std::string function = "Grid::lshell_to_qn_qs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t qN,qS;
  
  precision_t XyzDipoleLeft[3], XyzDipoleMid[3], XyzDipoleRight[3];
  precision_t XyzGeoLeft[3], XyzGeoMid[3], XyzGeoRight[3];
  precision_t rGeoLeft, rGeoMid, rGeoRight;
  precision_t LlrDipoleLeft[3], LlrDipoleMid[3], LlrDipoleRight[3];
  precision_t ThetaTilt, PhiTilt;
  precision_t Lat, Radius, rMin;
  // Named dimension constants
  static int Lon_= 0, Lat_= 1, Radius_= 2;
 
  //bound vars for bisection search
  precision_t ThetaRight, ThetaLeft, ThetaMid;
  precision_t rDipoleLeft,rDipoleMid,rDipoleRight;
  
  //Stopping condition for bisection search
  precision_t DeltaTheta;
  precision_t Tolerance = 1e-4;
    
  // status vars for bisection search
  int iStatusLeft, iStatusRight, iStatusMid;
  // note we normalize Lshell by equatorial radius
  precision_t RadiusEq = planet.get_radius(0.0);


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
      LlrDipoleLeft[Lon_] = Lon;
      LlrDipoleLeft[Lat_] = 0.5*cPI-ThetaLeft;
      LlrDipoleLeft[Radius_] = rDipoleLeft;
      transform_llr_to_xyz(LlrDipoleLeft, XyzDipoleLeft);
      
      LlrDipoleMid[Lon_] = Lon;
      LlrDipoleMid[Lat_] = 0.5*cPI-ThetaMid;
      LlrDipoleMid[Radius_] = rDipoleMid;
      transform_llr_to_xyz(LlrDipoleMid, XyzDipoleMid);
      
      LlrDipoleRight[Lon_] = Lon;
      LlrDipoleRight[Lat_] = 0.5*cPI-ThetaRight;
      LlrDipoleRight[Radius_] = rDipoleRight;
      transform_llr_to_xyz(LlrDipoleRight, XyzDipoleRight);
      
      // Transform to XyzGeo and unnormalize
      convert_dipole_geo_xyz(planet, XyzDipoleLeft, XyzGeoLeft);
      convert_dipole_geo_xyz(planet, XyzDipoleMid, XyzGeoMid);
      convert_dipole_geo_xyz(planet, XyzDipoleRight, XyzGeoRight);
      
      //cout << "XyzGeoLeft[0]" << XyzGeoLeft[0] << endl;
      //cout << "XyzGeoLeft[1]" << XyzGeoLeft[1] << endl;
      //cout << "XyzGeoLeft[2]" << XyzGeoLeft[2] << endl;

      XyzGeoLeft[0]=XyzGeoLeft[0]*RadiusEq;
      XyzGeoLeft[1]=XyzGeoLeft[1]*RadiusEq;
      XyzGeoLeft[2]=XyzGeoLeft[2]*RadiusEq;

      abort;
      
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
    rDipoleMid    = Lshell * pow(sin(ThetaMid),2.0);
    if (iQ ==0){
      qN = pow(rDipoleMid,-2.0)*cos(ThetaMid);
      // cout << "!!! For L = " << Lshell << endl;
      // cout << "!!! qN = " << qN << endl;
      // cout << "!!! ThetaMid = " << ThetaMid*cRtoD << endl;
    }else{
      qS = pow(rDipoleMid,-2.0)*cos(ThetaMid);
      // cout << "!!! qS = " << qS << endl;
    }
  }

  report.exit(function);
  return {qN,qS};
  
}

// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//  
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, precision_t XyzDipole[3], precision_t XyzGeo[3]) {
  precision_t XyzRemoveShift[3];
  precision_t XyzRemoveTilt[3];
  precision_t XyzRemoveRot[3];

  // get planetary parameters, use radius at equator for Lshell reference
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t radius = planet.get_radius(0.0);

  
  // get the dipole shift, but normalize it to equatorial radius 
  precision_t dipole_center[3];
  std::vector<precision_t> temp_dipole_center = planet.get_dipole_center();
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

// ----------------------------------------------------------------------
// Routine to fill in the q values for a particular L and lon
// using equations 7-8 from Huba et al 2000
// ----------------------------------------------------------------------
void Grid::fill_dipole_q_line(precision_t qN, precision_t qS, precision_t Gamma, int nZ, precision_t Lshell, precision_t Lon, double *q) {
  std::string function = "Grid::fill_dipole_q_line";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
  double x[nZ];
  double r[nZ];
  double theta[nZ];
  double Dx;
  precision_t Llr[3], Xyz[3];
  int DoTestLine = 0;
  
  //open test file for writing the grid data for plotting
  std::fstream gridfile;
  if (DoTestLine==1){
    gridfile.open ("grid.dat",ios::out);
  }
  
  // set the Dx (c in equation 8 from Huba et al 2000)
  // Note this equation has a typo in it. The proper
  // version  would be found by defining the bounds of
  // x from equation 7, and then dividing that into
  // equal segments.
  //Dx = 2.0*(1.0-sinh(Gamma*qN))/((static_cast<float>(nZ)-1.0)*sinh(Gamma*qS));

  Dx = (sinh(Gamma*qS)-sinh(Gamma*qN))/((static_cast<precision_t>(nZ)-1.0)*sinh(Gamma*qS));
  
  //Dx = 2.0/(static_cast<float>(nZ)-1.0);
  //Dx = (static_cast<double>(qN)-static_cast<double>(qS))/(static_cast<double>(nZ)-1.0);

  //cout << "Dx = " << Dx << endl;
  //cout << "Gamma = " << Gamma << endl;
  //cout << "nZ = " << nZ << endl;
  //cout << "qN = " << qN << endl;
  //cout << "qS = " << qS << endl;

  // set initial x[0] value using eq. 7  from Huba et al with qi=qN
  x[0] =  sinh(Gamma*qN)/sinh(Gamma*qS);
  q[0] =  asinh(x[0]*sinh(Gamma*qS))/Gamma;
    
  // fill x_i=x_(i-1)+Dx
  for(int i = 1; i < nZ; i++){
    x[i] = x[i-1]+Dx;
    q[i] = asinh(x[i]*sinh(Gamma*qS))/Gamma;

    // For given q and p we can now find cooresponding r and theta in
    // dipole coord. Starty by numerically solving for r (normalized) using
    // equation 4 of Huba et al 2000 given by q^2r^4+1/p r -1 = 0
     auto rtheta = p_q_to_r_theta(Lshell, q[i]);
     
     r[i] = rtheta.first;
     theta[i] = rtheta.second;
     
     //cout << "i, q " << i << "  " << q[i] << endl;
     //cout << "i, x " << i << "  " << x[i] << endl;
     //cout << "i, r, theta " << i << "  " << r[i]<<" "<<theta[i] << endl << endl;

     Llr[0] = Lon;
     Llr[1] = 0.5*cPI-theta[i];
     Llr[2] = r[i];
       
     transform_llr_to_xyz(Llr, Xyz);

     if (DoTestLine==1){
       gridfile << Xyz[0] <<" " 
                << Xyz[1] <<" "
                << Xyz[2] <<" "
                << Lshell <<  endl;

     }
  }

  if (DoTestLine==1){
    gridfile.close();
  }
  
  report.exit(function);
  return ;
}



// ----------------------------------------------------------------------
// Routine to convert p and q to r and theta. Appraoch is to first solve
// for r using eq 4 from Huba et al 2000. q^2*r^4+1/q*r-1=0
// This is solved numerically using Newton-Raphson (NR) technique.
// Once we know r we can reover theta from p=r*1/(sin theta)^2.
// note r here is normalized to planet radius. 
// 
// ----------------------------------------------------------------------
std::pair<precision_t,precision_t> Grid::p_q_to_r_theta(precision_t p, precision_t q) {
  //return quanties
  precision_t r, theta;
  // function value and derivative for NR method
  precision_t Func, dFunc;
  // tolerance for root finding
  precision_t Tolerance = 0.00001;

  // initial guess for r
  r = 100.0;

  Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
  dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;
  
  // cout<< "p,q="<<p<<" "<<q << endl;
  // cout<< Func<<" "<<dFunc;
  // cout<<endl;

  int itr=0;
  int maxItr=100;

  // apply NR iterations to get 
  while( abs(Func/dFunc) > Tolerance) { 
    try {
      Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
      dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;
      
      // in NR method r(i+1)=r(i)-f/f' for each iteration
      
      r = r - Func/dFunc;

      if (++itr > maxItr){ throw(itr);}
    }
    catch (int itr){
        cout<<"WARN: exceeded max #iterations.. exiting ";
        exit(10);
    }
    // cout << r << " " << Func << " "<< dFunc << endl;
  }
  
  // now that r is determined we can solve for theta
  //theta = asin(sqrt(r/p));
  theta = acos(q*pow(r,2.0));
  
  
  //cout << "for p,q = " << p <<" "<< q << endl;
  //cout << "  r     = " << r << endl;
  //cout << "  theta = " << theta << endl;
  //cout << endl;
  
  return {r,theta};
}

arma_vec Grid::get_r3_spacing(precision_t lat, precision_t rMin, 
                        precision_t rMax, int64_t nPts, int64_t nGcs) {

  precision_t rMaxReal = rMax;

  precision_t lShell = get_lshell(lat, rMin);
  if (lShell < rMaxReal) {
    rMaxReal = lShell;
    std::cout << "Limiting rMaxReal from " << rMax << " to " << rMaxReal << "\n";
  }
  precision_t rMin3 = pow(rMin, 1.0/3.0);
  precision_t rMax3 = pow(rMaxReal, 1.0/3.0);
  precision_t dr3 = (rMax3 - rMin3) / (nPts-nGcs*2);
  arma_vec r(nPts);
  for (int64_t iPt = 0; iPt < nPts; iPt++) {
    r(iPt) = pow(rMin3 + dr3 * (iPt - nGcs), 3);
  }
  return r;
}

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_dipole_grid(Quadtree quadtree, Planets planet) {
  
  std::string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
// turn the switch on! 
  IsGeoGrid = false;

  int64_t iLon, iLat, iAlt;

  // This is just an example:  
  report.print(3, "Getting mgrid_inputs inputs in dipole grid");

  Inputs::grid_input_struct grid_input = input.get_grid_inputs("MagGrid");

  report.print(3, "Setting inputs in dipole grid");
  precision_t min_apex = grid_input.min_apex;
  precision_t min_alt = grid_input.alt_min;
  precision_t max_alt = grid_input.alt_max;
  precision_t planetRadius = planet.get_radius(0.0);
  precision_t min_lshell = (min_apex + planetRadius)/planetRadius;
  precision_t min_r = (min_alt + planetRadius)/planetRadius;
  precision_t max_r = (max_alt + planetRadius)/planetRadius;
  precision_t min_lat = get_lat_from_r_and_lshell(min_r, min_lshell);
  precision_t stretch = (cPI/2 - min_lat) / (cPI/2);
  report.print(3, "Done setting inputs in dipole grid");

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  geoLon_scgc = magLon_scgc;

  precision_t dlat = size_up_norm(1) * cPI / (nLats - 2 * nGCs);
  precision_t lat0 = lower_left_norm(1) * cPI;

  arma_vec lat1d(nLats);

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++) {
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;
    std::cout << "Original : " <<  lat1d(iLat) << " ";
    if (lat1d(iLat) >= 0) {
      lat1d(iLat) = min_lat + lat1d(iLat) * stretch;
    } else {
      lat1d(iLat) = -min_lat + lat1d(iLat) * stretch;
    }
    std::cout << "Final : " <<  lat1d(iLat) << "\n";

  }
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }

  arma_vec rNorm1d, lat1dAlong;
  arma_cube r3d(nLons, nLats, nAlts);

  precision_t lShell;

  for (iLat = 0; iLat < nLats; iLat++) {
    lat0 = lat1d(iLat);
    if (lat0 > cPI/2) lat0 = cPI - lat0;
    if (lat0 < -cPI/2) lat0 = -cPI - lat0;
    lShell = get_lshell(lat0, min_r);
    std::cout << "iLat : " << iLat << " " << nLats << "\n";
    std::cout << "lShell : " << lat0 * cRtoD << " " << lShell << "\n";
    std::cout << "min_r : " << min_r << " " << max_r << " " << nAlts << " " << nGCs << "\n";
    rNorm1d = get_r3_spacing(lat0, min_r, max_r, nAlts, nGCs);
    lat1dAlong = get_lat_from_r_and_lshell(rNorm1d, lShell);
    if (lat0 < 0)
      lat1dAlong = -1.0 * lat1dAlong;
    for (iLon = 0; iLon < nLons; iLon++) {
      r3d.tube(iLon, iLat) = rNorm1d * planetRadius;
      magLat_scgc.tube(iLon, iLat) = lat1dAlong;
    }
  }

  geoLat_scgc = magLat_scgc;
  magAlt_scgc = r3d - planetRadius;

  std::vector<arma_cube> llr, xyz, xyzRot1, xyzRot2;
  llr.push_back(magLon_scgc);
  llr.push_back(magLat_scgc);
  llr.push_back(r3d);
  xyz = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planetRadius;

  // Calculate the radius, etc:
  fill_grid_radius(planet);
  calc_rad_unit(planet);
  calc_gravity(planet);

  std::cout << "magLon : " << magLon_scgc(12,10,5) * cRtoD << " "
              << magLat_scgc(12,10,5) * cRtoD << " "
              << magAlt_scgc(12,10,5) / 1000.0 << "\n";

  std::cout << "geoLon : " << geoLon_scgc(12,10,5) * cRtoD << " "
              << geoLat_scgc(12,10,5) * cRtoD << " "
              << geoAlt_scgc(12,10,5) / 1000.0 << "\n";

  report.exit(function);
  return;

}

