// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/sizes.h"

// Assume that lambda and front are scaled by radius squared!

int solver_conduction(float value[nGeoAltsG],
		      float lambda[nGeoAltsG],
		      float front[nGeoAltsG],
		      float dt,
		      float dalt_lower[nGeoAltsG],
		      float *conduction) {

  int iErr = 0;

  int iAlt;

  float r[nGeoAltsG];
  float di[nGeoAltsG];
  float m[nGeoAltsG];
  float du[nGeoAltsG];
  float dl[nGeoAltsG];
  float du12[nGeoAltsG];
  float du22[nGeoAltsG];
  float lou[nGeoAltsG];
  float a[nGeoAltsG];
  float b[nGeoAltsG];
  float c[nGeoAltsG];
  float d[nGeoAltsG];
  float cp[nGeoAltsG];
  float dp[nGeoAltsG];
  float result[nGeoAltsG];
  
  // Only solve for conduction in the middle section.
  conduction[0] = 0.0;
  conduction[nGeoAltsG-1] = 0.0;
    
  for (iAlt=1; iAlt < nGeoAltsG-1; iAlt++) {

    di[iAlt] = lambda[iAlt];
    m[iAlt] = dt / front[iAlt];

    // This is the dalt upper (lower @ iAlt+1)
    du[iAlt] = dalt_lower[iAlt+1];
    // This is the dalt lower
    dl[iAlt] = dalt_lower[iAlt];

    r[iAlt] = du[iAlt]/dl[iAlt];

    du12[iAlt] = du[iAlt]*du[iAlt] * (1+r[iAlt])*(1+r[iAlt]);
    du22[iAlt] = 0.5 * (dl[iAlt]*du[iAlt] + du[iAlt]*du[iAlt]);
    lou[iAlt] = di[iAlt]/du22[iAlt];

  }
    
  for (iAlt=2; iAlt < nGeoAltsG-2; iAlt++) {

    dl[iAlt] =
      di[iAlt+1] -
      di[iAlt-1] * r[iAlt] * r[iAlt] -
      di[iAlt] * (1.0-r[iAlt]*r[iAlt]);

    a[iAlt] = di[iAlt]/du22[iAlt]*r[iAlt] - dl[iAlt]/du12[iAlt]*r[iAlt]*r[iAlt];
    c[iAlt] = di[iAlt]/du22[iAlt] + dl[iAlt]/du12[iAlt];
    b[iAlt] =
      - 1.0/m[iAlt]
      - di[iAlt]/du22[iAlt]*(1.0+r[iAlt])
      - dl[iAlt]/du12[iAlt]*(1.0-r[iAlt]*r[iAlt]);
    d[iAlt] = -1.0*value[iAlt]/m[iAlt];

  }

  // Lower BCs:
  a[1] = 0.0;
  b[1] = -1.0;
  c[1] = 0.0;
  d[1] = -1.0*value[1];

  // Upper BCs:
  // This assumes a constant-gradient BC:
  //     (For neutral temperature, this isn't really needed, since it is iso-thermal
  //      in the upper thermosphere, but in the ionosphere, the electron and
  //      ion temperatures are typically sloped.)
  iAlt = nGeoAltsG-2;
  a[iAlt] = 1.0*( r[iAlt]*(1.0+r[iAlt])*di[iAlt]*m[iAlt]/du22[iAlt]);
  b[iAlt] = -1.0*( 1.0 + r[iAlt]*(1+r[iAlt])*di[iAlt]*m[iAlt]/du22[iAlt]);
  c[iAlt] = 0.0;
  d[iAlt] = -1.0*value[iAlt];

  cp[1] = c[1]/b[1];
  for (iAlt=2; iAlt <= nGeoAltsG-2; iAlt++) 
    cp[iAlt] = c[iAlt]/(b[iAlt]-cp[iAlt-1]*a[iAlt]);

  dp[1] = d[1]/b[1];
  for (iAlt=2; iAlt <= nGeoAltsG-2; iAlt++) 
    dp[iAlt] = (d[iAlt]-dp[iAlt-1]*a[iAlt])/(b[iAlt]-cp[iAlt-1]*a[iAlt]);

  result[nGeoAltsG-2] = dp[nGeoAltsG-2];
  for (iAlt=nGeoAltsG-3; iAlt>0; iAlt--)
    result[iAlt] = dp[iAlt] - cp[iAlt]*result[iAlt+1];
  
  conduction[0] = 0.0;
  for (iAlt=1; iAlt<nGeoAltsG-2; iAlt++) {
    conduction[iAlt] = result[iAlt] - value[iAlt];
    //cout << "Conduction : " << iAlt << " " <<  conduction[iAlt]*seconds_per_day << " (deg/day)\n";
  }
  conduction[nGeoAltsG-1] = conduction[nGeoAltsG-2];

  return iErr;

}

    
