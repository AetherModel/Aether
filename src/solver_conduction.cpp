// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <armadillo>
#include <iostream>
using namespace arma;

// Assume that lambda and front are scaled by radius squared!


fvec solver_conduction_new(fvec value,
			   fvec lambda,
			   fvec front,
			   float dt,
			   fvec dx) {

  // dx is lower
  
  long nPts = value.n_elem;
  
  fvec di = lambda;
  fvec m = dt / front;

  // These are to allow for a stretched grid:

  fvec du(nPts);
  du(span(0,nPts-2)) = dx(span(1,nPts-1));
  du(nPts-1) = du(nPts-2);
  fvec dl = dx;
  fvec r = du/dl;
  fvec du12 = du % du % (1 + r) % (1 + r);
  fvec du22 = 0.5 * (dl % du + du % du);
  fvec lou = di/du22;

  fvec conduction(nPts);
  conduction.zeros();
  
  long i;
  
  for (i = 2; i < nPts-2; i++) 
    dl(i) = di(i+1) - di(i-1) * r(i) * r(i) - di(i) * (1.0-r(i)*r[i]);
    
  fvec a = di / du22 % r - dl / du12 % r % r;
  fvec c = di / du22 + dl / du12;
  fvec b = -1.0 / m - di / du22 % (1.0 + r) - dl / du12 % (1.0 - r % r);
  fvec d = -1.0 * value / m;

  // Lower BCs (fixed value):
  a(1) = 0.0;
  b(1) = -1.0;
  c(1) = 0.0;
  d(1) = -1.0*value(1);

  // Upper BCs:
  // This assumes a constant-gradient BC:
  //     (For neutral temperature, this isn't really needed, since it is iso-thermal
  //      in the upper thermosphere, but in the ionosphere, the electron and
  //      ion temperatures are typically sloped.)

  i = nPts-2;
  a(i) = 1.0*( r(i)*(1.0+r(i))*di(i)*m(i)/du22(i));
  b(i) = -1.0*( 1.0 + r(i)*(1+r(i))*di(i)*m(i)/du22(i));
  c(i) = 0.0;
  d(i) = -1.0*value(i);

  fvec cp(nPts,fill::zeros);
  fvec dp(nPts,fill::zeros);
    fvec result(nPts,fill::zeros);

  cp(1) = c(1)/b(1);
  for (i = 2; i <= nPts-2; i++) 
    cp(i) = c(i)/(b(i)-cp(i-1)*a(i));

  dp(1) = d(1)/b(1);
  for (i = 2; i <= nPts-2; i++) 
    dp(i) = (d(i)-dp(i-1)*a(i))/(b(i)-cp(i-1)*a(i));

  result(nPts-2) = dp(nPts-2);
  for (i = nPts-3; i > 0; i--)
    result(i) = dp(i) - cp(i)*result(i+1);
  
  conduction = result - value;
  conduction(0) = 0.0;
  conduction(1) = 0.0;
  conduction(nPts-2) = 0.0;
  conduction(nPts-1) = 0.0;

  return conduction;

}
