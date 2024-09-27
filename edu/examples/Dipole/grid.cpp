
#include <iostream>
#include <math.h>
#include <armadillo>
#include <fstream>

using arma_vec = arma::Col<double>;
using arma_mat = arma::Mat<double>;

inline double deg2rad(double degrees)
{
  // function compiled inline to convert degrees to radians
  static const double pi_on_180 = 4.0 * atan(1.0) / 180.0;
  return degrees * pi_on_180;
}

inline double rad2deg(double degrees)
{
  // function compiled inline to convert degrees to radians
  static const double pi_on_180 = 4.0 * atan(1.0) / 180.0;
  return degrees / pi_on_180;
  // return degrees;
}

double qp_solve(double q, double p)
{
  double term0 = 256.0 / 27.0 * pow(q, 2.0) * pow(p, 4.0);
  double term1 = pow((1.0 + sqrt(1.0 + term0)), 2.0 / 3.0);
  double term2 = pow(term0, 1.0 / 3.0);
  double term3 = 0.5 * pow(((pow(term1, 2) + term1 * term2 + pow(term2, 2)) / term1), 3.0 / 2.0);
  double new_r = p * (4.0 * term3) / (1.0 + term3) / (1.0 + sqrt(2.0 * term3 - 1.0));
  return new_r;
}

int main()
{
  int nF = 100;
  int nZ = 200;
  float alt_min = 90.0;
  float min_blat = 2.25;
  float max_blat = 89.9;

  // user-defined constants
  double gams = 0.2;
  double baselat_spacing = 6.0;

  // constants
  int iErr = 0;
  float Re = 6371.0;
  float pi = 3.14159;
  float pio2 = pi / 2.0;
  bool didWork = true;

  int nFby2 = nF / 2;
  int nZby2 = nZ / 2;
  double altmin_inRe = (alt_min + Re) / Re;

  // outputs
  arma_mat bLats(nF, nZ), bLons(nF, nZ), bAlts(nF, nZ);

  arma_vec baseLats(nF);
  //, bLons(nF, nZ), bAlts(nF, nZ);

  // TODO: REFACTOR FROM HERE TO **1 ?
  // Lay down baseLat spacing according to an exponential factor:
  double del_lat, blat_min_, blat_max_, tmp_lat;
  blat_min_ = cos(deg2rad(pow(min_blat, 1.0 / baselat_spacing)));
  blat_max_ = cos(deg2rad(pow(max_blat, 1.0 / baselat_spacing)));
  del_lat = (blat_max_ - blat_min_) / (nF - 1.0);

  for (int i = 0; i < nF; i++)
  {
    // first put down "linear" spacing
    tmp_lat = blat_min_ + del_lat * i;
    // then scale it according to the exponent & convert back to deg
    tmp_lat = pow(rad2deg(acos(tmp_lat)), baselat_spacing);
    // place values in array backwards, S => N hemis
    baseLats(nF - i - 1) = -tmp_lat;
  }

  // Find L-Shell for each baseLat
  // using L=R/sin2(theta), where theta is from north pole
  arma_vec Lshells(nF);
  for (int i = 0; i < nF; i++)
  {
    Lshells(i) = ((alt_min + Re) / Re) / pow(sin(deg2rad(90.0 - baseLats(i))), 2.0);
  }

  // allocate & calculate some things outside of the main loop
  // fa, fb, fc are factors to make the code easier to read
  double q_S, q_N, delqp, fb;
  double qp0, fb0, ft, delq, qp2;
  arma_vec exp_q_dist(nZ), q_vals(nZ);

  for (int i = 0; i < nZ; i++)
  {
    exp_q_dist(i) = gams + (1 - gams) * exp(-pow(((i - nZby2) / (nZ / 10.0)), 2.0));
  }

  for (int i_nF = 0; i_nF < nF; i_nF++)
  {

    // min/max q
    q_N = cos(deg2rad(90.0 + baseLats(i_nF))) / pow((alt_min + Re) / Re, 2.0);
    q_S = cos(deg2rad(90 - baseLats(i_nF))) / pow((alt_min + Re) / Re, 2.0);

    // calculate const. stride similar to sami2/3 (huba & joyce 2000)
    // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
    // first loop for southern hemisphere, second for north.
    for (int i_nZ = 0; i_nZ < nZby2; i_nZ++)
    {
      delqp = (q_N - q_S) / nZ;
      qp0 = q_S + i_nZ * (delqp);
      delqp = altmin_inRe * delqp;
      fb0 = (1 - exp_q_dist(i_nZ)) / exp(-q_S / delqp - 1);
      ft = exp_q_dist(i_nZ) - fb0 + fb0 * exp(-(qp0 - q_S) / delqp);

      delq = qp0 - q_S;
      qp2 = q_S + ft * delq;

      bAlts(i_nF, i_nZ) = qp_solve(qp2, Lshells(i_nF));
      bLats(i_nF, i_nZ) = rad2deg(asin(qp2 * pow(bAlts(i_nF, i_nZ), 2.0)));

      //  test mirroring across hemi's

      bAlts(i_nF, nZ - i_nZ - 1) = qp_solve(-qp2, Lshells(i_nF));
      bLats(i_nF, nZ - i_nZ - 1) = -bLats(i_nF, i_nZ);
    }
  }

  // calculate const. stride similar to sami2/3 (huba & joyce 2000)
  // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
  // first loop for southern hemisphere, second for north.
  // for (int i_nZ=nZ; i_nZ>nZby2; i_nZ++){
  // delqp = (q_N - q_S)/nZ;
  // qp0 = q_S + i_nZ*(delqp);
  // delqp = altmin_inRe * delqp;
  // fb0 = (1-exp_q_dist(i_nZ)) / exp(-q_S/delqp - 1);
  // fa = exp_q_dist(i_nZ) - fb0;
  // ft = fa + fb0 * exp(-(qp0 - q_S)/delqp);
  //
  // delq = qp0 - q_S;
  // qp2 = q_S + ft * delq;
  ////fb = (exp_q_dist(i_nZ) - fa) + fa * exp(-(exp_q_dist(i_nZ) - q_S) / delqp);
  ////fb = (exp_q_dist(i_nZ) - fa) + fa * exp(-(q_S + i_nZ * delqp/altmin_inRe) / i_nZ);
  // bAlts(i_nF, i_nZ) = qp_solve(qp2,  Lshells(i_nF));
  // bLats(i_nF, i_nZ) = rad2deg(asin(qp2 * pow(bAlts(i_nF, i_nZ), 2.0)));
  // }
  //  }
  //}
  //}
  //}

  std::ofstream fout;
  fout.open("grid.csv");
  fout << "nf,nz,lat,alt,baselat\n";
  for (int fi = 0; fi < nF; fi++)
  {
    for (int zi = 0; zi < nZ; zi++)
    {
      fout << fi << "," << zi << "," << bLats(fi, zi) << "," << bAlts(fi, zi) << "," << baseLats(fi) << "\n";
    }
  }
  fout.close();

  return iErr;
}
