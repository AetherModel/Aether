#ifndef INCLUDE_CUBESPHERE_H_
#define INCLUDE_CUBESPHERE_H_

#include "aether.h"
#include <armadillo>

/*************************************************
 * \brief A class to handle quad sphere grids.
 *************************************************/
class CubeSphere {
  public:
    /// Constructor of a quad sphere grid with side of one size being \f$ n \f$.
    CubeSphere(const int n);
    /// Returns size of one side \f$ n \f$.
    const unsigned int size_side();
    /// Returns size of cube sphere \f$ n \times n \times 6 \f$.
    const unsigned int size();
    /** Returns x, y and z locations of a cube from -1 to +1.
     *
      */
    const arma_mat map_cube2sphere(const arma_mat& cube) {
      const arma_mat sphere(arma::size(cube));
      const arma_vec x2_2 = (cube(0)%cube(0))/2;
      const arma_vec y2_2 = (cube(1)%cube(1))/2;
      const arma_vec z2_2 = (cube(2)%cube(2))/2;

      sphere(0) = c2s(cube(0), y2_2, z2_2);
      sphere(1) = c2s(cube(1), z2_2, x2_2);
      sphere(2) = c2s(cube(2), x2_2, y2_2);

      return sphere;
    }

  private:
    /// The size of one side of a cube sphere.
    unsigned int n;

    const arma_vec c2s(const arma_vec& x,
                       const arma_vec& y2_2,
                       const arma_vec& z2_2) {
      return x%(arma::sqrt(arma::ones<arma_vec>-y2_2-z2_2+y2_2%z2_2*(2/3)));
    }
};

#endif INCLUDE_CUBESPHERE_H_
