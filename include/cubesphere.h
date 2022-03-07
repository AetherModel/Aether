#ifndef INCLUDE_CUBESPHERE_H_
#define INCLUDE_CUBESPHERE_H_

#include "aether.h"
#include <memory>
#include <armadillo>

/*************************************************
 * \brief A namespace withe all quad sphere grid logic.
 *************************************************/
namespace CubeSphere {

/** Returns x, y and z locations of a cube from -1 to +1.
 **/
arma_mat map_cube2sphere(const arma_mat& cube);

/** Helper function the transforms cartesian coordinates to cube sphere grid
 * vertices in cartesian coordinates.
 */
arma_vec cube2sphere(const arma_vec& x, const arma_vec& y, const arma_vec& z);

/** Helper function the transforms cartesian coordinates to cube sphere grid
 * vertices in cartesian coordinates.
 */
precision_t cube2sphere(const precision_t& x, const precision_t& y, const precision_t& z);

/** Helper function that transforms cartesian coordinates of cube to cube sphere
 * using an x, y, z vector as input.
 */
arma_vec cube2sphere(const arma_vec& p);

/** Returns normalized x, y, z locations based on number of subdivisions per
 * face
 */
arma_cube sphere_face_locations(const uint& face, const uint& subdivisions);

/** Returns the normalized x, y, z coordinate locations of a cell in the
 * cubesphere grid.
 */
arma_vec cartesian_location(const uint& face, const uint& subdivisions,
                            const uint& row, const uint& column);

/** Quad Tree Node base class.
 */
class QTNode {
      unsigned short m_pixnum;
      unsigned int m_max_depth = 1;
      unsigned int m_depth = 0;
      std::shared_ptr<QTNode> m_children[4];
      /// Position in cube.
      arma_vec m_cube_position = {0, 0, 0};
      /// Center position in cube sphere.
      arma_vec m_sphere_position = {0, 0, 0};
      // other details to add

      public:
      /// Constructor
      QTNode(Grid grid,
             const uint& face,
             const arma_vec& cube_position,
             const uint& depth,
             const uint& max_depth);
      /// Destructor
      ~QTNode();
      /// Returns if the node is a leaf based on resolution.
      bool is_leaf();
      };

/// The normalized origins of each face of the cube (i.e. corner)
static const arma_mat ORIGINS = {
    {-1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0},
    {0.0, 0.0, -1.0},
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
};

/// Normalized right steps in cube
static const arma_mat RIGHTS = {
		{2.0, 0.0, 0.0},
		{0.0, 0.0, 2.0},
		{-2.0, 0.0, 0.0},
		{0.0, 0.0, -2.0},
		{2.0, 0.0, 0.0},
		{2.0, 0.0, 0.0}
};

/// Normalized right steps in cube
static const arma_mat UPS = {
		{0.0, 2.0, 0.0},
		{0.0, 2.0, 0.0},
		{0.0, 2.0, 0.0},
		{0.0, 2.0, 0.0},
		{0.0, 0.0, 2.0},
		{0.0, 0.0, -2.0}
};

} // CubeSphere::

class CubeSphereGrid : public Grid {
};

#endif  // INCLUDE_CUBESPHERE_H_
