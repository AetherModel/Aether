// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_QUADTREE_H_
#define INCLUDE_QUADTREE_H_

#include "aether.h"

/**************************************************************
 * \class Quadtree
 *
 * \brief Defines the quadtree for blocks 
 * 
 * Aether is logically an i, j, k grid structure.  Aether does domain
 * decomposition on the (i, j) coordinates and each processor works on
 * the full domain of the k dimension.  (e.g., nn spherical
 * coordinates, i = longitude, j = latitude, and k = altitude.)  The
 * quadtree takes the (i, j) dimensions and makes blocks out of them.
 * 
 * In the quadtree there are a number of root nodes, which are then
 * divided into 2 x 2 blocks.  Each of those can then be subdivided
 * into 2 x 2 blocks.  Each block resides on a separate processor.
 * The depth of the tree is then related to the number of processors
 * as N_roots * 4 ^ Depth, where Depth = 0 are the root nodes alone.
 *
 * For the sphere, N_roots = 1, while for the cubesphere N_roots = 6.
 *
 * In the logical (i, j) block, left/right is i, which down/up is j.
 * conventions:
 * du = down / up
 * lr = left / right
 *
 * \author Aaron Ridley
 *
 * \date 2022/07/05 
 *
 **************************************************************/

class Quadtree {

public:

  /// number of blocks in each direction:
  const uint64_t nLR = 2;
  const uint64_t nDU = 2;

  struct qtnode {

    /// Vector of children for the quadtree
    std::vector<qtnode> children;

    /// The depth of the certain leaf:
    uint64_t depth;

    /// normalized coordinate of lower left corner:
    arma_vec lower_left_norm = {0.0, 0.0, 0.0};

    /// normalized sizes:
    arma_vec size_right_norm = {0.0, 0.0, 0.0};
    arma_vec size_up_norm = {0.0, 0.0, 0.0};

    /// The processor of the leaf/node:
    uint64_t iProcNode;
    /// The processor in the down direction:
    uint64_t iProc_down;
    /// The processor in the up direction:
    uint64_t iProc_up;
    /// The processor in the left direction:
    uint64_t iProc_left;
    /// The processor in the right direction:
    uint64_t iProc_right;
    /// For the CubeSphere, this defines the side of the cube:
    uint64_t iSide;

    /// The root node in the down direction:
    uint64_t iRoot_down;
    /// The root node in the up direction:
    uint64_t iRoot_up;
    /// The root node in the left direction:
    uint64_t iRoot_left;
    /// The root node in the right direction:
    uint64_t iRoot_right;
    /// The root node of the current node / leaf:
    uint64_t iRoot;
  };

  /// Number of root nodes:
  int64_t nRootNodes;
  
  /// The quadtree root nodes:
  std::vector<struct qtnode> root_nodes;

  /// Maximum depth of the tree, defined by the number of processors:
  int64_t max_depth;

  /// For the given processor, the normalized limits of the block:
  arma_vec limit_low = {0.0, 0.0, 0.0};
  arma_vec limit_high = {0.0, 0.0, 0.0};
  /// For the given processor, the side that it is on:
  uint64_t iSide = -1;
  
  /**********************************************************************
     \brief Initializes the quadtree
   **/
  Quadtree();

  /**********************************************************************
     \brief Builds the quadtree
   **/
  void build(); 

  /**********************************************************************
     \brief Makes a new node on the quadtree, recursively
     \param lower_left_norm_in The normalized coord of the lower left corner
     \param size_right_norm_in The normalized size in the right/left direc.
     \param size_up_norm_in The normalized size in the up/down direction.
     \param iProc_in_out Which processor contains the leaf / node
     \param depth_in the current depth of the leaf / node
     \param iSide basically the root node, or the side of the cubesphere
   **/
  qtnode new_node(arma_vec lower_left_norm_in,
		  arma_vec size_right_norm_in,
		  arma_vec size_up_norm_in,
		  uint64_t &iProc_in_out,
		  uint64_t depth_in,
		  uint64_t iSide);

  /**********************************************************************
     \brief Get different vectors from the node
     \param node which node to get the vector from
     \param which defines the vector to get:
                  LL = lower left;
                  SR = size in the right/left direction; 
                  SU = size in the up/down direction;
                  MID = mid point of the node;
   **/
  arma_vec get_vect(qtnode node, std::string which);
  arma_vec get_vect(std::string which);

  /**********************************************************************
     \brief finds the node of the given normalized point.  Returns
            -1 if not found on the given node.
     \param point the x, y, z normalized coordinate of the point.
     \param node search on the given node for the point.
   **/
  int64_t find_point(arma_vec point);
  int64_t find_point(arma_vec point, qtnode node);

  /**********************************************************************
     \brief finds the root node of the given point, mainly used to figure
            out which side of a cubesphere the point is on.
     \param point the x, y, z normalized coordinate of the point.
   **/
  int64_t find_root(arma_vec point);

  /**********************************************************************
     \brief If the point is outside of the normalized limits of the 
            quadtree, this tries to put the point back into the domain
     \param point the x, y, z normalized coordinate of the point.
   **/
  arma_vec wrap_point_sphere(arma_vec point);

  /**********************************************************************
     \brief If the point is outside of the normalized limits of the 
            quadtree, this tries to put the point back into the domain
     \param point the x, y, z normalized coordinate of the point.
   **/
  arma_vec wrap_point_cubesphere(arma_vec point);
  
  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  
  bool is_ok();

private:

  /// Defines whether the quadtree state is ok:
  bool IsOk = true;
  /// Defines whether the quadtree is a sphere or not:
  bool IsSphere = false;
  /// Defines whether the quadtree is a cubesphere or not:
  bool IsCubeSphere = false;
  
};

#endif  // INCLUDE_QUADTREE_H_
