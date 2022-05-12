// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_QUADTREE_H_
#define INCLUDE_QUADTREE_H_

#include "aether.h"

/*************************************************
 * QTNode
 *
 * conventions:
 * du = down / up
 * lr = left / right


 *************************************************/

class Quadtree {

public:
  
  const uint64_t nLR = 2;
  const uint64_t nDU = 2;

  struct qtnode {

    std::vector<qtnode> children;
  
    uint64_t depth;

    /// normalized coordinate of lower left corner:
    arma_vec lower_left_norm = {0.0, 0.0, 0.0};

    /// normalized sizes:
    arma_vec size_right_norm = {0.0, 0.0, 0.0};
    arma_vec size_up_norm = {0.0, 0.0, 0.0};

    uint64_t iProcNode;
    uint64_t iProc_down;
    uint64_t iProc_up;
    uint64_t iProc_left;
    uint64_t iProc_right;

  };

  int64_t nRootNodes;
  std::vector<struct qtnode> root_nodes;
  int64_t max_depth;

  arma_vec limit_low = {0.0, 0.0, 0.0};
  arma_vec limit_high = {0.0, 0.0, 0.0};
  
  /// Constructor

  Quadtree(Inputs input, Report report);
  void build(Inputs input, Report report); 
  qtnode new_node(arma_vec lower_left_norm_in,
		  arma_vec size_right_norm_in,
		  arma_vec size_up_norm_in,
		  uint64_t &iProc_in_out,
		  uint64_t depth_in);

  arma_vec get_vect(qtnode node, std::string which);
  arma_vec get_vect(std::string which);

  int64_t find_point(arma_vec point);
  int64_t find_point(arma_vec point, qtnode node);

  arma_vec wrap_point_sphere(arma_vec point);
  arma_vec wrap_point_cubesphere(arma_vec point);
  
  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  
  bool is_ok();

private:

  bool IsOk = true;
  bool IsSphere = false;
  bool IsCubeSphere = false;
  
};

#endif  // INCLUDE_QUADTREE_H_
