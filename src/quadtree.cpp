
#include "aether.h"

int64_t iProcQuery = -1;

Quadtree::Quadtree() {
  if (input.get_is_cubesphere())
    nRootNodes = 6;
  else
    nRootNodes = 1;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Quadtree::is_ok() {
  return IsOk;
}

// --------------------------------------------------------------------------
// build quadtree
// --------------------------------------------------------------------------

void Quadtree::build() {

  arma_mat origins;
  arma_mat rights;
  arma_mat ups;

  if (input.get_is_cubesphere()) {
    origins = CubeSphere::ORIGINS;
    rights = CubeSphere::RIGHTS;
    ups = CubeSphere::UPS;
    IsCubeSphere = true;
  } else {
    origins = Sphere::ORIGINS;
    rights = Sphere::RIGHTS;
    ups = Sphere::UPS;
    IsSphere = true;
  }

  arma_vec o(3), r(3), u(3);

  // This captures the limits of the sphere, independent of what the
  // user asks for
  for (int64_t iNode = 0; iNode < nRootNodes; iNode++) {
    for (int i = 0; i < 3; i++) {
      o(i) = origins(iNode, i);
      r(i) = rights(iNode, i);
      u(i) = ups(iNode, i);

      if (o(i) < limit_low(i))
        limit_low(i) = o(i);

      if (o(i) + r(i) > limit_high(i))
        limit_high(i) = o(i) + r(i);

      if (o(i) + u(i) > limit_high(i))
        limit_high(i) = o(i) + u(i);
    }
  }

  // Before we build the quadtree, we need to allow the user to
  // restrict the domain.  This will only work for the spherical
  // grid so far:

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  if (grid_input.lon_min > 0.0 ||
      grid_input.lon_max < 2.0 * cPI ||
      grid_input.lat_min > -cPI / 2.0 ||
      grid_input.lat_max < cPI / 2.0) {
    // We are dealing with less than the whole Earth...
    origins(0) = grid_input.lon_min / cPI;
    origins(1) = grid_input.lat_min / cPI;
    rights(0) = (grid_input.lon_max - grid_input.lon_min) / cPI;
    ups(1) = (grid_input.lat_max - grid_input.lat_min) / cPI;
  }

  qtnode tmp;

  for (uint64_t iNode = 0; iNode < nRootNodes; iNode++) {
    if (report.test_verbose(2))
      std::cout << "Making quadtree node : " << iNode << "\n";

    uint64_t iP = iNode * pow(4, max_depth);
    uint64_t iDepth = 0;

    if (report.test_verbose(2))
      std::cout << "  iProcessor Start : " << iP << "\n";

    for (int i = 0; i < 3; i++) {
      o(i) = origins(iNode, i);
      r(i) = rights(iNode, i);
      u(i) = ups(iNode, i);
    }

    tmp = new_node(o, r, u, iP, iDepth, iNode);
    root_nodes.push_back(tmp);
  }
}

// --------------------------------------------------------------------------
// make a new node in the quadtree
//   iSide_in indicates the root node, since for cubesphere, these are same
// --------------------------------------------------------------------------

Quadtree::qtnode Quadtree::new_node(arma_vec lower_left_norm_in,
                                    arma_vec size_right_norm_in,
                                    arma_vec size_up_norm_in,
                                    uint64_t &iProc_in_out,
                                    uint64_t depth_in,
                                    uint64_t iSide_in) {

  qtnode tmp, tmp_child;

  tmp.lower_left_norm = lower_left_norm_in;
  tmp.size_right_norm = size_right_norm_in;
  tmp.size_up_norm = size_up_norm_in;
  tmp.depth = depth_in;
  tmp.iSide = iSide_in;

  // refine further
  if (tmp.depth < max_depth) {

    arma_vec size_right_norm_child = tmp.size_right_norm / nLR;
    arma_vec size_up_norm_child = tmp.size_up_norm / nDU;

    // refine the quad tree
    uint64_t iChild = 0;

    for (int iDU = 0; iDU < nDU; iDU++) {
      for (int iLR = 0; iLR < nLR; iLR++) {

        arma_vec lower_left_norm_child =
          tmp.lower_left_norm +
          iLR * size_right_norm_child +
          iDU * size_up_norm_child;

        tmp_child = new_node(lower_left_norm_child,
                             size_right_norm_child,
                             size_up_norm_child,
                             iProc_in_out,
                             tmp.depth + 1,
                             iSide_in);

        tmp.children.push_back(tmp_child);

        if (iProc == iProcQuery) {
          std::cout << "      depth_in : " << depth_in << " "
                    << max_depth << " "
                    << iDU << " "
                    << iLR << " "
                    << iSide_in << " "
                    << tmp_child.iProcNode << "\n";
        }

        iChild++;
      }  // nDU
    }  // nLR
  }  else {
    tmp.iProcNode = iProc_in_out;

    if (iGrid == iProc_in_out)
      iSide = iSide_in;

    iProc_in_out++;
  }

  return tmp;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

arma_vec Quadtree::get_vect(Quadtree::qtnode node, std::string which) {

  arma_vec vect(3, fill::value(-1e32));

  if (node.children.size() > 0) {
    for (uint64_t iChild = 0; iChild < 4; iChild++) {
      vect = get_vect(node.children[iChild], which);

      if (vect(0) > -1e32)
        break;
    }
  } else {
    if (node.iProcNode == iGrid) {
      if (which == "LL")
        vect = node.lower_left_norm;

      if (which == "SR")
        vect = node.size_right_norm;

      if (which == "SU")
        vect = node.size_up_norm;

      if (which == "MID")
        vect = node.lower_left_norm + node.size_right_norm / 2 + node.size_up_norm / 2;
    }
  }

  return vect;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

arma_vec Quadtree::get_vect(std::string which) {

  arma_vec vect(3, fill::value(-1e32));

  for (int64_t iNode = 0; iNode < nRootNodes; iNode++) {
    vect = get_vect(root_nodes[iNode], which);

    if (vect(0) > -1e32)
      break;
  }

  return vect;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

int64_t Quadtree::find_point(arma_vec point, Quadtree::qtnode node) {

  int64_t iNode = -1;

  if (iProc == iProcQuery) {
    std::cout << "find_point - depth : " << node.depth << " "
              << node.children.size() << "\n";
    std::cout << "point : ";
    display_vector(point);
  }

  if (node.children.size() > 0) {
    for (uint64_t iChild = 0; iChild < 4; iChild++) {
      iNode = find_point(point, node.children[iChild]);

      if (iProc == iProcQuery)
        std::cout << "iNode : " << iNode << " " << iChild << "\n";

      if (iNode > -1)
        break;
    }
  } else {
    int iFound = 0;

    int iLR_ = -1;
    int iDU_ = -1;
    int iCO_ = -1;

    for (int i = 0; i < 3; i++) {
      if (node.size_right_norm(i) == 0 &&
          node.size_up_norm(i) == 0) {
        iCO_ = i;

        if (iProc == iProcQuery)
          std::cout << "found co : " << iCO_ << "\n";
      }

      if (node.size_right_norm(i) != 0 &&
          node.size_up_norm(i) == 0) {
        iLR_ = i;

        if (iProc == iProcQuery)
          std::cout << "found lr : " << iLR_ << "\n";
      }

      if (node.size_right_norm(i) == 0 &&
          node.size_up_norm(i) != 0) {
        iDU_ = i;

        if (iProc == iProcQuery)
          std::cout << "found du : " << iDU_ << "\n";
      }
    }

    if (fabs(node.lower_left_norm(iCO_) - point(iCO_)) < 1.0e-6)
      iFound++;

    if ((point(iLR_) >= node.lower_left_norm(iLR_) &&
         point(iLR_) < node.lower_left_norm(iLR_) + node.size_right_norm(iLR_)) ||
        (point(iLR_) >= node.lower_left_norm(iLR_) + node.size_right_norm(iLR_) &&
         point(iLR_) < node.lower_left_norm(iLR_)))
      iFound++;

    if ((point(iDU_) >= node.lower_left_norm(iDU_) &&
         point(iDU_) < node.lower_left_norm(iDU_) + node.size_up_norm(iDU_)) ||
        (point(iDU_) >= node.lower_left_norm(iDU_) + node.size_up_norm(iDU_) &&
         point(iDU_) < node.lower_left_norm(iDU_)))
      iFound++;

    if (iFound == 3)
      iNode = node.iProcNode;

    if (iProc == iProcQuery && iFound == 3) {
      std::cout << "Found on iNode : " << iNode << "\n   point     : ";
      display_vector(point);
      std::cout << "   lower left : ";
      display_vector(node.lower_left_norm);
      std::cout << "   size_right : ";
      display_vector(node.size_right_norm);
      std::cout << "   size_up    : ";
      display_vector(node.size_up_norm);
    }

  }

  return iNode;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

arma_vec Quadtree::wrap_point_sphere(arma_vec point) {

  arma_vec wrap_point = point;

  int64_t iEdge_ = -1;

  for (int i = 0; i < 3; i++) {
    if (point(i) < limit_low(i) ||
        point(i) > limit_high(i)) {
      iEdge_ = i;
      break;
    }
  }

  if (iEdge_ > -1) {
    if (iProc == iProcQuery) {
      std::cout << " point out of bounds!  Wrapping : " << point << "\n";
      std::cout << " edge : " << iEdge_ << "\n";
    }

    if (iEdge_ == 0) {
      // E/W wrap
      if (point(0) < limit_low(0))
        wrap_point(0) = limit_high(0) + point(0);
      else
        wrap_point(0) = limit_low(0) + (point(0) - limit_high(0));
    }

    if (iEdge_ == 1) {
      // N/S connection
      if (wrap_point(1) < limit_low(1))
        wrap_point(1) = 2 * limit_low(1) - point(1);
      else
        wrap_point(1) = 2 * limit_high(1) - point(1);

      wrap_point(0) = point(0) + (limit_high(0) - limit_low(0)) / 2.0;

      if (wrap_point(0) >= limit_high(0))
        wrap_point(0) = wrap_point(0) - limit_high(0);
    }

    if (iProc == iProcQuery)
      std::cout << " wrap_point : " << wrap_point << "\n";
  }

  return wrap_point;

}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

arma_vec Quadtree::wrap_point_cubesphere(arma_vec point) {

  arma_vec wrap_point = point;

  int64_t iComp_ = -1;
  double delta = 0.0, limit = 0.0;

  // First, determine if we go off the current face.
  //  - check the limits in the each direction to see if they are exceeded
  //  - if they are, calculate how far we go off, since we will wrap that much
  //  - store limit - we can use that to determine face to move to
  for (int i = 0; i < 3; i++) {
    if (point(i) < limit_low(i)) {
      iComp_ = i;
      delta = limit_low(i) - point(i);
      wrap_point(i) = limit_low(i);
      limit = limit_low(i);
      break;
    }

    if (point(i) > limit_high(i)) {
      iComp_ = i;
      delta = point(i) - limit_high(i);
      wrap_point(i) = limit_high(i);
      limit = limit_high(i);
      break;
    }
  }

  if (iComp_ > -1) {

    // Try to figure out which face we are going to land on:

    int64_t iFaceFrom_ = iSide;
    int64_t iFaceTo_ = -1;
    int64_t iFace = 0;

    for (iFace = 0; iFace < 6; iFace++) {
      if (CubeSphere::RIGHTS(iFace, iComp_) == 0 &&
          CubeSphere::UPS(iFace, iComp_) == 0 &&
          CubeSphere::ORIGINS(iFace, iComp_) == limit)
        iFaceTo_ = iFace;
    }

    // in theory, the line should now be touching the limits of
    // the face. How to really test this?

    int64_t iCompTo_ = -1;
    double sn = 1.0;

    for (int i = 0; i < 3; i++) {
      if (CubeSphere::RIGHTS(iFaceFrom_, i) == 0 &&
          CubeSphere::UPS(iFaceFrom_, i) == 0) {
        iCompTo_ = i;

        if (CubeSphere::ORIGINS(iFaceFrom_, i) > 0)
          sn = -1.0;
      }
    }

    // move away from the edge now:
    wrap_point(iCompTo_) = wrap_point(iCompTo_) + sn * delta;

    if (iProc == iProcQuery) {
      std::cout << " point out of bounds!  Wrapping : ";
      display_vector(point);
      std::cout << "   delta : " << delta << "\n";
      std::cout << "   limit : " << limit << "\n";
      std::cout << "   face from : " << iFaceFrom_ << "\n";
      std::cout << "   face to : " << iFaceTo_ << "\n";
      std::cout << "   comp from : " << iComp_ << "\n";
      std::cout << "   comp to : " << iCompTo_ << "\n";
      std::cout << "     sign : " << sn << "\n";
      std::cout << "   new point : ";
      display_vector(wrap_point);
    }
  }

  return wrap_point;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

int64_t Quadtree::find_point(arma_vec point) {

  arma_vec wrap_point;

  if (IsSphere)
    wrap_point = wrap_point_sphere(point);

  if (IsCubeSphere)
    wrap_point = wrap_point_cubesphere(point);

  int64_t iNode = -1;

  for (int64_t iRoot = 0; iRoot < nRootNodes; iRoot++) {
    iNode = find_point(wrap_point, root_nodes[iRoot]);

    if (iNode > -1)
      break;
  }

  return iNode;
}

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

int64_t Quadtree::find_root(arma_vec point) {

  arma_vec wrap_point;

  if (IsSphere)
    wrap_point = wrap_point_sphere(point);

  if (IsCubeSphere)
    wrap_point = wrap_point_cubesphere(point);

  int64_t iNode = -1, iRoot;

  for (iRoot = 0; iRoot < nRootNodes; iRoot++) {
    if (iProc == iProcQuery)
      std::cout << "Root node : " << iRoot << " of " << nRootNodes << "\n";

    iNode = find_point(wrap_point, root_nodes[iRoot]);

    if (iNode > -1)
      break;
  }

  if (iRoot == nRootNodes)
    iRoot = -1;

  return iRoot;
}
