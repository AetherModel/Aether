// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// Want to be able to call the calculations with fvecs and vectors,
// but only want to do calculation in one function:

fvec calc_bin_edges(fvec centers) {
  std::vector<float> centers_vec = make_vector_from_fvec(centers);
  std::vector<float> edges_vec = calc_bin_edges(centers_vec);
  fvec edges = make_fvec_from_vector(edges_vec);
  return edges;
}

std::vector<float> calc_bin_edges(std::vector<float> centers) {

  std::vector<float> edges;
  float dc;
  int64_t nPts = centers.size();

  if (nPts < 2) {
    edges.push_back(-1);
    return edges;
  }

  dc = centers[1] - centers[0];
  // edges[0] is to the left of centers[0]
  // edges[1] is between centers[0] and centers[1]
  edges.push_back(centers[0] - dc / 2);
  edges.push_back(centers[0] + dc / 2);

  // This assumes that:
  //    c1 - e1 = e2 - c1 => e2 = 2 * c1 - e1
  //    c2 - e2 = e3 - c2 -> e3 = 2 * c2 - e2
  //    etc...
  for (int64_t iPt = 1; iPt < nPts; iPt++)
    edges.push_back(2 * centers[iPt] - edges[iPt]);

  return edges;
}

fvec calc_bin_widths(fvec centers) {
  std::vector<float> centers_vec = make_vector_from_fvec(centers);
  std::vector<float> widths_vec = calc_bin_widths(centers_vec);
  fvec widths = make_fvec_from_vector(widths_vec);
  return widths;
}

std::vector<float> calc_bin_widths(std::vector<float> centers) {

  std::vector<float> edges = calc_bin_edges(centers);
  std::vector<float> widths;

  int64_t nPts = edges.size();

  // Widths are simply the differences between the edges:
  for (int64_t iPt = 1; iPt < nPts; iPt++)
    widths.push_back(edges[iPt] - edges[iPt - 1]);

  return widths;
}
