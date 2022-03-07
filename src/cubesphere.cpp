#include "aether.h"
#include <cmath>

arma_vec CubeSphere::cube2sphere(const arma_vec& x, const arma_vec& y, const arma_vec& z) {
    const arma_vec y2_2 = (y%y) / 2;
    const arma_vec z2_2 = (y%y) / 2;
    const arma_vec one = arma::ones<arma_vec>(x.size());
    return x%(arma::sqrt(one-y2_2-z2_2+y2_2%z2_2*(2/3)));
    }

precision_t CubeSphere::cube2sphere(const precision_t& x, const precision_t& y, const precision_t& z) {
    const precision_t y2_2 = y*y/2;
    const precision_t z2_2 = z*z/2;
    return x*(sqrt(1-y2_2-z2_2+y2_2*z2_2*(2/3)));
    }

arma_vec CubeSphere::cube2sphere(const arma_vec& p){
    const arma_vec p2 = p%p;
    arma_vec ret_val = {0, 0, 0};
    static const arma_mat others = {
        {0, 1, 1},
        {1, 0, 1},
        {0, 1, 1}
        };
    const arma_mat last_term = {
        {0, 0, p2(1)},
        {p2(2), 0, 0},
        {0, p2(0), 0}
        };
    const arma_vec one = arma::ones<arma_vec>(3);
    ret_val = p%arma::sqrt(one-0.5*others*p2+last_term/3.0);

    /* ret_val(0) = p(0)*sqrt(1.0-0.5*(p2(1)+p2(2))+(1.0/3.0)*(p2(1)*p2(2))); */
    /* ret_val(1) = p(1)*sqrt(1.0-0.5*(p2(2)+p2(0))+(1.0/3.0)*(p2(2)*p2(0))); */
    /* ret_val(2) = p(2)*sqrt(1.0-0.5*(p2(0)+p2(1))+(1.0/3.0)*(p2(0)*p2(1))); */

    return ret_val;
    }

CubeSphere::QTNode::QTNode(Grid grid,
                           const uint& face,
                           const arma_vec& cube_position,
                           const uint& depth,
                           const uint& max_depth){
    m_cube_position = cube_position;
    m_depth = depth;

    // refine further
    if (max_depth != depth) {
        const arma_vec right = CubeSphere::RIGHTS.row(face);
        const arma_vec up = CubeSphere::UPS.row(face);

        // refine the quad tree
        uint child = 0;
        for (int j=0; j<2; j++) {
            for (int i=0; i<2; i++) {
                // new corner position
                const arma_vec position = m_cube_position + 2.0*(right*i + up*j)*pow(0.5, depth);
                m_children[child] = std::shared_ptr<QTNode>(
                        new QTNode(grid, face, position, depth+1, max_depth));
                child++;
            }
        }

    } else {
        // create grid
    }
}

CubeSphere::QTNode::~QTNode(){
    // Place holder just in case.
}

bool CubeSphere::QTNode::is_leaf(){
    for (int i=0; i<4; i++) {
        if (m_children[i] != nullptr) { return false; }
    }
    return true;
}
