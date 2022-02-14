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

CubeSphere::QTNode::QTNode(const uint& face,
                           const arma_vec& position,
                           const uint& resolution,
                           const uint& depth){
    m_position = position;
    m_depth = depth;

    if (resolution != depth) {
        // m_children = std::make_shared<QTNode>(position-
        }
    }

CubeSphere::QTNode::~QTNode(){
    // Place holder just in case.
}

bool CubeSphere::QTNode::is_leaf(){
    return m_resolution == m_depth;
}
