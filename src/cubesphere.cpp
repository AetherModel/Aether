#include "aether.h"

arma_cube map_cube2sphere(const arma_mat& cube) {
    const arma_cube sphere(arma::size(cube));

    sphere(0) = cube2sphere(cube.slice(0), cube.slice(1), cube.slice(2));
    sphere(1) = cube2sphere(cube.slice(1), cube.slice(2), cube.slice(1));
    sphere(2) = cube2sphere(cube.slice(2), cube.slice(0), cube.slice(1));

    return sphere;
    }

arma_vec cube2sphere(const arma_vec& x, const arma_vec& y, const arma_vec& z) {
    const arma_vec y2_2 = (y*y) % 2;
    const arma_vec z2_2 = (y*y) % 2;
    return x%(arma::sqrt(arma::ones<arma_vec>-y2_2-z2_2+y2_2%z2_2*(2/3)));
}

precision_t cube2sphere_float(const precision_t x, const precision_t y, const precision_t z) {
    const precision_t y2_2 = y*y/2;
    const precision_t z2_2 = z*z/2;
}
