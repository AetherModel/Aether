// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#include <armadillo>
using namespace arma;


float solver_chemistry(float old_density, float source, float loss, float dt) {

  float new_density;

  // This isn't a great way of doing things probably, but, what the heck:

  // if sources > losses, take an explicit time-step:

  if (source > loss) {
    new_density = old_density + dt * (source - loss);
  } else {

    // take implicit time-step:

    float normalized_loss = loss / (old_density + 1e-6);
    new_density = (old_density + dt * source) / (1.0 + dt * normalized_loss);

  }

  return new_density;
  
}

fcube solver_chemistry_new(fcube density, fcube source, fcube loss, float dt) {

  std::cout << "in solver\n";

  //  loss = loss / (density + 1e-6);
  
//  std::cout << "in solver 0\n";
//  fcube normalized_loss = loss;
//  std::cout << "in solver 1\n";
//  normalized_loss = normalized_loss / (density + 1e-6);
//  std::cout << "in solver 2\n";
//  fcube new_density = (density + dt * source) / (1.0 + dt * normalized_loss);
  std::cout << "in solver 3\n";

  return density;

}

    
