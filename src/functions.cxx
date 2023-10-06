#include "functions.hxx"
#include <math.h>

double Functions::AngularDistance(const double& phi1, const double& phi2){
  // taking the differnce between two angles requires careful treatment
  return math::min(math::abs(phi1 - phi2), 2*math::pi - math::abs(phi1 - phi2));
};


double Functions::CambridgeAachenDistance2(const double& rapidity1, const double& phi1, 
                                           const double& rapidity2, const double& phi2){
  double delta_rapidity = math::abs(rapidity1 - rapidity2);
  // taking the differnce between two angles requires careful treatment
  double delta_phi = Functions::AngularDistance(phi1, phi2);
  return delta_rapidity*delta_rapidity + delta_phi*delta_phi;
};

std::vector<std::vector<double>> Functions::Affinities(
        const std::vector<double>& rapidities, const std::vector<double>& phis,
        const double& sigma){
  std::vector<std::vector<double>> affinities;
  int n_particles = rapidities.size();
  // Lower triangle
  for (int row_n = 0; row_n < n_particles; row_n++){
    std::vector<double> row;
    // Calculate up to the diagonal
    for (int col_n = 0; col_n < row_n; col_n++){
      double distance2 = Functions::CambridgeAachenDistance2(
              rapidities[row_n], phis[row_n], rapidities[col_n], phis[col_n]);
      row.push_back(math::exp(-distance2/sigma));
    };
    // The diagnal is zero
    row.push_back(0.);
    affinities.push_back(row);
  };
  // Fill in the other half of the triangle by copying
  for (int row_n = 0; row_n < n_particles; row_n++){
    std::vector<double> row = affinities[i];
    for (int col_n = row_n+1; col_n < n_particles; col_n++){
      row.push_back(affinities[col_n][row_n]);
    };
  };
  return affinities;
};

std::vector<std::vector<double>> Functions::Laplacien(
        const std::vector<double>& rapidities, const std::vector<double>& phis,
        const double& sigma, const bool& normalised=true){
  int n_particles = rapidities.size();
  std::vector<std::vector<double>> laplacien = Functions::Affinities(rapidities, phis, sigma);
  // Unnormalised Laplacien
  for (int row_n = 0; row_n < n_particles; row_n++){
    double sum = 0.;
    for (int col_n = 0; col_n < n_particles; col_n++){
      sum += laplacien[row_n][col_n];
      laplacien[row_n][col_n] *= -1;
    };
    laplacien[row_n][row_n] = sum;
  };

  if (normalised){
    // Normalised Laplacien
    double diagonal, inv_sqrt_diagonal;
    for (int row_n = 0; row_n < n_particles; row_n++){
      double diagonal = laplacien[row_n][row_n];
      if (diagonal == 0.){
        inv_sqrt_diagonal = 0.;
      } else {
        inv_sqrt_diagonal = 1./math::sqrt(diagonal);
      };
      for (int col_n = 0; col_n < n_particles; col_n++){
        laplacien[row_n][col_n] *= inv_sqrt_diagonal;
        laplacien[col_n][row_n] *= inv_sqrt_diagonal;
      };
    };
  };

  // TODO
  double max_eigenvalue = Functions::MaxEigenvalue(laplacien);
  // TODO
  Functions::RescaleLaplacian(laplacien, max_eigenvalue);

  return laplacien;
};
