#include "functions.hxx"
#include <math.h>

void Functions::RescaleLaplacian(std::vector<std::vector<double>>& laplacien){
  double scale_factor = 2./2.;
  int size = laplacien.size();
  for (int row=0; row<size; row++){
    for (int col=0; col<size; col++){
      laplacien[row][col] *= scale_factor;
    };
    laplacien[row][row] -= 1.;
  };
};

/**
 * @brief Compute Chebyshev coefficients of a kernal.
 * @param kernal The kernal to compute the coefficients of.
 * @param max_coefficients The number of coefficients to compute.
 * @param grid_order The order of the grid to use. If -1, then max_coefficients+1 is used.
 * @param approx_interval_min lower bound of the interval of approximation.
 * @param approx_interval_max upper bound of the interval of approximation.
 * @return The Chebyshev coefficients of the kernal.
 **/
std::vector<double> Functions::ChebyshevCoefficients(double (*kernal)(const double&),
                                          const int& max_coefficients, const int& grid_order=-1,
                                          const double& approx_interval_min=-1.0, const double& approx_interval_max=1.0){
  if (grid_order == -1){
    grid_order = max_coefficients + 1;
  };
  double half_interval = (approx_interval_max - approx_interval_min)/2.;
  double center = (approx_interval_max + approx_interval_min)/2.;
  double grid_value, kernal_value;
  std::vector<double> grid;
  std::vector<double> kernal_values;
  double pi_over_grid_order = math::pi/grid_order;
  for (int i=1; i<grid_order+1; i++){
    grid_value = (i - 0.5) * pi_over_grid_order;
    grid.push_back(grid_value);
    kernal_value = kernal(half_interval * math::cos(grid_value) + center);
    kernal_values.push_back(kernal_value);
  };
  double coefficient_value;
  std::vector<double> chebyshev_coefficients;
  double two_over_grid_order = 2./grid_order;
  for (int i=0; i<grid_order+1; i++){
    coefficient_value = 0.;
    for (int j=0; j<grid_order; j++){
      coefficient_value += kernal_values[j] * math::cos(grid[i] * j) * two_over_grid_order;
    };
    chebyshev_coefficients.push_back(coefficient_value);
  };

  return chebyshev_coefficients;
};
    
std::vector<double> Functions::ChebyshevCoefficients(const int& max_coefficients, const int& grid_order=-1,
                                                     const double& approx_interval_min=-1.0,
                                                     const double& approx_interval_max=1.0){
  return ChebyshevCoefficients([](const double& x){return math::exp(-x);}, max_coefficients, grid_order,
                               approx_interval_min, approx_interval_max);
};


/**
 * @brief Using Chebyshev polynomials, find points with similar eigenvector components to the center_idx.
 * @param laplacian The laplacian.
 * @param chebyshev_coefficients The Chebyshev coefficients.
 * @param center_idx The index to center the approximation around.
 * @param interval The interval to approximate in.
 * @return The approximation of the eigenvalues.
 **/
std::vector<double> LaplacianWavelet(const std::vector<std::vector<double>> &laplacian, const std::vector<double> &chebyshev_coefficients,
                                     const int& center_idx, const std::pair<double, double>& interval) const {
  //TODO
}


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

double Functions::GeneralisedKtDistance(const double& pt1, const double& rapidity1, const double& phi1, 
                                        const double& pt2, const double& rapidity2, const double& phi2,
                                        const double& exponent){
  double ca_distance2 = CambridgeAachenDistance2(rapidity1, phi1, rapidity2, phi2);
  return math::pow(math::min(pt1, pt2), exponent) * math::sqrt(ca_distance2);
};

std::vector<std::vector<double>> GeneralisedKtDistanceMatrix(const std::vector<double>& pts,
                                                             const std::vector<double>& rapidities,
                                                             const std::vector<double>& phis,
                                                             const double& exponent){

};

enum JetMetrics {cambridge_aachen, kt, antikt};
/**
 * @brief Distance between two particles in a named jet metric.
 * @param pt1 The transverse momentum of the first particle.
 * @param rapidity1 The rapidity of the first particle.
 * @param phi1 The azimuthal angle of the first particle.
 * @param pt2 The transverse momentum of the second particle.
 * @param rapidity2 The rapidity of the second particle.
 * @param phi2 The azimuthal angle of the second particle.
 * @param the enum value of the metric.
 * @return The distance between the two particles.
 **/
double NamedDistance(const double& pt1, const double& rapidity1, const double& phi1, 
                     const double& pt2, const double& rapidity2, const double& phi2,
                     const JetMetrics& metric) const;

/**
 * @brief Distance matrix between a set of particles in a named jet metric.
 * @param pts The transverse momenta of the particles.
 * @param rapidities The rapidities of the particles.
 * @param phis The azimuthal angles of the particles.
 * @param metric The enum value of the metric.
 * @return The matrix of distances.
 **/
std::vector<std::vector<double>> NamedDistanceMatrix(
    const std::vector<double>& pts, const std::vector<double>& rapidities, const std::vector<double>& phis,
    const JetMetrics& metric) const;


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

std::vector<std::vector<double>> Functions::Laplacian(
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

  //Functions::RescaleLaplacian(laplacien);

  return laplacien;
};

std::vector<double> Functions::PxPyPz(const double& energy, const double& pt,
                                      const double& rapidity, const double& phi){
  double px = pt*math::cos(phi);
  double py = pt*math::sin(phi);
  double rapidity_factor = math::exp(2*rapidity);
  double pz = energy*(1-rapidity_factor)/(1+rapidity_factor);
  return {px, py, pz};
};

std::vector<double> Functions::PtRapPhi(const double& energy, const double& px,
                                        const double& py, const double& pz){
  double pt = math::sqrt(px*px + py*py);
  double phi = math::atan2(py, px);
  double rapidity = 0.5*math::log((energy+pz)/(energy-pz));
  return {pt, rapidity, phi};
};




