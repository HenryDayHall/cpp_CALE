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


static std::vector<double> Functions::VectorAddition(const std::vector<double> vector1, const std::vector<double> vector2){
  int size = vector1.size();
  std::vector<double> result(size, 0.);
  for (int i=0; i<size; i++){
    result[i] = vector1[i] + vector2[i];
  };
  return result;
};

static std::vector<double> Functions::VectorAddition(const double& factor1, const std::vector<double> vector1,
                                                     const double& factor2, const std::vector<double> vector2){
  int size = vector1.size();
  std::vector<double> result(size, 0.);
  for (int i=0; i<size; i++){
    result[i] = factor1*vector1[i] + factor2*vector2[i];
  };
  return result;
};


static void Functions::VectorAdditionInPlace(std::vector<double> vector1, const std::vector<double> vector2){
  int size = vector1.size();
  for (int i=0; i<size; i++){
    vector1[i] += vector2[i];
  };
};

static void Functions::VectorAdditionInPlace(const double& factor1, std::vector<double> vector1,
                                             const double& factor2, const std::vector<double> vector2){
  int size = vector1.size();
  for (int i=0; i<size; i++){
    vector1[i] *= factor1;
    vector1[i] += factor2 * vector2[i];
  };
};


void Functions::RescaleMatrixInPlace(const double& multipler, std::vector<std::vector<double>>& matrix){
  int n_rows = matrix.size();
  for (int row=0; row<n_rows; row++){
    int n_cols = matrix[row].size();
    for (int col=0; col<n_cols; col++){
      matrix[row][col] *= multipler;
    };
  };
};

static std::vector<std::vector<double>> Functions::RescaleMatrix(const double& multipler, const std::vector<std::vector<double>>& matrix){
  int n_rows = matrix.size();
  std::vector<std::vector<double>> result(n_rows, std::vector<double>());
  for (int row=0; row<n_rows; row++){
    int n_cols = matrix[row].size();
    for (int col=0; col<n_cols; col++){
      result[row].push_back(matrix[row][col] * multipler);
    };
  };
  return result;
};

static std::vector<double> Functions::MatrixDotVector(const std::vector<std::vector<double>>& matrix,
                                                      const std::vector<double>& vector){
  int size = matrix.size();
  std::vector<double> result(size, 0.);
  for (int row=0; row<size; row++){
    for (int col=0; col<size; col++){
      result[row] += matrix[row][col] * vector[col];
    };
  };
  return result;
};

std::vector<double> Functions::LaplacianWavelet(const std::vector<std::vector<double>> &laplacian, const std::vector<double> &chebyshev_coefficients,
                                     const int& center_idx, const std::pair<double, double>& interval) const {
  int n_rows = laplacian.size();
  // An all zero laplacian is a special case
  bool all_zero = true;
  for (int row = 0; row < n_rows; row++) {
    for (int col = 0; col < n_rows; col++) {
      if (laplacian[row][col] != 0.) {
        all_zero = false;
        break;
      };
    };
    if (!all_zero) break;
  };
  if (all_zero) {
    std::vector<double> result(n_rows, 0.);
    return result;
  };

  // Note, while this takes the place of make_wavelets, it's mostly cheby_op, as the starting point is the chebyshev coefficients.
  // basically equation 66 of https://inria.hal.science/hal-01943589/document
  int n_coeffients = chebyshev_coefficients.size();
  double half_interval = (interval.second - interval.first)/2.;
  double inverse_half_interval = 1./half_interval;
  double inverse_interval = 2./half_interval;
  double center = (interval.second + interval.first)/2.;
  // We will need two vectors for taking repeated dot products onto the laplacian
  // Last iteration
  std::vector<double> fourier_transform_old(n_rows, 0.);
  fourier_transform_old[center_idx] = 1.;
  // Current iteration
  std::vector<double> fourier_transform_new = VectorAddition(MatrixDotVector(Laplacian, fourier_transform_old),
                                                             RescaleMatrix(-center, fourier_transform_old))
  RescaleMatrixInPlace(inverse_half_interval, fourier_transform_new);
  // Placeholder for swapping them
  std::vector<double> fourier_transform_placeholder(n_rows, 0.);

  // We also need a place to store the growing sum
  std::vector<double> results = VectorAddition(RescaleMatrix(0.5*chebyshev_coefficients[0], fourier_transform_old),
                                               RescaleMatrix(chebyshev_coefficients[1], fourier_transform_new));

  for (int k=2; k < n_coeffients; k++){
    // update the old fourier transform to make a new one
    VectorAdditionInPlace(-1, fourier_transform_old, inverse_interval,
                          VectorAddition(1, MatrixDotVector(Laplacian, fourier_transform_new),
                                         -center, fourier_transform_new));
    // switch the old and new fourier transforms
    fourier_transform_placeholder = fourier_transform_old;
    fourier_transform_old = fourier_transform_new;
    fourier_transform_new = fourier_transform_placeholder;
    // update the results
    VectorAdditionInPlace(1, results, chebyshev_coefficients[k], fourier_transform_new);
  };

  return results;
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




