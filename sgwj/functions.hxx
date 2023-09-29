#ifndef FUNCTIONS_HXX
#define FUNCTIONS_HXX

#include <iostream>
#include <math>
#include <vector>

// Making a class of this, in case future modifications require some
// clever use of existing calculations.
class Functions
{
  public:

    /**
     * @brief Calculates the maximum eigenvalue of a symmetrix matrix.
     * @param symmetrix_matrix The symmetrix matrix.
     */
    double MaxEigenvalue(const std::vector<std::vector<double>>& symmetrix_matrix);

    /**
     * @brief Approximate the maximum eigenvalue of a symmetrix matrix.
     * @param symmetrix_matrix The symmetrix matrix.
     */
    double RoughMaxEigenvalue(const std::vector<std::vector<double>>& symmetrix_matrix);

    /**
     * @brief Rescale eigenvalues of the laplacien to -1 to 1.
     * @param laplacien The laplacien.
     */
    void RescaleLaplacian(std::vector<std::vector<double>>& laplacien);


    /**
     * @brief Compute a set of wavelet scales for spectrum bound.
     *
     * Scales are spaced logarithmically between the mimumum and
     * maximum effective scales.
     * @param min_eigenvalue The minimum non-zero eigenvalue of the laplacien.
     * @param max_eigenvalue The maximum eigenvalue of the laplacien.
     * @param number_of_scales The number of scales to compute.
     */
    std::vector<double> ChoseScales(const double& min_eigenvalue, const double& max_eigenvalue, const int& number_of_scales);

    /**
     * @brief Compute the Mh exponential kernal at inputs.
     * @param input The location at which to compute.
     **/
    std::vector<double> MhKernal(const std::vector<double>& input);

    /**
     * @brief Compute the Abspline kernal at inputs.
     * @param input The location at which to compute.
     **/
    std::vector<double> AbsplineKernal(const std::vector<double>& input, const double& param_a, const double& param_b,
                                       const double& param_t1, const double& param_t2);

    /**
     * @brief types of kernal avaliable.
     **/
    enum KernalType {abspline, mh};
        
    /**
     * @brief Compute Chebyshev coefficients of a kernal.
     * @param kernal The kernal to compute the coefficients of.
     * @param max_coefficients The number of coefficients to compute.
     * @param grid_order The order of the grid to use. If -1, then max_coefficients+1 is used.
     **/
    std::vector<double> ChebyshevCoefficients(std::vector<double> (*kernal)(const std::vector<double>&),
                                              const int& max_coefficients, const int& grid_order=-1,
                                              const double& approx_interval_min=-1.0, const double& approx_interval_max=1.0);

    /**
     * @brief Distance between to angles.
     * @param phi1 The first angle.
     * @param phi2 The second angle.
     **/
    double AngularDistance(const double& phi1, const double& phi2);

    /**
     * @brief Distance between two particles in the Cambridge-Aachen metric.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     **/
    double CambridgeAachenDistance2(const double& rapidity1, const double& phi1, 
                                    const double& rapidity2, const double& phi2);

    /**
     * @brief Graph laplacien for a set of particles.
     * @param rapidities The rapidities of the particles.
     * @param phis The azimuthal angles of the particles.
     * @param normalised Whether to normalise the laplacien.
     **/
    double Laplacien(const std::vector<double>& rapidities, const std::vector<double>& phis,
                     const bool& normalised=true);

};

#endif // FUNCTIONS_HXX
