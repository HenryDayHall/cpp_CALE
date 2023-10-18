#include <pybind11/pybind11.h>

#include <pybind11/stl.h>
#include "functions.hxx"
#include "cluster.hxx"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

std::vector<std::vector<float>> empty_matrix_for_testing = {};

PYBIND11_MODULE(CALE_pybind, m) {
    m.doc() = R"pbdoc(
        Spectral Graph Wavelet Jets
        -----------------------
        Written in c++, bound in python for your convenience.

        .. currentmodule:: CALE_pybind

        .. autosummary::
           :toctree: _generate

    )pbdoc";


  // Functions ~~~~~~~~~~~~~~~~~~~~~~
    m.def("ChebyshevCoefficients",
    py::overload_cast<const int&, const int&, const double&, const double&>(&Functions::ChebyshevCoefficients),
        R"pbdoc(
     * @brief Compute Chebyshev coefficients of a negative exponent kernal.
     * @param max_coefficients The number of coefficients to compute.
     * @param grid_order The order of the grid to use. If -1, then max_coefficients+1 is used.
     * @param approx_interval_min lower bound of the interval of approximation.
     * @param approx_interval_max upper bound of the interval of approximation.
     * @return The Chebyshev coefficients of the kernal.
     )pbdoc",
        py::arg("max_coefficients"), py::arg("grid_order")=-1, py::arg("approx_interval_min")=-1.0, py::arg("approx_interval_max")=1.0
        );

    m.def("VectorAddition",
        py::overload_cast<const std::vector<double>, const std::vector<double>>(&Functions::VectorAddition),
        R"pbdoc(
     * @brief Sum two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @return vector1 + vector2.
     )pbdoc",
        py::arg("vector1"), py::arg("vector2")
        );

    m.def("VectorAddition",
        py::overload_cast<const double&, const std::vector<double>, const double&, const std::vector<double>>(&Functions::VectorAddition),
        R"pbdoc(
     * @brief Sum two vectors with rescaling.
     * @param factor1 The factor to multiply vector1 by.
     * @param vector1 The first vector.
     * @param factor2 The factor to multiply vector2 by.
     * @param vector2 The second vector.
     * @return vector1 + vector2.
     )pbdoc",
        py::arg("factor1"), py::arg("vector1"), py::arg("factor2"), py::arg("vector2")
        );


    m.def("MatrixDotVector", &Functions::MatrixDotVector, R"pbdoc(
     * @brief Dot product between a matrix and a vector
     * @param matrix The matrix.
     * @param vector The vector.
     * @return The dot product.
     )pbdoc",
        py::arg("matrix"), py::arg("vector")
        );


    m.def("LaplacianWavelet", &Functions::LaplacianWavelet, R"pbdoc(
     * @brief Using Chebyshev polynomials, find points with similar eigenvector components to the center_idx.
     * @param laplacian The laplacian.
     * @param chebyshev_coefficients The Chebyshev coefficients.
     * @param center_idx The index to center the approximation around.
     * @param interval The interval to approximate in.
     * @return The approximation of the eigenvalues.
     )pbdoc",
        py::arg("laplacian"), py::arg("chebyshev_coefficients"), py::arg("center_idx"), py::arg("interval")
        );

    m.def("AngularDistance", &Functions::AngularDistance, R"pbdoc(
     * @brief Distance between to angles.
     * @param phi1 The first angle.
     * @param phi2 The second angle.
     * @return The distance between the two angles.
     )pbdoc",
        py::arg("phi1"), py::arg("phi2")
        );

    m.def("CambridgeAachenDistance2", &Functions::CambridgeAachenDistance2, R"pbdoc(
     * @brief Distance squared between two particles in the Cambridge-Aachen metric.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @return The distance squared between the two particles.
     )pbdoc",
        py::arg("rapidity1"), py::arg("phi1"), py::arg("rapidity2"), py::arg("phi2")
        );


    m.def("GeneralisedKtDistance2", &Functions::GeneralisedKtDistance2, R"pbdoc(
     * @brief Distance squared between two particles in the Generalised-Kt metric.
     * @param pt1 The transverse momentum of the first particle.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param pt2 The transverse momentum of the second particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @param the exponent of the metric.
     * @return The distance squared between the two particles.
     )pbdoc",
        py::arg("pt1"), py::arg("rapidity1"), py::arg("phi1"), py::arg("pt2"), py::arg("rapidity2"), py::arg("phi2"), py::arg("exponent")
        );


    m.def("GeneralisedKtDistance2Matrix", &Functions::GeneralisedKtDistance2Matrix, R"pbdoc(
     * @brief Distance squared matrix between a set of particles in the Generalised-Kt metric.
     * Note the distance is squared, not the matrix.
     * @param pts The transverse momenta of the particles.
     * @param rapidities The rapidities of the particles.
     * @param phis The azimuthal angles of the particles.
     * @param the exponent of the metric.
     * @return The distance squared matrix between the particles.
     )pbdoc",
        py::arg("pts"), py::arg("rapidities"), py::arg("phis"), py::arg("exponent")
        );

    py::enum_<Functions::JetMetrics>(m, "JetMetrics")
        .value("cambridge_aachen", Functions::JetMetrics::cambridge_aachen)
        .value("kt", Functions::JetMetrics::kt)
        .value("antikt", Functions::JetMetrics::antikt);

    m.def("NamedDistance2", &Functions::NamedDistance2, R"pbdoc(
     * @brief Distance squared between two particles in a named jet metric.
     * @param pt1 The transverse momentum of the first particle.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param pt2 The transverse momentum of the second particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @param the enum value of the metric.
     * @return The distance squared between the two particles.
     )pbdoc",
        py::arg("pt1"), py::arg("rapidity1"), py::arg("phi1"), py::arg("pt2"), py::arg("rapidity2"), py::arg("phi2"), py::arg("metric")
        );

    m.def("NamedDistance2Matrix", &Functions::NamedDistance2Matrix, R"pbdoc(
     * @brief Distance squared matrix between a set of particles in a named jet metric.
     * Note the distances are squared, not the matrix.
     * @param pts The transverse momenta of the particles.
     * @param rapidities The rapidities of the particles.
     * @param phis The azimuthal angles of the particles.
     * @param metric The enum value of the metric.
     * @return The matrix of squared distances.
     )pbdoc",
        py::arg("pts"), py::arg("rapidities"), py::arg("phis"), py::arg("metric")
        );

    m.def("Affinities", &Functions::Affinities, R"pbdoc(
     * @brief Affinity matrix for a set of particles.
     * @param distances2 The distances squared in the chosen metric between particles.
     * @param sigma The sigma in the denominator of the exponent of the affinity.
     * @return The matrix of affinites.
     )pbdoc",
        py::arg("distances2"), py::arg("sigma")
        );

    m.def("Laplacian", &Functions::Laplacian, R"pbdoc(
     * @brief Graph laplacien for a set of particles.
     * @param distances2 The distances squared in the chosen metric between particles.
     * @param sigma The sigma in the denominator of the exponent of the affinity.
     * @param normalised Whether to normalise the laplacien.
     * @return The laplacien.
     )pbdoc",
        py::arg("distances2"), py::arg("sigma"), py::arg("normalised")=true
        );


    m.def("PxPyPz", &Functions::PxPyPz, R"pbdoc(
     * @brief From detector coordinates, calculate cartesien coordinates
     * @param energy The energy.
     * @param pt The transverse momentum.
     * @param rapidity The rapidity.
     * @param phi The azimuthal angle.
     * @return The px, py, pz.
    )pbdoc",
        py::arg("energy"), py::arg("pt"), py::arg("rapidity"), py::arg("phi")
        );

    m.def("PtRapPhi", &Functions::PtRapPhi, R"pbdoc(
     * @brief From cartesien coordinates, calculate detector coordinates
     * @param energy The energy.
     * @param px The x momentum.
     * @param py The y momentum.
     * @param pz The z momentum.
     * @return The pt, rapidity, phi.
     )pbdoc",
        py::arg("energy"), py::arg("px"), py::arg("py"), py::arg("pz")
        );

  // Cluster ~~~~~~~~~~~~~~~~~

    py::class_<Cluster>(m, "Cluster")
        .def(py::init<const double&, const double&, const int&>(),
            R"pbdoc(
     * @brief Constructor for a clustering algorithm.
     *
     * @param sigma
     * @param cutoff
     * @param n_rounds
     )pbdoc",
            py::arg("sigma"), py::arg("cutoff"), py::arg("n_rounds")
            )
    .def("SetInputs", &Cluster::SetInputs,
        R"pbdoc(
     * @brief Chose a new set of particles to start clustering.   
     *
     * @param labels labels of the particles to be clustered
     * @param energies energies of the particles to be clustered
     * @param pts transverse momenta of the particles to be clustered
     * @param rapidites rapidities of the particles to be clustered
     * @param phis azimuthal angles of the particles to be clustered
     )pbdoc",
        py::arg("labels"),
        py::arg("energies"),
        py::arg("pts"),
        py::arg("rapidites"),
        py::arg("phis")
        )
    .def("GetLabels", &Cluster::GetLabels,
        R"pbdoc(
     * @brief Return the labels of the current pseudojets.
     * @return vector of labels
     )pbdoc",
        py::return_value_policy::reference_internal
        )
    .def("GetEnergies", py::overload_cast<>(&Cluster::GetEnergies, py::const_),
        R"pbdoc(
     * @brief Return the energies of the current pseudojets.
     * @return vector of energies
     )pbdoc",
        py::return_value_policy::reference_internal
        )
    .def("GetEnergies", py::overload_cast<const std::vector<int>&>(&Cluster::GetEnergies, py::const_),
        R"pbdoc(
     * @brief Return the energies of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of energies
     )pbdoc",
       py::arg("labels")
        ) 
    .def("GetPts", py::overload_cast<>(&Cluster::GetPts, py::const_),
        R"pbdoc(
     * @brief Return the transverse momenta of the current pseudojets.
     * @return vector of transverse momenta
     )pbdoc",
       py::return_value_policy::reference_internal
        ) 
    .def("GetPts", py::overload_cast<const std::vector<int>&>(&Cluster::GetPts, py::const_),
        R"pbdoc(
     * @brief Return the transverse momenta of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of transverse momenta
     )pbdoc",
       py::arg("labels")
        ) 
    .def("GetRapidites", py::overload_cast<>(&Cluster::GetRapidites, py::const_),
        R"pbdoc(
     * @brief Return the rapidities of the current pseudojets.
     * @return vector of rapidities
     )pbdoc",
       py::return_value_policy::reference_internal
        ) 
    .def("GetRapidites", py::overload_cast<const std::vector<int>&>(&Cluster::GetRapidites, py::const_),
        R"pbdoc(
     * @brief Return the rapidities of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of rapidities
     )pbdoc",
       py::arg("labels")
        ) 
    .def("GetPhis", py::overload_cast<>(&Cluster::GetPhis, py::const_),
        R"pbdoc(
     * @brief Return the azimuthal angles of the current pseudojets.
     * @return vector of azimuthal angles
     )pbdoc",
       py::return_value_policy::reference_internal
        ) 
    .def("GetPhis", py::overload_cast<const std::vector<int>&>(&Cluster::GetPhis, py::const_),
        R"pbdoc(
     * @brief Return the azimuthal angles of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of azimuthal angles
     )pbdoc",
       py::arg("labels")
        ) 
    .def("GetNextMerge", &Cluster::GetNextMerge,
        R"pbdoc(
     * @brief Say which particles form the next jet.
     * Returned as a list of labels.
     * @return vector of labels
     )pbdoc")
    .def("DoMerge", py::overload_cast<std::vector<int>>(&Cluster::DoMerge),
        R"pbdoc(
     * @brief Merge a set of pesudojets into a complete jet.
     * Calculate the kinematics of the new jet internally, and
     * give it the next free label.
     *
     * @param labels labels of the pseudojets to merge
     )pbdoc",
       py::arg("labels")
        ) 
    .def("DoMerge", py::overload_cast<std::vector<int>, int, double,
        double, double, double>(&Cluster::DoMerge),
        R"pbdoc(
     * @brief Merge a set of pesudojets into a complete jet.
     * Merge the two jets, and assign the combined jet externally
     * provided label and kinematics.
     *
     * @param labels labels of the pseudojets to merge
     * @param label_new label of the new jet
     * @param energy energy of the new jet
     * @param pt transverse momentum of the new jet
     * @param rapidity rapidity of the new jet
     * @param phi azimuthal angle of the new jet
     )pbdoc",
       py::arg("labels"), py::arg("label_new"),
       py::arg("energy"), py::arg("pt"), py::arg("rapidity"), py::arg("phi")
        )
    .def("DoAllMerges", &Cluster::DoAllMerges,
        R"pbdoc(
     * @brief Automatically complete the merging.
     *)pbdoc")
    .def("IsFinished", &Cluster::IsFinished,
        R"pbdoc(
     * @brief Check if all jets are complete.
     * @return true if all jets are complete
     *)pbdoc")
    .def("GetJets", &Cluster::GetJets,
        R"pbdoc(
     * @brief Get a list of currently existing jets.
     *
     * Mid cluster will return the pseudojets.
     * Order matches order of GetJetConstituents, but
     * otherwise is arbitrary.
     *
     * @return each item in the vector is the label or one completed jet.
     *)pbdoc",
        py::return_value_policy::reference_internal
          )
    .def("GetJetConstituents", &Cluster::GetJetConstituents,
        R"pbdoc(
     * @brief Get a list of particles in each jet.
     *
     * Mid cluster will return the pseudojets.
     * Order matches order of GetJets, but otherwise is arbitrary.
     *
     * @return each item in the vector is all the labels of particles in one jet
     *)pbdoc")
    .def("GetNextFreeLabel", &Cluster::GetNextFreeLabel,
        R"pbdoc(
     * @brief Get the next free label.
     * @return next free label
     *)pbdoc")
     ;


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
