#ifndef CLUSTER_HXX
#define CLUSTER_HXX

#include <utility>
#include <vector>

class Cluster
{
  public:
    /**
     * @brief Constructor
     *
     * @param labels labels of the particles to be clustered
     * @param energies energies of the particles to be clustered
     * @param pts transverse momenta of the particles to be clustered
     * @param rapidites rapidities of the particles to be clustered
     * @param phis azimuthal angles of the particles to be clustered
     */
    Cluster(const std::vector<int>& labels,
            const std::vector<double>& energies,
            const std::vector<double>& pts,
            const std::vector<double>& rapidites,
            const std::vector<double>& phis);

    /**
     * @brief Say which jets should be merged next.
     * Returned as a pair of labels.
     * If both numbers are the same, that jet is finished.
     * @return pair of labels
     */
    std::pair<int, int> GetNextMerge();
    /**
     * @brief Merge two jets.
     * Calculate the kinematics of the new jet internally, and
     * give it the next free label.
     *
     * @param label_1 label of the first jet
     * @param label_2 label of the second jet
     */
    void DoMerge(int label_1, int label_2);
    /**
     * @brief Merge two jets.
     * Merge the two jets, and assign the combined jet externally
     * provided label and kienmatics.
     *
     * @param label_1 label of the first jet
     * @param label_2 label of the second jet
     * @param energy energy of the new jet
     * @param pt transverse momentum of the new jet
     * @param rapidity rapidity of the new jet
     * @param phi azimuthal angle of the new jet
     */
    void DoMerge(int label_1, int label_2, int label_new,
                 double energy, double pt, double rapidity, double phi);

  private:
    /**
     * @brief A numeric label > 0 for each jet.
     **/
    std::vector<int> m_labels;
    /**
     * @brief The energy of each jet.
     **/
    std::vector<double> m_energies;
    /**
     * @brief The pt of each jet.
     **/
    std::vector<double> m_pts;
    /**
     * @brief The rapidity of each jet.
     **/
    std::vector<double> m_rapidites;
    /**
     * @brief The phi of each jet.
     **/
    std::vector<double> m_phis;

    /**
     * @brief Alternative storage of the kinematics of each jet.
     **/
    std::vector<double> m_pxs;
    std::vector<double> m_pys;
    std::vector<double> m_pzs;

    /**
     * @brief Indicates if a jet is avaliable for further joining.
     **/
    std::vector<bool> m_avaliable;
    /**
     * @brief Indicates if a jet has been merged with the beam.
     **/
    std::vector<bool> m_finished;
    /**
     * @brief The label of the parent jet for each jet.
     * Jets with no parents get a label of -1.
     **/
    std::vector<int> m_parentLabels;

    /**
     * @brief Calculate the kinematics of merging two jets.
     * @param label_1 label of the first jet
     * @param label_2 label of the second jet
     * @return vector of kinematics, energy, pt, rapidity, phi
     **/
    std::vector<double> GetMergedKinematics(int label_1, int label_2);

    /**
     * @brief Get the next free label.
     * @return next free label
     **/
    int GetNextFreeLabel();
};

#endif // CLUSTER_HXX
