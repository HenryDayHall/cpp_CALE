#include "cluster.hxx"
#include <iostream>
#include <vector>

int main(){
  Cluster my_algo = Cluster(0.1, 0., 15);

  std::vector<int> labels = {0, 1, 2, 3, 4};
  std::vector<double> energies = {10, 10, 10, 10, 10};
  std::vector<double> pts = {10, 10, 10, 10, 10};
  std::vector<double> rapidites = {0.1, 0.1, 0.3, 0.4, 0.4};
  std::vector<double> phis = {0.1, 0.1, 0.3, 0.4, 0.4};

  my_algo.SetInputs(labels, energies, pts, rapidites, phis);

  my_algo.DoAllMerges();

  std::vector<int> completed_jets = my_algo.GetJets();
  std::vector<std::vector<int>> constituents = my_algo.GetJetConstituents();

  int n_jets = completed_jets.size();
  std::cout << "There are " << n_jets << " jets." << std::endl;
  const std::vector<double> jet_energies = my_algo.GetEnergies(completed_jets);
  const std::vector<double> jet_pts = my_algo.GetPts(completed_jets);
  const std::vector<double> jet_rapidites = my_algo.GetRapidites(completed_jets);
  const std::vector<double> jet_phis = my_algo.GetPhis(completed_jets);
  for (int jet_n = 0; jet_n < n_jets; jet_n++){
    std::cout << "Jet " << jet_n << " (" << jet_energies[jet_n] << ", " << jet_pts[jet_n] << ", " << jet_rapidites[jet_n] << ", " << jet_phis[jet_n] << ") has constituents: " << std::endl;
    const std::vector<double> constitutent_energies = my_algo.GetEnergies(constituents[jet_n]);
    const std::vector<double> constitutent_pts = my_algo.GetPts(constituents[jet_n]);
    const std::vector<double> constitutent_rapidites = my_algo.GetRapidites(constituents[jet_n]);
    const std::vector<double> constitutent_phis = my_algo.GetPhis(constituents[jet_n]);
    for (int i = 0; i < constituents[jet_n].size(); i++){
      std::cout << "  (" << constitutent_energies[i] << ", " << constitutent_pts[i] << ", " << constitutent_rapidites[i] << ", " << constitutent_phis[i] << ")" << std::endl;
    }
    std::cout << std::endl;
  }
  
}

