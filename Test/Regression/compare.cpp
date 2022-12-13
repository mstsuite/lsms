//
// Created by F.Moitzi on 30.11.2022.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

double REL_TOL = 1.0e-7;
double ABS_TOL = 1.0e-12;

bool isclose(double ref, double val) {
  return std::fabs(ref - val) <= (ABS_TOL + REL_TOL * std::fabs(val));
}


int main(int argc, char **argv) {

  if (argc < 3) { /* validate one argument given for filename */
    std::cerr << "error: insufficient arguments.\n"
                 "usage: " << argv[0] << " filename.\n";
    return 1;
  }


  int counter = 0;
  double total_energy = 0.0;
  double fermi_level = 0.0;
  double mag_mom = 0.0;
  double rms = 0.0;

  std::ifstream fin(argv[1]);

  if (!fin.is_open()) {
    return 1;
  }

  std::ifstream fin2(argv[2]);

  if (!fin2.is_open()) {
    return 1;
  }

  std::vector<double> tot1;
  std::vector<double> fermi1;
  std::vector<double> mag_mom1;
  std::vector<double> rms1;

  std::vector<double> tot2;
  std::vector<double> fermi2;
  std::vector<double> mag_mom2;
  std::vector<double> rms2;

  while (fin >> counter >> total_energy >> fermi_level >> mag_mom >> rms) {

    tot1.emplace_back(total_energy);
    fermi1.emplace_back(fermi_level);
    mag_mom1.emplace_back(mag_mom);
    rms1.emplace_back(rms);

  }

  while (fin2 >> counter >> total_energy >> fermi_level >> mag_mom >> rms) {

    tot2.emplace_back(total_energy);
    fermi2.emplace_back(fermi_level);
    mag_mom2.emplace_back(mag_mom);
    rms2.emplace_back(rms);

  }



  std::vector<bool> result(tot1.size());

  std::transform(tot1.begin(), tot1.end(), tot2.begin(), result.begin(), isclose);
  if (!std::all_of(result.begin(), result.end(), [](bool v) { return v; })) {
    return EXIT_FAILURE;
  }

  std::transform(fermi1.begin(), fermi1.end(), fermi2.begin(), result.begin(), isclose);
  if (!std::all_of(result.begin(), result.end(), [](bool v) { return v; })) {
    return EXIT_FAILURE;
  }
  std::transform(mag_mom1.begin(), mag_mom1.end(), mag_mom2.begin(), result.begin(), isclose);
  if (!std::all_of(result.begin(), result.end(), [](bool v) { return v; })) {
    return EXIT_FAILURE;
  }

  std::transform(rms1.begin(), rms1.end(), rms2.begin(), result.begin(), isclose);
  if (!std::all_of(result.begin(), result.end(), [](bool v) { return v; })) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

