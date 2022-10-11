
#ifndef MUST_IO_HPP
#define MUST_IO_HPP

#include <fstream>
#include <iostream>
#include <string>

template <class T1, class T2>
void write_to_file(T1 R, T2 value, const std::string &filename) {
  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < R.size(); i++) {
    file << std::setw(24) << std::setprecision(16) << std::scientific << R[i]
         << " ";
    file << std::setw(24) << std::setprecision(16) << std::scientific
         << value[i] << std::endl;
  }

  file.close();
}

#endif  // MUST_IO_HPP
