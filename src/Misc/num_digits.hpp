//
// Created by F.Moitzi on 18.10.2022.
//

#ifndef LSMS_NUM_DIGITS_HPP
#define LSMS_NUM_DIGITS_HPP

template <class T>
static int num_digits(T number) {
  int digits = 0;

  if (number < 0) digits = 1;  // remove this line if '-' counts as a digit

  while (number) {
    number /= 10;
    digits++;
  }
  return digits;
}

#endif  // LSMS_NUM_DIGITS_HPP
