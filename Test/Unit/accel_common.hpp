//
// Created by F.Moitzi on 28.09.2022.
//

#ifndef LSMS_ACCEL_COMMON_HPP
#define LSMS_ACCEL_COMMON_HPP

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || \
    defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)

#include "Accelerator/DeviceStorage.hpp"

DeviceStorage *deviceStorage;
DeviceConstants deviceConstants;

#endif

#endif  // LSMS_ACCEL_COMMON_HPP
