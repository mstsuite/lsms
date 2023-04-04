//
// Created by F.Moitzi on 20.12.2021.
//

#include <iostream>

#define DEBUG

#if defined(DEBUG)
#define DEBUG_PRINT(args) std::cout << args << std : endl;
#else
#define DEBUG_PRINT(args)
#endif
