//
// Created by F.Moitzi on 21.12.2022.
//

#ifndef LSMS_UTILS_H
#define LSMS_UTILS_H

#include <type_traits>

namespace lsms {

template<typename E>
constexpr auto convert(E e)
{
  return static_cast<typename std::underlying_type<E>::type>(e);
}

}

#endif // LSMS_UTILS_H
