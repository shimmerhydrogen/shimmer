#ifndef __COMMON_UTILITIES_H
#define __COMMON_UTILITIES_H

namespace GERG
{
  /// \brief Tells the compiler the parameter is unused
  template<class T>
  static void Unused(const T&) { }
}

#endif // __COMMON_UTILITIES_H
