
#include "cvxcanon/util/Utils.hpp"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include <cstdlib>
#include <memory>
#include <string>

std::string string_printf(const std::string fmt_str, ...) {
  int final_n, n = fmt_str.size() * 2;
  std::string str;
  std::unique_ptr<char[]> formatted;
  va_list ap;

  while (1) {
    formatted.reset(new char[n]);
    strcpy(&formatted[0], fmt_str.c_str());  // NOLINT(runtime/printf)
    va_start(ap, fmt_str);
    final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
    va_end(ap);
    if (final_n < 0 || final_n >= n)
      n += abs(final_n - n + 1);
    else
      break;
  }

  return std::string(formatted.get());
}
