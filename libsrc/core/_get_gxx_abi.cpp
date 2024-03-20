#include <iostream>

int main() {
  if (__GXX_ABI_VERSION >= 2000 || __GXX_ABI_VERSION < 1000) return 1;
  std::cout << (__GXX_ABI_VERSION % 100);
  return 0;
}
