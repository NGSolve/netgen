#include <iostream>

int main() {
  #ifdef _GLIBCXX_USE_CXX11_ABI
  if(_GLIBCXX_USE_CXX11_ABI)
    std::cout << 1;
  else
    std::cout << 0;
  #else // _GLIBCXX_USE_CXX11_ABI
    std::cout << 0;
  #endif // _GLIBCXX_USE_CXX11_ABI
  return 0;
}
