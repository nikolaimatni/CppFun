#include <iostream>

short CountBits(unsigned int x) {
  short num_bits = 0;
  while (x) {
    num_bits += x & 1;
    x >>= 1;
  }
  return num_bits;
}

int main() {
  using namespace std;

  unsigned int x = 1234;
  cout << "num bits in " << x << " is " << CountBits(x) << "\n";
  return 0;
}
