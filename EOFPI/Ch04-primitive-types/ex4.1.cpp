#include <cmath>
#include <iostream>
// we want to compute the parity of a *large number* of 64 bit words

short Parity(unsigned long x) {
  // brute force
  short parity = 0;
  while (x) {
    parity ^= (x & 1);
    x >>= 1;
  }
  return parity;
}

short Parity2(unsigned long x) {
  // erase lowest set bit over and over, XOR'ing 1 to parity each time = parity
  // = sum(1's mod 2)
  short parity = 0;
  while (x) {
    parity ^= 1;
    x &= (x - 1); // erase lowest set bit
  }
  return parity;
}

short Parity3(unsigned long x) {
  // XORs are associative and commutative, so can group things however we want
  x ^= x >> 32; // b63-b32 with b31-b0... now only b31-b0 are relevant
  x ^= x >> 16; // b31-b16 with b15-b0... now only b15-0 are relevant...
  x ^= x >> 8;
  x ^= x >> 4;
  x ^= x >> 2;
  x ^= x >> 1; // now only the last bit is relevant, and is the parity!
  return x & 0b1;
}

int Reverse(int x) {

  int result = 0;
  int x_remaining = std::abs(x);
  while (x_remaining) {
    result = result * 10 + x_remaining % 10;
    x_remaining /= 10;
  }

  return x < 0 ? -result : result;
}

double Power(double x, int y) {
  double result = 1.0;
  long long power = y;
  if (y < 0) {
    power = -power;
    x = 1.0 / x;
  }
  while (power) {
    if (power & 1) { // if power is odd
      result *= x;
      std::cout << "updating result by *=x " << x << "\n";
    }

    x *= x;
    power >>= 1;
    std::cout << x << " , " << power << "\n";
  }
  return result;
}

int main() {
  using namespace std;

  unsigned long x = 1234568;

  if (Parity(x) == Parity3(x))
    cout << "The parity of the word is: " << Parity(x) << " \n";

  cout << "The reverse of -146 is " << Reverse(-146) << "\n";

  cout << "2^5 = " << Power(2.0, 5) << "\n";

  return 0;
}
