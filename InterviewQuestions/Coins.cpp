#include <iostream>

using namespace std;

int PenniesAndNickels(int total) {
  // assumption is total >=10;
  int ways = 0;
  for (int i = 1; i * 5 <= total; i++) {
    int remainder = total - 5 * i;
    ways++;
  }
  cout << "Pennies and nickels " << ways << "\n";
  return ways;
}

int PenniesNickelsAndDimes(int total) {
  // assumption is total >= 10;
  int ways = 0;
  // loop through possibilities using different number of dimes
  for (int i = 1; i * 10 <= total; i++) {
    int remainder = total - 10 * i;
    cout << "dimes remainder " << remainder << "\n";
    if (remainder < 5)
      ways++;
    else if (remainder < 10)
      ways += 2;
    else
      ways += PenniesAndNickels(remainder);
  }

  cout << "And dimes " << ways << "\n";
  return ways;
}

int PenniesNickelsDimesAndQuarters(int total) {
  // assupmtion is total >= 25;
  int ways = 0;
  // loop thorugh possibilities using different number of quarters
  for (int i = 1; 25 * i <= total; i++) {
    int remainder = total - 25 * i;
    if (remainder < 5)
      ways++;
    else if (remainder < 10)
      ways += 2;
    else
      ways += PenniesNickelsAndDimes(remainder);
  }
  cout << "Quarter ways " << ways << "\n";
  return ways;
}

int HowManyWays(int total) {
  // lets knock out some base cases
  // if total < 5, ways = 1;
  if (total < 5)
    return 1;
  // if total < 10, ways = 2; (all pennies, or 1 nickel + pennies)
  if (total < 10)
    return 2;

  // general approach will be to first find the different number of ways using
  // only pennies and nickels.  Then augment that to incorporate dimes, then
  // augment that to incorporate quarters.
  int ways = 1; // pennies only
  ways += PenniesAndNickels(total);
  ways += PenniesNickelsAndDimes(total);
  if (total >= 25)
    ways += PenniesNickelsDimesAndQuarters(total);

  return ways;
}

int main() {
  //
  int total = 26;

  // all pennies
  // one nickel, 6 pennies
  // two nickels, 1 penny
  // 1 dime, 1 penny
  cout << HowManyWays(total) << "\n";
  return 0;
}
