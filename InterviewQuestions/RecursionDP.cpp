#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

using namespace std;

int fibonacci_dumb(int i) {
  if (i == 0 || i == 1)
    return i;
  return fibonacci_dumb(i - 1) + fibonacci_dumb(i - 2);
}

int fibonacci_dp(int n) {
  if (n == 0 || n == 1)
    return n;

  int im2 = 0;
  int im1 = 1;
  int i = im1 + im2;
  int j = 2;
  while (j < n) {
    im2 = im1;
    im1 = i;
    i = im1 + im2;

    j++;
  }
  return i;
}

int fib_memo(int i, int *memo) {
  if (i == 0 || i == 1)
    return i;

  if (memo[i] == 0)
    memo[i] = fib_memo(i - 1, memo) + fib_memo(i - 2, memo);

  return memo[i];
}

int fibonacci_memo(int n) {
  int *memo = new int[n + 1];
  int fib_n = fib_memo(n, memo);
  delete[] memo;
  return fib_n;
}

int steps_memo(int i, int *memo) {
  if (i <= 1)
    return i;
  if (i == 2)
    return 3;
  if (i == 3)
    return 4;

  if (memo[i] == 0)
    memo[i] = steps_memo(i - 1, memo) + steps_memo(i - 2, memo) +
              steps_memo(i - 3, memo);

  return memo[i];
}

int steps_memo(int n) {
  int *memo = new int[n + 1];

  int ways = steps_memo(n, memo);
  delete[] memo;

  return ways;
}

struct Point {
  int r_, c_;
  Point(int r, int c) : r_{r}, c_{c} {};
  Point up() { return Point(r_, c_ - 1); }
  Point left() { return Point(r_ - 1, c_); }
  bool operator==(const Point &p) { return r_ == p.r_ && c_ == p.c_; }
};

template <size_t r, size_t c> using Map = array<array<int, c>, r>;
template <size_t r, size_t c>
bool GetPath(Map<r, c> world, Point start, vector<Point> &path,
             vector<Point> &failed) {
  // check for failure of this path direction
  if (start.r_ < 0 || start.c_ < 0 || !world[start.r_][start.c_] ||
      (find(failed.begin(), failed.end(), start)) != failed.end()) {
    return false;
  }

  // check if we've reached the top left corner yet
  if (start.r_ == 0 && start.c_ == 0) {
    // we've found a path! so let's add this to our list and recurse up
    path.push_back(start);
    return true;
  }

  if (GetPath(world, start.up(), path, failed) ||
      GetPath(world, start.left(), path, failed)) {
    // if going up or left gets a path, push back this point and recurse up
    path.push_back(start);
    return true;
  }

  // if not this is a failure point, so add it to failure list and return false
  failed.push_back(start);
  return false;
}

int main() {
  cout << "Fib_memo(20) " << fibonacci_memo(20) << "\n";
  cout << "Fib_dp(20) " << fibonacci_dp(20) << "\n";

  cout << "Steps(4) " << steps_memo(4) << "\n";

  vector<vector<int>> world = {
      {1, 1, 1, 1}, {1, 1, 1, 1}, {0, 1, 0, 1}, {1, 1, 0, 1}, {1, 1, 1, 1}};

  Map<5, 4> staticworld;
  staticworld[0] = {1, 1, 1, 1};
  staticworld[1] = {1, 1, 1, 1};
  staticworld[2] = {0, 1, 0, 1};
  staticworld[3] = {1, 1, 0, 1};
  staticworld[4] = {1, 1, 1, 1};

  vector<Point> path{};
  vector<Point> failed{};

  if (GetPath(staticworld, Point(5 - 1, 4 - 1), path, failed)) {
    for (int i = 0; i < path.size(); i++)
      cout << "(" << path[i].r_ << "," << path[i].c_ << "), ";
    cout << "\n";
  } else
    cout << "no path exists!\n";

  return 0;
}
