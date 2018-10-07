#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

template <typename T, size_t rows, size_t cols>
using StaticMatrix = array<array<T, cols>, rows>;

template <typename T> using DynamicMatrix = vector<vector<T>>;

template <size_t cols> void print_array(int(array[][cols])) {
  for (int i = 0; i < 5; i++)
    cout << array[i][0] << "\n";
}

template <typename T> void print(T array) {
  for (int i = 0; i < 5; i++)
    cout << array[i][0] << "\n";
}

class TestClass {
public:
  vector<int> vec_;
  int *array_;
  int n_;

  TestClass(vector<int> vec, int n) : vec_(vec), array_(new int[n]), n_(n) {}

  TestClass(const TestClass &other) = delete;
  TestClass &operator=(const TestClass &other) = delete;

  TestClass(TestClass &&other)
      : vec_(move(other.vec_)), array_(other.array_), n_(other.n_) {
    other.array_ = nullptr;
  }

  TestClass &operator=(TestClass &&other) {
    if (this == &other)
      return *this;

    vec_ = move(other.vec_);
    swap(array_, other.array_);
    n_ = other.n_;
  }
};

int main() {

  int arr[5][1]{{1}, {2}, {3}, {4}, {5}};
  print_array(arr);

  StaticMatrix<int, 5, 1> matrix;
  for (int row = 0; row < 5; row++)
    matrix[row][0] = row + 1;

  print(arr);
  print(matrix);

  auto pr = [](array<int, 5> arr) {
    for (int i = 0; i < 5; i++)
      cout << arr[i] << " ";
    cout << "\n";
  };

  array<int, 5> myarr{1, 2, 3, 4, 5};

  pr(myarr);

  for (auto r = myarr.crbegin(); r != myarr.crend(); r++)
    cout << *r << " ";
  cout << "\n";

  vector<int> vec(5, 10);

  for (auto it = vec.begin(); it != vec.end(); it++)
    cout << *it << " ";
  cout << "\n";

  for (int i = 0; i < 5; ++i)
    cout << i << " ";
  cout << "\n";

  int i = 10;
  while (i--)
    cout << i << " ";
  cout << "\n";

  vector<int> myvec{2, -10, 11, 5, 55, 8};
  int sorted = 0;
  while (sorted < myvec.size()) {
    int min_pos = sorted;
    for (int i = sorted + 1; i < myvec.size(); i++)
      min_pos = myvec[i] < myvec[min_pos] ? i : min_pos;
    swap(myvec[sorted], myvec[min_pos]);
    sorted++;
  }

  cout << myvec.size() << "\n";

  for (auto it = myvec.begin(); it != myvec.end(); it++)
    cout << *it << " ";
  cout << "\n";
  for_each(myvec.begin(), myvec.end(), [](int val) { cout << val << " "; });
  cout << "\n";

  srand(time(NULL));
  cout << double(rand()) / RAND_MAX << "\n";

  auto prt = [](const vector<int> &vec) {
    for_each(vec.begin(), vec.end(), [](int val) { cout << " " << val; });
    cout << "\n";
  };

  array<int, 5> arr5{1, 2, 3, 4, 5};
  array<int, 5> arr5cpy(move(arr5));

  pr(arr5);
  pr(arr5cpy);

  vector<int> newvec(move(myvec));
  prt(newvec);
  prt(myvec);

  int m = 5;
  unique_ptr<int[]> p(new int[m]{1, 2, 3, 4, 5});
  // auto pp = make_unique<int[]>(m);
  int *ppp = new int[m]{1, 2, 3, 4, 5};
  unique_ptr<int[]> pppp(ppp);

  TestClass test(newvec, 10);
  delete[] test.array_;
  if (test.array_)
    cout << "test not empty after delete\n";

  TestClass newtest(move(test));
  prt(newtest.vec_);
  prt(test.vec_);
  if (!newtest.array_)
    cout << "nullptr\n";
  return 0;
}
