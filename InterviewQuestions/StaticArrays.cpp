#include <array>
#include <iostream>
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

int main() {
  //
  int array[5][1]{{1}, {2}, {3}, {4}, {5}};
  print_array(array);

  StaticMatrix<int, 5, 1> matrix;
  for (int row = 0; row < 5; row++)
    matrix[row][0] = row + 1;

  print(array);
  print(matrix);

  return 0;
}
