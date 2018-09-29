#include <iostream>

using namespace std;
// homebrew implementation of vector<T>
template <typename T> class ArrayList {
private:
  T *data_;
  int size_;
  int capacity_;

  void Expand() {

    // double capacity and copy over
    T *new_data = new T[2 * capacity_ + 1];
    capacity_ = 2 * capacity_ + 1;

    // copy over existing data
    for (int i = 0; i < size_; ++i)
      new_data[i] = data_[i];

    // clean up old data
    delete[] data_;

    // reassign data_ to point to new_data
    data_ = new_data;
  }

public:
  ArrayList() : data_(nullptr), size_(0), capacity_(0) {}
  ArrayList(size_t size) : data_(new T[size]), size_(size), capacity_(size) {}
  size_t size() { return size_; }
  void push_back(T item) {
    if (size_ == capacity_)
      Expand();

    data_[size_] = item;
    size_++;
  }

  void print() {
    for (int i = 0; i < size_; ++i)
      cout << data_[i] << " ";
    cout << "\n";
  }

  T *data() { return data_; }

  T &operator[](const int i) { return data_[i]; }

  ~ArrayList() { delete[] data_; }
};

template <typename T> bool IsUnique1(T *array, size_t size) {
  for (int idx = 0; idx < size; idx++) {
    for (int i = idx + 1; i < size; i++) {
      if (array[i] == array[idx])
        return false;
    }
  }
  return true;
}

template <typename T> class Matrix {
private:
  T **matrix_;
  size_t N_;

public:
  Matrix(const size_t N) : matrix_(new T *[N]), N_(N) {
    for (int i = 0; i < N_; i++)
      matrix_[i] = new T[N];
  }
  Matrix(T **init, size_t N) : matrix_(init), N_(N) {}

  T *row(int i) { return matrix_[i]; }

  T val(int i, int j) const { return matrix_[i][j]; }

  // when we call this on a non-const object, pass back reference to allow for
  // modifications
  T &operator()(int i, int j) { return matrix_[i][j]; }

  // when we call this on a cont object, only pass back a copy
  T operator()(int i, int j) const { return matrix_[i][j]; }

  void print() const {
    for (int i = 0; i < N_; i++) {
      cout << "[";
      for (int j = 0; j < N_; j++) {
        cout << matrix_[i][j] << ", ";
      }
      cout << "]\n";
    }
    cout << " --------------------------\n";
  }

  T **matrix() const { return matrix_; }

  void Rotate() {
    // starting at outside layer, split each layer of the matrix into a top,
    // right, bottom, and left vector of length layer_size -1, then swap them ou
    // element by element, then move inwards.
    int n = N_;
    int layer = n;
    for (int off = 0; off < n / 2; off++) {
      layer--;
      // cout << "layer subvec length is " << layer << "\n";
      for (int i = 0; i < layer; i++) {
        // temp = top[i];
        T temp = matrix_[off][off + i];
        //  cout << "Top[i] pre swap is " << temp << "\n";

        // top[i] = left[i];
        matrix_[off][off + i] = matrix_[n - 1 - off - i][off];
        //  cout << "Left[i] pre swap is " << matrix_[n - 1 - off - i][off] <<
        //  "\n";

        // left[i] = bottom[i]
        matrix_[n - 1 - off - i][off] = matrix_[n - 1 - off][n - 1 - off - i];
        // cout << "Bottom[i] pre swap is "
        //   << matrix_[n - 1 - off][n - 1 - off - i] << "\n";

        // bottom[i] = right[i]
        matrix_[n - 1 - off][n - 1 - off - i] = matrix_[off + i][n - 1 - off];
        // cout << "Right[i] pre swap is " << matrix_[off + i][n - 1 - off]
        //    << "\n";

        // right[i] = temp
        matrix_[off + i][n - 1 - off] = temp;
      }
    }
  }

  ~Matrix() {
    for (int i = 0; i < N_; i++)
      delete[] matrix_[i];

    delete matrix_;
  }
};

template <typename T> Matrix<T> Rotate(const Matrix<T> &input, const size_t N) {
  Matrix<T> rotated(N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      rotated(j, N - 1 - i) = input(i, j); // val(i, j);
    }
  }

  return rotated;
}

int main() {
  ArrayList<int> my_array;
  for (int i = 0; i < 10; i++) {
    my_array.push_back(i);
  }
  cout << "\n";

  my_array.print();

  my_array[5] = 324;

  my_array.print();

  my_array[5] = 3;
  cout << IsUnique1(my_array.data(), my_array.size()) << "\n";

  int N = 4;
  auto print_line = []() { cout << "------------------\n"; };

  print_line();

  int **matrix = new int *[N];
  for (int i = 0; i < N; i++) {
    matrix[i] = new int[N];
    cout << "[";
    for (int j = 0; j < N; j++) {
      matrix[i][j] = i + 2 * j * j;
      cout << matrix[i][j] << ", ";
    }
    cout << "]\n";
  }

  print_line();

  Matrix<int> original(matrix, N);
  // original.print();

  original(0, 0) = 100;

  int &q = original(2, 2);
  q = -111;
  original.print();

  print_line();
  Matrix<int> rotated = Rotate(original, N);

  rotated.print();
  print_line();
  original.Rotate();
  original.print();

  return 0;
}
