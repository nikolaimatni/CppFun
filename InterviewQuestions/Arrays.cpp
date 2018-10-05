#include <iostream>
#include <memory> // for std::unique_ptr

using namespace std;

// homebrew implementation of vector<T> using unique_ptr
template <typename T> class SmartArrayList {
private:
  unique_ptr<T> data_;
  int size_;
  int capacity_;

  void Expand() {

    // double capacity and copy over
    unique_ptr<T> new_data(new T[2 * capacity_ + 1]);

    capacity_ = 2 * capacity_ + 1;

    // copy over existing data
    for (int i = 0; i < size_; ++i)
      new_data.get()[i] = data_.get()[i];

    // reassign data_ to point to new_data
    data_ = move(new_data);
  }

public:
  SmartArrayList() : data_(unique_ptr<T>(nullptr)), size_(0), capacity_(0) {}
  SmartArrayList(size_t size)
      : data_(new T[size]), size_(size), capacity_(size) {}
  size_t size() { return size_; }
  void push_back(T item) {
    if (size_ == capacity_)
      Expand();

    data_.get()[size_] = item;
    size_++;
  }

  void print() {
    for (int i = 0; i < size_; ++i)
      cout << data_.get()[i] << " ";
    cout << "\n";
  }

  T *data() { return data_.get(); }

  T &operator[](const int i) { return data_.get()[i]; }

  ~SmartArrayList() {}
};

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

  // copy semantics
  Matrix(const Matrix &m) {
    cout << "copy constructor!\n";
    // perform a deep copy
    N_ = m.N_;
    matrix_ = new T *[N_];
    for (int i = 0; i < N_; i++)
      matrix_[i] = new T[N_];

    for (int i = 0; i < N_; i++)
      for (int j = 0; j < N_; j++)
        matrix_[i][j] = m.matrix_[i][j];
  }

  Matrix &operator=(const Matrix &m) {
    cout << "copy assignment!\n";
    // check for self assignment;
    if (&m == this)
      return *this;

    N_ = m.N_;
    // now perform a deep copy
    if (!matrix_) {
      // if memory hasn't been allocated for the matrix yet, allocate it
      matrix_ = new T *[N_];
      for (int i = 0; i < N_; i++)
        matrix_[i] = new T[N_];
    }

    // now copy data
    for (int i = 0; i < N_; i++)
      for (int j = 0; j < N_; j++)
        matrix_[i][j] = m.matrix_[i][j];

    return *this;
  }

  // move semantics
  Matrix(Matrix &&m) {
    cout << "move constructor!\n";
    // initialize local variable;
    N_ = m.N_;

    // claim ownership of source's matrix_
    matrix_ = m.matrix_;

    // set source Matrix to ``default'' state
    m.matrix_ = nullptr;
    m.N_ = 0;
  }

  Matrix &operator=(Matrix &&m) {
    // check for self-assignment
    cout << "move assignment!\n";
    if (&m == this)
      return *this;

    // initilize size
    N_ = m.N_;

    // claim ownership of source's matrix_
    matrix_ = m.matrix_;

    // set source matrix to default state
    m.matrix_ = nullptr;
    m.N_ = 0;
  }

  size_t size() { return N_; }

  T *row(int i) { return matrix_[i]; }

  T val(int i, int j) const { return matrix_[i][j]; }

  T **data() { return matrix_; }

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
  SmartArrayList<int> my_array;
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
  Matrix<int> rotated = Rotate((original), N);

  rotated.print();
  print_line();
  original.Rotate();
  original.print();

  Matrix<int> copy(move(rotated));
  Matrix<int> empty_vessel(N);
  empty_vessel = move(copy);
  print_line();
  empty_vessel.print();

  if (rotated.size()) {
    print_line();
    cout << "printing rotated\n";
    rotated.print();
  }

  if (copy.size()) {
    print_line();
    cout << "printing copy \n";
    copy.print();
  }

  return 0;
}
