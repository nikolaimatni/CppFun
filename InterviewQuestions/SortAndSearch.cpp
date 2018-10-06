#include <iostream>
#include <memory>
#include <vector>

using namespace std;

void print(const vector<int> &A) {
  for (auto a : A)
    cout << a << ",";
  cout << "\n";
}

template <typename T>
typename vector<T>::iterator BinaryFind(vector<T> &A, T val) {
  auto first = A.begin();
  auto last = A.end();
  decltype(A.begin()) mid;
  while (first <= last) {
    mid = first + distance(first, last) / 2;
    if (val < *mid) {
      last = mid - 1;
    } else if (val > *mid) {
      first = mid + 1;
    } else
      return mid;
  }
  return A.end();
}

auto BinarySearch(vector<int> &A, int val) -> decltype(A.begin()) {

  auto first = A.begin();
  auto last = A.end();
  auto mid = first + distance(first, last) / 2;

  while (first <= last) {
    if (val < *mid) {
      last = mid - 1;
      mid = first + distance(first, mid) / 2;
    }
    if (val > *mid) {
      first = mid + 1;
      mid = last - distance(mid, last) / 2;
    }
    if (val == *mid)
      return mid;
  }

  return mid;
}

vector<int>::iterator BinarySearch(vector<int>::iterator first,
                                   vector<int>::iterator last, int val) {

  cout << "Searching between values " << *first << " and " << *last << "\n";
  vector<int>::iterator mid; // = first + distance(first, last) / 2;

  while (first <= last) {
    mid = first + distance(first, last + 1) / 2;
    cout << "mid " << *mid << "\n";
    cout << "bad mid " << *(first + distance(first, last) / 2) << "\n";
    if (val < *mid) {
      last = mid - 1;
    } else if (val > *mid) {
      first = mid + 1;
    } else if (val == *mid)
      return mid;
  }

  return mid;
}

void SortedMerge(vector<int> &A, vector<int> &B) {
  // assumption A.size() !=
  // start with brute force approach, then update to implement binary search to
  // find insertion point
  auto idx =
      A.begin(); // this is the last spot that we inserted an element of B in.
                 // Since they are both sorted, we only need to backwards in A

  for (int b = 0; b < B.size(); b++) {
    cout << "B[b] = " << B[b] << "\n";
    bool inserted = false;
    idx = BinarySearch(idx, A.end(), B[b]);
    cout << "BinSearch says insert at " << *idx << "\n";
    if (idx < A.end()) {
      idx = A.insert(idx + 1, B[b]);

      print(A);
    } else {
      cout << "shouldn't get here\n";
      A.insert(A.end(), B.begin() + b, B.end());
      break;
    }
  }
}

int BinaryFind(const vector<int> &arr, int first, int val) {
  int last = arr.size() - 1;
  int mid;
  while (first <= last) {
    mid = (first + last) / 2;
    cout << "first " << first << " and last " << last << "\n";
    if (val < arr[mid]) {
      last = mid - 1;
    } else if (val > arr[mid]) {
      first = mid + 1;
    } else {
      return mid;
    }
  }
  return -1;
}

int BinaryFind(const vector<int> &arr, int first, int last, int val) {
  int mid;
  while (first <= last) {
    mid = (first + last) / 2;
    cout << "first " << first << " and last " << last << "\n";
    if (val < arr[mid]) {
      last = mid - 1;
    } else if (val > arr[mid]) {
      first = mid + 1;
    } else {
      return mid;
    }
  }
  return -1;
}

int RotatedSearch(const vector<int> &arr, int val) {
  int pivot = 0;
  int n = arr.size();
  for (int i = 0; i < arr.size() - 1; i++) {
    if (arr[i] == val)
      return i;
    if (arr[n - 1 - i] == val)
      return n - 1 - i;
    if (arr[i] > arr[i + 1]) {
      pivot = i + 1;
      break;
    }
    if (arr[n - 1 - i] < arr[n - 2 - i]) {
      pivot = n - 1 - i;
      break;
    }
  }

  cout << "pivot value " << arr[pivot] << "\n";
  return BinaryFind(arr, pivot, val);
}

int ListySearch(const vector<int> &A, int val) {
  // find end of list
  int end = 1;
  while (A[end - 1] >= 0) // O(logn)
    end *= 2;

  // now i know the real end is between end and end/2
  // let's binary search on something O(n) => O(logn)
  int first = end / 2;
  int last = end;
  int mid;
  while (first <= last) {
    mid = (first + last) / 2;
    if (A[mid] < 0) { // mid is too far to the right
      last = mid - 1;
    }
    if (A[mid] >= 0) { // mid might be too far to the left
      first = mid + 1;
    }
  }
  // mid should be pointing at the very last element
  cout << "mid " << mid << " points to " << A[mid] << "\n";

  return BinaryFind(A, 0, mid - 1, val);
}

int main() {
  vector<int> A{1, 5, 8, 23};
  vector<int> B{3, 4, 7, 10, 11};

  SortedMerge(A, B);
  print(A);

  vector<int> C{3, 7, 8, 23, 25};
  cout << *BinarySearch(C.begin(), C.end(), 10) << "\n";
  cout << *BinarySearch(C, 10) << "\n";
  cout << *(C.begin() + 1 + distance(C.begin() + 1, C.end()) / 2) << "\n";
  cout << *BinaryFind(C, 7) << "\n";
  // auto insert_point = BinarySearch(A.begin(), A.end(), 6);
  // cout << *insert_point << "\n";
  // A.insert(insert_point + 1, 6);
  // print(A);

  vector<int> D{21, 23, 25, 1, 3, 4, 5, 7, 10, 14, 15, 16, 17, 18, 19, 20};
  cout << "5 is in position " << RotatedSearch(D, 5) << "\n";

  vector<int> E(1000, -1);
  for (int i = 0; i < 10; i++)
    E[i] = 2 * i;

  cout << "ListySearch(20) = " << ListySearch(E, 4) << "\n";
  cout << "BinaryFind(20) = " << BinaryFind(E, 0, 10, 4) << "\n";
  return 0;
}
