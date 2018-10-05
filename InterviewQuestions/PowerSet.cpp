#include <iostream>
#include <vector>

using namespace std;

using Subsets = vector<vector<int>>;

Subsets ComputePowerSet(const vector<int> my_set) {
  if (my_set.size() == 1)
    return Subsets{{my_set[0]}};
  // assumption |my_set| >=2;
  Subsets powerset{};
  // start by computing power of first two positions of my_set;
  int a = my_set[0];
  int b = my_set[1];

  powerset.push_back({});
  powerset.push_back({a});
  powerset.push_back({b});
  powerset.push_back({a, b});

  //  int c = my_set[2];
  // now let's build the powerset for {my_set[0],my_set[1],my_set[2]} from the
  // power set of {my_set[0],my_set[1]}
  // need to append singleton {c}
  // powerset.push_back({c});

  // we can construct the new powerset by appending the next element to all
  // previous subsets and adding those (this is why we need to include the empty
  // set, so we dont' miss the singleton!)
  for (int n = 2; n < my_set.size(); n++) {
    int new_element = my_set[n];
    int current_size = powerset.size(); // how many subsets we need to append to
    for (int i = 0; i < current_size; i++) {
      vector<int> subset = powerset[i]; // get the existing subset
      subset.push_back(new_element);    // append the new element
      powerset.push_back(subset);       // add the new subset to the powerset
    }
  }

  return powerset;
}

void print(const Subsets &powerset) {
  for (int i = 0; i < powerset.size(); i++) {
    cout << "{";
    for (int j = 0; j < powerset[i].size(); j++) {
      cout << powerset[i][j] << ", ";
    }
    cout << "},";
  }
  cout << "\n";
}

int main() {
  vector<int> my_set{1, 2, 3, 4, 5, 6, 7};
  // we know the answer is {1}, {2}, {3}, {1,2}, {1,3}, {2,3}, {1,2,3}
  Subsets powerset = ComputePowerSet(my_set);
  print(powerset);

  return 0;
}
