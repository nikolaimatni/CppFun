#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

struct Box {
  double h_, w_, d_; // width, height, depth
  vector<int>
      stack_; // index of boxes allowed to be stacked directly on this box
  Box(double h, double w, double d) : h_(h), w_(w), d_(d), stack_{} {}
};

using Boxes = vector<Box>;

double MaxCostToGo(vector<int> stack, vector<double> V) {
  double ctg = 0;
  for (int i : stack) {
    cout << V[i] << " " << ctg << "\n";
    ctg = max(ctg, V[i]);
  }
  return ctg;
}

double DPSolution(vector<double> &Vm1, vector<double> &V, const Boxes &boxes) {
  for (int layer = boxes.size() - 1; layer >= 0; layer--) {
    for (int i = 0; i < boxes.size(); i++) {
      Vm1[i] = boxes[i].h_ +
               MaxCostToGo(boxes[i].stack_,
                           V); // find cost to go for box i in current layer
    }
    V = Vm1; // update cost to go and go down a layer
  }
  return (*max_element(V.begin(), V.end()));
}

void InitializeStack(Box &box, Boxes &boxes, int j) {
  for (int i = 0; i < boxes.size(); i++) {
    Box &next = boxes[i];
    if (i != j && box.h_ >= next.h_ && box.w_ >= next.w_ && box.d_ >= next.d_) {
      box.stack_.push_back(i);
      cout << "adding box " << i << " to  stack of box " << j << "\n";
    }
  }
}

double MaxHeight(Boxes &boxes) {
  // Use Dynamic Progamming
  // V_{l-1}(x) = x.h_ + max_{z \in x.stack_} V_{l}(z); l = layer
  // V_{top}(x) = x.h_;
  int n = boxes.size();
  vector<double> Vm1(n, 0); // use to store V_{l-1}
  vector<double> V(n, 0);   // use to store V_{l}
  for (int i = 0; i < n; i++) {
    // initiliaze terminal value function, and populates stack_ for each box
    V[i] = boxes[i].h_;
    InitializeStack(boxes[i], boxes, i);
  }
  return DPSolution(Vm1, V, boxes);
}

int main() {
  Boxes boxes;

  boxes.push_back(Box(2, 5, 8));
  boxes.push_back(Box(0.1, 0.1, 0.1));
  boxes.push_back(Box(5, 2, 8));
  boxes.push_back(Box(8, 8, 8));
  boxes.push_back(Box(2, 2, 2));

  cout << "The maximum height for our stack is " << MaxHeight(boxes) << ".\n";

  return 0;
}
