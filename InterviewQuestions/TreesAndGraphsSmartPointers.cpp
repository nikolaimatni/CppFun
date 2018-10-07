#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <memory>
#include <queue>
#include <vector>

using namespace std;

template <typename Value> struct Node {
  shared_ptr<Node> left_;
  shared_ptr<Node> right_;
  Value val_;
  bool seen_ = false;
  int dist_ = 0;
  Node(Value val) : val_(val), left_{nullptr}, right_{nullptr} {}

  ~Node() { cout << "destroying node with value " << val_ << "\n"; }
};

using NodePtr = shared_ptr<Node<int>>;

NodePtr BinTree(vector<int>::iterator first, vector<int>::iterator last) {
  auto mid = first + (last - first) / 2;
  int mid_val = *mid;

  NodePtr root(new Node<int>(mid_val));
  if (first < mid) {

    root->left_ = BinTree(first, mid);
  }
  if (mid + 1 < last) {

    root->right_ = BinTree(mid + 1, last);
  }
  return root;
}

void InOrderTraversal(const NodePtr &root) {
  cout << " " << root->val_;
  if (root->left_)
    InOrderTraversal(root->left_);
  if (root->right_)
    InOrderTraversal(root->right_);
}

class Adjacency {
private:
  vector<vector<int>> adj_;

public:
  Adjacency(int n) : adj_(vector<vector<int>>(n, vector<int>{})) {}
  void AddEdge(int src, int dest) {
    adj_[src].push_back(dest);
    adj_[dest].push_back(src);
  }
  vector<int> get_neighbors(int n) { return adj_[n]; }
  void print() {
    for (int i = 0; i < adj_.size(); i++) {
      cout << "Node " << i << " has neighbors:";
      for (int j = 0; j < adj_[i].size(); j++)
        cout << " " << adj_[i][j];
      cout << "\n";
    }
  }
};

template <typename T, typename Compare = less<T>> class SearchTree {
private:
  vector<shared_ptr<Node<T>>> nodes_;
  function<bool(const T &, const T &)> comp_;
  Adjacency adj_;

public:
  shared_ptr<Node<T>> get_node(const T &val) {
    auto node = find_if(
        nodes_.begin(), nodes_.end(),
        [&val](shared_ptr<Node<T>> const &n) { return n->val_ == val; });
    if (node != nodes_.end())
      return *node;

    cout << "Making a new node with val " << val << "\n";
    auto new_node = shared_ptr<Node<T>>(new Node<T>(val));
    nodes_.push_back(new_node);
    return new_node;
  }

  SearchTree(Adjacency adj, Compare comp = Compare())
      : nodes_{}, comp_(comp), adj_(adj) {}

  int BFSDistBetween(const T &src, const T &dest) {
    if (src == dest)
      return 0;

    auto root = get_node(src);

    queue<decltype(root)> q;
    q.push(root);
    root->seen_ = true;
    root->dist_ = 0;

    while (q.size()) {
      auto node = q.front();
      q.pop();

      for (int i : adj_.get_neighbors(node->val_)) {
        auto nbr = get_node(i);
        if (nbr->seen_ == false) {
          nbr->seen_ = true;
          nbr->dist_ = node->dist_ + 1;
          if (nbr->val_ == dest)
            return nbr->dist_;
          q.push(nbr);
        }
      }
    }
    return -1;
  }
};

int main() {
  //
  vector<int> input{1, 2, 3, 4, 5, 18, 21, 28};

  NodePtr bintree = BinTree(input.begin(), input.end());

  InOrderTraversal(bintree);
  cout << "\n";

  Adjacency adj(8);
  adj.AddEdge(0, 1);
  adj.AddEdge(0, 3);
  adj.AddEdge(1, 2);
  adj.AddEdge(3, 4);
  adj.AddEdge(3, 7);
  adj.AddEdge(4, 5);
  adj.AddEdge(4, 6);
  adj.AddEdge(4, 7);
  adj.AddEdge(5, 6);
  adj.AddEdge(6, 7);

  adj.print();

  SearchTree<int> searchtree(adj);
  cout << searchtree.BFSDistBetween(5, 2) << "\n";

  // Implement a BFS search tree for path planning using smart pointers
  return 0;
}
