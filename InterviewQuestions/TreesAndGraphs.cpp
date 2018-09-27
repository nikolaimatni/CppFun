#include <array>
#include <iostream>
#include <list>

using namespace std;

struct Node {
  int val_;
  Node *left_;
  Node *right_;

  Node(int val) : val_(val) {}
};

template <size_t N>
Node *BinTree(const array<int, N> &input, int start, int end) {
  int mid_pos = start + (end - start) / 2;
  int mid_val = input[mid_pos];

  cout << "Making new node with value: " << mid_val << ".\n";

  Node *node = new Node(mid_val);

  // parition array into left component [start, mid_pos) and right component
  // (mid_pos,end)
  int end_left = mid_pos;
  int start_right = mid_pos + 1;

  // check if we are done going to the left, if not, recurse to populate left
  // child
  if (end_left > start) {
    cout << "Going left!\n";
    node->left_ = BinTree(input, start, end_left);
  }

  // check if we are done going to the right, if not, recurse to populate right
  // child
  if (start_right < end) {
    cout << "Going right!\n";
    node->right_ = BinTree(input, start_right, end);
  }

  return node;
}

bool BFS(Node *src, Node *dest) {
  if (src == dest)
    return true;

  if (src == nullptr || dest == nullptr)
    return false;

  list<Node *> q;
  q.push_back(src);

  while (q.size()) {
    Node *n = q.front();
    q.pop_front(); // pop node off front of the queue

    if ((n->left_ == dest) || (n->right_ == dest))
      return true;

    if (n->left_)
      q.push_back(n->left_);
    if (n->right_)
      q.push_back(n->right_);
  }
  return false;
}

// We know that n1 isn't a descendent of n2, or vice versa, so we can do the
// recursive search and guarantee that it terminates
Node *FCAHelper(Node *root, Node *n1, Node *n2) {

  bool n1_left = BFS(root->left_, n1);
  bool n2_right = BFS(root->right_, n2);

  if ((n1_left && n2_right) || (!n1_left && !n2_right))
    return root; // if n1 and n2 are in different subtrees, this is the FCA

  if (n1_left && !n2_right)
    return FCAHelper(
        root->left_, n1,
        n2); // if everyone is in left subtree find FCA in left subtree

  if (!n1_left && n2_right)
    return FCAHelper(root->right_, n1,
                     n2); // ditto but if everyone is in right subtree

  // if none of the above happened, this means that one of the nodes don't
  // exist, so return null ptr

  return nullptr;
}

// Before entering into recursive search to find FCA of (n1,n2), make sure
// n1 and n2 exist in the graph, and then make sure that n1
// isn't an ancestor of n2, and vice versa.
Node *FCA(Node *root, Node *n1, Node *n2) {

  if (!BFS(root, n1) || !BFS(root, n2))
    return nullptr;

  if (BFS(n1, n2)) {
    cout << "n1 ancestor of n2\n";
    return n1;
  }

  if (BFS(n2, n1))
    return n2;

  cout << "entering FCA\n";

  return FCAHelper(root, n1, n2);
}

int main() {
  const int N = 7;
  array<int, N> input{1, 2, 3, 4, 5, 18, 21};

  Node *root = BinTree(input, 0, N);

  Node *fca = FCA(root, root->right_->left_, root->right_->right_);

  cout << (fca ? fca->val_ : -1111) << "\n";

  return 0;
}
