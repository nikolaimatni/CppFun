#include <array>
#include <iostream>
#include <list>

using namespace std;

struct Node {
  int val_;
  Node *left_;
  Node *right_;

  Node(int val) : val_(val), left_(nullptr), right_(nullptr) {}
};

struct HeapNode {
  int val_;
  HeapNode *left_;
  HeapNode *right_;
  HeapNode *parent_;

  HeapNode(int val)
      : val_(val), left_(nullptr), right_(nullptr), parent_(nullptr) {}
};

class Heap {
private:
  HeapNode *root_;
  HeapNode *insert_here_;
  HeapNode *last_;
  bool min_;

public:
  Heap(HeapNode *root, bool min_heap = true)
      : root_(root), insert_here_(root), last_(root), min_(min_heap) {}

  int last() { return last_->val_; }
  void InsertNode(HeapNode *node) {

    last_ = node; // note we don't shuffle nodes around, just their values, so
                  // this sets last_ pointing to the last position we inserted
                  // a node

    node->parent_ = insert_here_;

    if (insert_here_->left_) {
      insert_here_->right_ = node;
      UpdateInsertionNode();
      cout << "Updating insertion node to " << insert_here_->val_ << "\n";
    } else {
      insert_here_->left_ = node;
    }
    auto comp = [this](const HeapNode *const &a, const HeapNode *const &b) {
      return (min_ ? (a->val_ < b->val_) : (a->val_ > b->val_));
    };

    while ((node->parent_) && comp(node, node->parent_)) {
      SwapWithParent(node); // just swap the values
      node = node->parent_; // update the pointer
    }
  }

  void SwapWithParent(HeapNode *node) {
    int temp_val = node->val_;
    node->val_ = node->parent_->val_;
    node->parent_->val_ = temp_val;
  }

  HeapNode *LeftMostNode(HeapNode *root) {
    HeapNode *temp = root;
    while (temp->left_)
      temp = temp->left_;

    return temp;
  }

  HeapNode *RightMostNode(HeapNode *root) {
    HeapNode *temp = root;
    while (temp->right_)
      temp = temp->right_;

    return temp;
  }

  void UpdateInsertionNode() {
    // if we're calling this, this means that insert_here_ is full,
    // so let's first go up one level if we're not at root

    if (insert_here_ == root_) {
      insert_here_ = root_->left_;
      return;
    }

    HeapNode *up = insert_here_->parent_;

    cout << "up value is " << up->val_ << "\n";
    // climb up the tree until we aren't a right node of the parent
    // this indicates the right point to switch to a subtree that will have free
    // spots in its leftmost position
    while (insert_here_ == up->right_) {
      insert_here_ = up;
      if (up->parent_)
        up = up->parent_;
      else {
        // if we hit the root and we're still in the loop, it means its time to
        // add a new layer so insertion point is the left-most node in the tree
        insert_here_ = LeftMostNode(root_);

        return;
      }
    }

    // up indicates the switch point, so go one hop to the right, and then all
    // the way down to the left
    insert_here_ = LeftMostNode(up->right_);
  }

  HeapNode *root() { return root_; }

  void UpdateLast() {
    HeapNode *temp = last_;
    // do stuff
    HeapNode *up = last_->parent_;

    if (up && last_ == up->right_) {
      // if last was a right-node, then the update is easy
      last_ = up->left_;
      delete temp;
      return;
    }

    // otherwise we need to find the right subtree to go right on so
    // climb up the tree until we find the right spot to split
    while (up && last_ != up->right_) {
      if (up->parent_) {
        last_ = up;
        up = last_->parent_;
      } else {
        // if we hit the root without finding a match, then this is a corner
        // case for when we need to reset to the far right node
        last_ = RightMostNode(root_);
        delete temp;
        return;
      }
    }

    // otherwise take one step left from up to get into the subtree and go all
    // the way right
    last_ = RightMostNode(up->left_);

    delete temp;
  }

  HeapNode *ExtractRoot() {
    root_->val_ = last_->val_;
    // Update last (this also deletes the node we just swapped out)
    UpdateLast();

    // now bubble down
    HeapNode *node = root_;

    auto comp = [this](const HeapNode *const &a, const HeapNode *const &b) {
      return (min_ ? (a->val_ < b->val_) : (a->val_ > b->val_));
    };

    while (node->left_ || node->right_) {
      // while there are still children to swap with
      if (node->left_ && node->right_) {
        bool lcomp = comp(node->left_, node);
        bool rcomp = comp(node->right_, node);
        if (lcomp && rcomp) {
          // if both children can be swapped, swap with the extremal of the two
          HeapNode *swap =
              comp(node->left_, node->right_) ? node->left_ : node->right_;

          SwapWithParent(swap); // just swap the vals
          node = swap;          // now update the pointer

        } else if (lcomp || rcomp) {
          HeapNode *swap = lcomp ? node->left_ : node->right_;
          SwapWithParent(swap); // just swaps values;
          node = swap;          // now update the pointer

        } else {
          // none of our children come before us in the ordering, so we're done
          return root_;
        }
      }
    }

    return root_;
  }
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

  if ((n1_left && n2_right))
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

  HeapNode heap_root(10);
  Heap heap(&heap_root);

  for (int i = 9; i > 0; --i)
    heap.InsertNode(new HeapNode(i));

  heap.InsertNode(new HeapNode(11));
  heap.InsertNode(new HeapNode(12));
  heap.InsertNode(new HeapNode(13));
  heap.InsertNode(new HeapNode(14));
  heap.InsertNode(new HeapNode(15));
  heap.InsertNode(new HeapNode(16));

  HeapNode *hr = heap.root();

  cout << hr->val_ << " -l- " << hr->left_->val_ << "; -r- " << hr->right_->val_
       << "\n";

  cout << "This should equal 1238 " << hr->val_ << hr->left_->val_
       << hr->left_->right_->val_ << hr->left_->right_->left_->val_ << "\n";

  cout << heap.last() << " is last's val before extraction \n";
  heap.ExtractRoot();
  cout << hr->val_ << "\n";
  cout << heap.last() << " is last's val after extraction \n";
  cout << "This shoudl equal 23816 " << hr->val_ << hr->left_->val_
       << hr->left_->right_->val_ << hr->left_->right_->left_->val_ << "\n";

  return 0;
}
