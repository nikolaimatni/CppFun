#include <array>
#include <functional>
#include <iostream>
#include <list>
#include <vector>

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

/*this was an exercise in futility, as there's a much simpler way to do this
 * using a vector to represent the heap.  This was still good fun to play around
 * with the tree structure.  The motivation to use node based representations in
 * the other seetings (shortest path using BFS an DblBFS) was to avoid having to
 * hold the entire data-structure in memory in the event that there were local
 * rules to determine node's neighbors.  This motivation goes away for a
 * BinaryHeap because the it is a data-structure meant to be held entirely in
 * memory (e.g., to implement a priority queue) */
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

  if (n1_left)
    return FCAHelper(
        root->left_, n1,
        n2); // if everyone is in left subtree find FCA in left subtree

  if (n2_right)
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

template <typename Item> class ArrayHeap {
  using Compare = function<bool(const Item &, const Item &)>;

  // we leave position 0 empty, root is at 1, so then
  // node i's left child is at 2*i and 2*i + 1
  // node i's parent is i/2;

private:
  vector<Item> nodes_;
  Compare comp_;

public:
  ArrayHeap(Item root, Compare comp) : nodes_{vector<Item>(2)}, comp_(comp) {
    nodes_[1] = (root);
  }
  Item root() { return nodes_[1]; }
  Item get_node(int i) { return nodes_[i]; }

  void Insert(Item node) {
    nodes_.push_back(node);

    int i = nodes_.size() - 1;
    int i_parent = i / 2;

    while (i_parent > 0 && comp_(nodes_[i], nodes_[i_parent])) {
      // while we're not at the root and our parent precedes us, swap with
      // parent
      swap(nodes_[i], nodes_[i_parent]);
      i = i_parent;
      i_parent = i / 2;
    }
  }

  void ExtractRoot() {
    // swap root and last item, then delete last item
    swap(nodes_[1], nodes_.back());
    nodes_.pop_back();

    // now bubble down
    int i = 1;
    int left = 2;
    int right = 3;

    while (left < nodes_.size() || right < nodes_.size()) {
      // while we haven't reached the end of the graph yet
      bool goleft = comp_(nodes_[left], nodes_[i]) &&
                    left < nodes_.size(); // does left child precede i and is it
                                          // still in the graph
      bool goright = comp_(nodes_[right], nodes_[i]) &&
                     right < nodes_.size(); // does right child precede i and
                                            // is it still in the graph
      int smallest = 0;

      if (goleft)
        smallest = left;
      if (goright && comp_(nodes_[right], nodes_[left]))
        smallest = right;

      if (smallest) {
        swap(nodes_[i], nodes_[smallest]);
        i = smallest;
        left = 2 * i;
        right = 2 * i + 1;
      } else {
        return;
      }

      /*      if (goleft && goright) {
        // if both children precede, swap with minimum
        int gonext = comp_(nodes_[left], nodes_[right]) ? left : right;
        swap(nodes_[i], nodes_[gonext]);
        i = gonext;
        left = 2 * i;
        right = 2 * i + 1;
      } else if (goleft) {
        // if only left precedes, swap with left
        swap(nodes_[i], nodes_[left]);
        i = left;
        left = 2 * i;
        right = 2 * i + 1;
      } else if (goright) {
        // if only right precedes, swap with right
        swap(nodes_[i], nodes_[right]);
        i = right;
        left = 2 * i;
        right = 2 * i + 1;
      } else {
        // if neither children precedes, we're done!
        return;
        }*/
    }
  }

  int size() { return nodes_.size(); }
};

template <typename Item, typename Compare>
void HeapSort(vector<Item> &input, Compare comp) {
  ArrayHeap<Item> heap(input.back(), comp);
  input.pop_back();

  while (input.size()) {
    heap.Insert(input.back());
    input.pop_back();
  }

  while (heap.size() - 1) {
    input.push_back(heap.root());
    heap.ExtractRoot();
  }
}

int main() {
  const int N = 7;
  array<int, N> input{1, 2, 3, 4, 5, 18, 21};

  Node *root = BinTree(input, 0, N);

  Node *fca = FCA(root, root->right_->left_, root->right_->right_);

  cout << (fca ? fca->val_ : -1111) << "\n";

  auto comp = [](const int &a, const int &b) { return a < b; };
  ArrayHeap<int> heap(10, comp);

  auto print_heap = [&heap]() {
    for (int i = 1; i < heap.size(); i++)
      cout << heap.get_node(i) << " ";
    cout << "\n";
  };

  heap.Insert(5);
  heap.Insert(8);
  heap.Insert(4);
  heap.Insert(1);

  print_heap();

  heap.ExtractRoot();

  print_heap();

  heap.ExtractRoot();

  print_heap();

  vector<int> val{12, 23, 0, -12, 45, 1, 43, 8, 10, 4};
  HeapSort(val, [](const int a, const int b) { return a * a > b * b; });

  for (int i = 0; i < val.size(); i++)
    cout << val[i] << " ";
  cout << "\n";

  //  1 4 8 10 5
  //  4 5 8 10
  //  5 10 8

  /*HeapNode heap_root(10);
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

  cout << hr->val_ << " -l- " << hr->left_->val_ << "; -r- " <<
  hr->right_->val_
       << "\n";

  cout << "This should equal 1238 " << hr->val_ << hr->left_->val_
       << hr->left_->right_->val_ << hr->left_->right_->left_->val_ << "\n";

  cout << heap.last() << " is last's val before extraction \n";
  heap.ExtractRoot();
  cout << hr->val_ << "\n";
  cout << heap.last() << " is last's val after extraction \n";
  cout << "This shoudl equal 23816 " << hr->val_ << hr->left_->val_
       << hr->left_->right_->val_ << hr->left_->right_->left_->val_ << "\n";
  */

  return 0;
}
