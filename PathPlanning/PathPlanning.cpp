#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>

/* This project is inspired by the code found at
 * https://www.geeksforgeeks.org/shortest-path-unweighted-graph/ but goes well
 * beyond what's shown there.  I've implemented a Search Tree class that builds
 * up a Node based representation of a search tree as it performs either DFS,
 * BFS, or Double Sided BFS (DblBFS) to find and print the shortest path between
 * a source and destination node in a graph as specified by the adjacency list
 * encoded as an array of vectors, where adj[i] is a vector containing integers
 * that label a node i's neighbors*/

using namespace std;

struct Node {

  vector<Node *> nbrs_; // as defined by adjacency matrix
  int id_;              // just a label
  int dist_;            // distance from the source node
  Node *pred_;          // predecessor in our search algorithm
  bool visited_; // whether we've seen this node or not yet to prevent inf loops

  Node(int id, vector<Node *> nbrs = {}, int dist = numeric_limits<int>::max(),
       Node *pred = nullptr, bool visited = false)
      : id_{id}, nbrs_{nbrs}, dist_{dist}, pred_{pred}, visited_{visited} {
    // cout << "Creating Node " << id_ << "\n";
  }

  // Disallow copy construction
  Node(const Node &node) = delete;

  // Disallow copy assignment
  Node &operator=(const Node &node) = delete;

  ~Node() {} // cout << "Destroying node " << id_ << ".\n"; }
};

template <size_t V> using Adjacency = array<vector<int>, V>;

template <size_t V> class SearchTree {
private:
  unordered_map<int, Node *> nodes_;
  Adjacency<V> adj_;

public:
  SearchTree(int src_id, Adjacency<V> adj) : adj_{adj} {
    nodes_[src_id] = new Node(src_id);
  }

  Node *GetNode(int id) {
    if (nodes_.find(id) != nodes_.end())
      return nodes_[id];

    Node *new_node = new Node(id);
    nodes_[id] = new_node;
    return new_node;
  }

  void FillNeighbors(int id) {
    for (const auto &nbr_id : adj_[id])
      nodes_[id]->nbrs_.push_back(GetNode(nbr_id));
  }

  bool BFS(int src, int dest) {

    list<Node *> q;

    Node *src_node = GetNode(src);
    src_node->visited_ = true;
    src_node->dist_ = 0;
    FillNeighbors(src);

    // BFS initialization
    q.push_back(src_node);

    while (q.size()) {

      Node *node = q.front();
      q.pop_front(); // pop the first element of the top of the queue used to
                     // implement BFS.

      for (Node *nbr : node->nbrs_) {
        if (!(nbr->visited_)) {
          // standard BFS
          // cout << "Visiting node: " << nbr->id_ << "\n";
          nbr->visited_ = true;
          nbr->dist_ = node->dist_ + 1;
          nbr->pred_ = node;
          if (nbr->id_ == dest)
            return true;

          // if we didnt find our destination, let's keep going and populate
          // stuff
          FillNeighbors(nbr->id_);
          q.push_back(nbr);
        }
      }
    }
    return false;
  }

  // Run BFS from src to dest, and in reverse from dest to src, and wait until
  // we meet in the middle.  Improves worst case run time from O(b^d) to
  // O(b^{d/2})
  bool DblBFS(int src, int dest) {

    list<Node *> src_q;
    list<Node *> dest_q;

    typename SearchTree<V>::SearchTree rev_tree{dest, adj_};

    // Use this-> prefix here to distinguish between rev_tree to make sure we
    // are getting nodes and filling neighbors in the right tree
    Node *src_node = this->GetNode(src);
    Node *dest_node = rev_tree.GetNode(dest);

    src_node->visited_ = true;
    src_node->dist_ = 0;
    this->FillNeighbors(src);

    dest_node->visited_ = true;
    dest_node->dist_ = 0;
    rev_tree.FillNeighbors(dest);

    // BFS initialization
    src_q.push_back(src_node);
    dest_q.push_back(dest_node);

    while (src_q.size() && dest_q.size()) {

      // src BFS step
      Node *node = src_q.front();
      src_q.pop_front(); // pop the first element of the top of the queue used
                         // to implement BFS.

      for (Node *nbr : node->nbrs_) {
        if (!(nbr->visited_)) {
          // standard BFS
          // cout << "Visiting node: " << nbr->id_ << "\n";
          nbr->visited_ = true;
          nbr->dist_ = node->dist_ + 1;
          nbr->pred_ = node;

          this->FillNeighbors(nbr->id_);
          src_q.push_back(nbr);
        }
      }

      // check intersection of src_q and dest_q
      int middle = Intersect(src_q, dest_q);
      if (middle >= 0) {
        this->PrintPath(src, middle);
        rev_tree.PrintPath(dest, middle);
        return true;
      }

      // dest BFS step
      node = dest_q.front();
      dest_q.pop_front(); // pop the first element of the top of the queue used
                          // to implement BFS.

      for (Node *nbr : node->nbrs_) {
        if (!(nbr->visited_)) {
          // standard BFS
          // cout << "Visiting node: " << nbr->id_ << "\n";
          nbr->visited_ = true;
          nbr->dist_ = node->dist_ + 1;
          nbr->pred_ = node;

          rev_tree.FillNeighbors(nbr->id_);
          dest_q.push_back(nbr);
        }
      }

      // check intersection of src_q and dest_q
      middle = Intersect(src_q, dest_q);
      if (middle >= 0) {
        this->PrintPath(src, middle);
        rev_tree.PrintPath(dest, middle);
        return true;
      }
    }
    return false;
  }

  int Intersect(list<Node *> q1, list<Node *> q2) {
    auto comp = [](const Node *const &a, const Node *const &b) {
      return a->id_ > b->id_;
    };
    q1.sort(comp);
    q2.sort(comp);
    // sort(q2.begin(), q2.end(),
    //   [](const Node *a, const Node *b) { return a->id_ > b->id_; });
    for (auto it = q1.begin(); it != q1.end(); it++) {
      for (auto jt = q2.begin(); jt != q2.end(); jt++) {
        if ((*it)->id_ == (*jt)->id_)
          return (*it)->id_;
      }
    }
    return -1;
  }

  void DFS(int src, int dest, bool &success) {
    if (success)
      return;

    Node *src_node = GetNode(src);
    src_node->visited_ = true;
    FillNeighbors(src);

    for (const auto &nbr : src_node->nbrs_) {
      if (!(nbr->visited_) && !success) {
        nbr->visited_ = true;
        nbr->pred_ = src_node;

        if (nbr->id_ == dest) {
          success = true;
          return;
        }

        FillNeighbors(nbr->id_);
        DFS(nbr->id_, dest, success);
      }
    }
  }

  void PrintPath(int src, int dest) {
    vector<int> path;
    path.push_back(dest);
    Node *p = GetNode(dest)->pred_;

    while (p) {
      path.push_back(p->id_);
      p = p->pred_;
    }

    cout << "A path between them is given by ";
    for (auto it = path.end() - 1; it >= path.begin(); --it)
      cout << (*it) << " ";
    cout << ".\n";
  }

  void PrintPathFromNodeBFS(int src, int dest) {

    bool success = BFS(src, dest);
    if (!success) {
      cout << "There is no path between " << src << " and " << dest << ".\n";
      return;
    }

    cout << "Src node " << src << " is " << GetNode(dest)->dist_
         << " hops from dest node " << dest << ".\n";

    PrintPath(src, dest);
  }

  void PrintPathFromNodeDFS(int src, int dest) {
    bool success = false;
    DFS(src, dest, success);

    if (!success) {
      cout << "There is no path between " << src << " and " << dest << ".\n";
      return;
    }

    PrintPath(src, dest);
  }

  void Clear() {
    for (auto node : nodes_) {
      if (node.second)
        delete node.second;
    }

    nodes_.clear();
  }

  ~SearchTree() {
    // cout << "Search tree destructor\n";
    for (auto node : nodes_) {
      if (node.second)
        delete node.second;
    }
  }
};

// helper function to get a node's neighbors from adjacency matrix, this is a
// surrogate for a more realistic
template <size_t V>
vector<Node *> GetNeighbors(int id, const Adjacency<V> &adj,
                            map<int, Node *> &nodes_visited) {
  vector<Node *> nbrs;

  for (const auto &node : adj[id]) {
    if (nodes_visited.find(node) != nodes_visited.end())
      nbrs.push_back(nodes_visited[node]);
    else {
      // cout << "We are making a new neighbor node with id " << node << "\n";
      Node *new_nbr = new Node(node, {});
      nbrs.push_back(new_nbr);
      nodes_visited[node] = new_nbr;
    }
  }

  return nbrs;
}

template <size_t V> void AddEdge(Adjacency<V> &adj, int src, int dest) {
  adj[src].push_back(dest);
  adj[dest].push_back(src);
}

template <size_t V> void PrintAdjacency(const Adjacency<V> &adj) {
  for (int i = 0; i < V; ++i) {
    cout << "Node " << i << " has neighbors: ";
    for (const auto &nbr : adj[i])
      cout << nbr << " ";
    cout << "\n";
  }
}

template <size_t V>
bool BFS(const Adjacency<V> &adj, int src, int dest, array<int, V> &dist,
         array<int, V> &pred) {
  list<int> q;
  array<bool, V> visited;

  // initialize everything
  for (int i = 0; i < V; ++i) {
    dist[i] = numeric_limits<int>::max();
    pred[i] = -1;
    visited[i] = false;
  }

  // BFS initialization
  dist[src] = 0;       // src is 0 hops from src
  visited[src] = true; // we start at src, hence we've visited it;
  q.push_back(src);

  while (q.size()) {
    int node = q.front();
    q.pop_front(); // pop the first element of the top of the queue used to
                   // implement BFS.

    for (const auto &nbr : adj[node]) {
      if (!visited[nbr]) {
        // standard BFS
        visited[nbr] = true;
        q.push_back(nbr);
        // cout << "Visited node " << nbr << ".\n";
        // path planning related stuff
        dist[nbr] = dist[node] + 1;
        pred[nbr] = node;
        if (nbr == dest)
          return true;
      }
    }
  }

  return false;
}

// We only need the template and adjacency matrix because we need it to
// determine neighbors, nothing else
template <size_t V>
Node *NodeBasedBFS(const Adjacency<V> &adj, int src, int dest,
                   map<int, Node *> &nodes_visited) {

  list<Node *> q;

  Node *src_node =
      new Node(src, GetNeighbors(src, adj, nodes_visited), 0, nullptr, true);

  // BFS initialization
  q.push_back(src_node);
  nodes_visited[src] = src_node;

  // cout << "entering while loop\n";
  while (q.size()) {
    Node *node = q.front();
    q.pop_front(); // pop the first element of the top of the queue used to
                   // implement BFS.

    for (Node *nbr : node->nbrs_) {
      if (!(nbr->visited_)) {
        // standard BFS
        // cout << "Visiting node: " << nbr->id_ << "\n";
        nbr->visited_ = true;
        // nodes_visited[nbr->id_] = nbr;
        nbr->dist_ = node->dist_ + 1;
        nbr->pred_ = node;
        if (nbr->id_ == dest)
          return nbr;

        // if we didnt find our destination, let's keep going and populate
        // stuff

        nbr->nbrs_ = GetNeighbors(nbr->id_, adj, nodes_visited);

        q.push_back(nbr);
      }
    }
  }
  return nullptr;
}

template <size_t V>
void PrintShortestPath(const Adjacency<V> &adj, int src, int dest,
                       array<int, V> dist, array<int, V> pred) {
  bool success = BFS(adj, src, dest, dist, pred);
  if (!success) {
    cout << "No path exists between " << src << " and " << dest << ".\n";
    return;
  }

  // Let's back out the path that BFS found
  vector<int> path;

  path.push_back(dest);
  int p = pred[path.front()];
  // loop through list of node predecessors until we reach source, which has
  // pred[src] = -1;

  while (p != -1) {
    path.push_back(p);
    p = pred[p];
  }

  // Ok let's print everything nice and pretty
  cout << "Src node " << src << " is " << dist[dest] << " hops from dest node "
       << dest << ".\n";
  cout << "The path between them is given by ";
  for (auto it = path.end() - 1; it >= path.begin(); --it)
    cout << (*it) << " ";
  cout << ".\n";
}

void PrintPathFromNode(const Node *dest, int src) {
  vector<int> path;
  path.push_back(dest->id_);
  Node *p = dest->pred_;
  while (p) {
    path.push_back(p->id_);
    p = p->pred_;
  }
  cout << "Src node " << src << " is " << dest->dist_ << " hops from dest node "
       << dest->id_ << ".\n";
  cout << "The path between them is given by ";
  for (auto it = path.end() - 1; it >= path.begin(); --it)
    cout << (*it) << " ";
  cout << ".\n";
}

int main() {
  // no. of vertices
  const int v = 8;

  // array of vectors is used to store the graph
  // in the form of an adjacency list
  Adjacency<v> adj;

  // Creating graph given in the above diagram.
  // AddEdge function takes adjacency list, source
  // and destination vertex as argument and forms
  // an edge between them.
  AddEdge(adj, 0, 1);
  AddEdge(adj, 0, 3);
  AddEdge(adj, 1, 2);
  AddEdge(adj, 3, 4);
  AddEdge(adj, 3, 7);
  AddEdge(adj, 4, 5);
  AddEdge(adj, 4, 6);
  AddEdge(adj, 4, 7);
  AddEdge(adj, 5, 6);
  AddEdge(adj, 6, 7);

  PrintAdjacency(adj);

  int source = 5, dest = 2;

  array<int, v> dist, pred;

  PrintShortestPath(adj, source, dest, dist, pred);

  SearchTree<v> searchtree{source, adj};

  cout << "Calling BFS\n";
  searchtree.PrintPathFromNodeBFS(source, dest);

  searchtree.Clear();

  cout << "Calling DFS\n";
  searchtree.PrintPathFromNodeDFS(source, dest);

  searchtree.Clear();
  cout << "Calling DblBFS\n";

  searchtree.DblBFS(source, dest);

  /*map<int, Node *> nodes_visited;

    Node *destNode = NodeBasedBFS(adj, source, dest, nodes_visited);*/

  /*cout << "Did our node based implementation find a path? "
       << (destNode ? "yes!" : "no :(") << "\n";

  if (destNode)
    PrintPathFromNode(destNode, source);

  // memory management -- replace with shared_ptr

  for (auto &val : nodes_visited) {
    if (val.second)
      delete val.second;
      }*/

  return 0;
}
