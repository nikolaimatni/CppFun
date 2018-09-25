#include <array>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <vector>

/* This project is based on the code found at
 * https://www.geeksforgeeks.org/shortest-path-unweighted-graph/ I wanted to
 * play around with array and vector from STL, and interface them with
 * templates, hence the restriction to need to know the size of the graph at
 * compile time.  I could drop all of this by replacing the outer array with a
 * vector and reserving space of size V at creation, or using a raw array of
 * vector<int>'s.  Another limitation of this approach is that it requires O(V)
 * space to be pre-allocated, which isn't feasible for large state-spaces, so I
 * implemented a "node" based version of everything. */

using namespace std;

// using raw pointers because nodes are not responsible for memory management of
// other nodes
class Node {
public:
  vector<Node *> nbrs_;
  int id_;
  int dist_;
  Node *pred_;
  bool visited_;

  Node(int id, vector<Node *> nbrs, int dist = numeric_limits<int>::max(),
       Node *pred = nullptr, bool visited = false)
      : id_{id}, nbrs_{nbrs}, dist_{dist}, pred_{pred}, visited_{visited} {
    //  cout << "Creating a Node with id " << id_ << "\n";
  }

  Node(const Node &node) = delete;

  ~Node() { cout << "Calling Node destructor.\n"; }
};

template <size_t V> using Adjacency = array<vector<int>, V>;

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
        cout << "Visited node " << nbr << ".\n";
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

  cout << "entering while loop\n";
  while (q.size()) {
    Node *node = q.front();
    q.pop_front(); // pop the first element of the top of the queue used to
                   // implement BFS.

    for (Node *nbr : node->nbrs_) {
      if (!(nbr->visited_)) {
        // standard BFS
        cout << "Visiting node: " << nbr->id_ << "\n";
        nbr->visited_ = true;
        // nodes_visited[nbr->id_] = nbr;
        nbr->dist_ = node->dist_ + 1;
        nbr->pred_ = node;
        if (nbr->id_ == dest)
          return nbr;

        // if we didnt find our destination, let's keep going and populate stuff

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

  int source = 2, dest = 7;

  array<int, v> dist, pred;

  PrintShortestPath(adj, source, dest, dist, pred);

  map<int, Node *> nodes_visited;

  Node *destNode = NodeBasedBFS(adj, source, dest, nodes_visited);

  cout << "Did our node based implementation find a path? "
       << (destNode ? "yes!" : "no :(") << "\n";

  if (destNode)
    PrintPathFromNode(destNode, source);

  // memory management -- replace with shared_ptr

  for (auto &val : nodes_visited) {
    if (val.second)
      delete val.second;
  }

  return 0;
}
