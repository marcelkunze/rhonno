#ifndef _GRAPH_H_
#define _GRAPH_H_
// C++ Implementation of Kosaraju's algorithm to print all SCCs
// M.Kunze, Heidelberg University, 2018

#include <list>
#include <stack>
#include <vector>

class Graph
{
    int V;    // No. of vertices
    std::vector<std::list<int>> adj;    // Adjacency lists
 
    // Fills Stack with vertices (in increasing order of finishing
    // times). The top element of stack has the maximum finishing 
    // time
    void fillOrder(int v, std::vector<bool> &visited, std::stack<int> &Stack);
 
    // A recursive function to print DFS starting from v
    void DFSUtil(int v, std::vector<bool> &visited);
public:
    Graph(int V);
    void addEdge(int v, int w);
 
    // The main function that finds and prints strongly connected
    // components
    int printSCCs();
 
    // Function that returns reverse (or transpose) of this graph
    Graph getTranspose();
};

#endif

