#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <list>
#include <stack>

class Graph
{
    int V;    // No. of vertices
    std::list<int> *adj;    // An array of adjacency lists
 
    // Fills Stack with vertices (in increasing order of finishing
    // times). The top element of stack has the maximum finishing 
    // time
    void fillOrder(int v, bool visited[], std::stack<int> &Stack);
 
    // A recursive function to print DFS starting from v
    void DFSUtil(int v, bool visited[]);
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

