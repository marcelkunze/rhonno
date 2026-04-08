// C++ Implementation of Kosaraju's algorithm to print all SCCs
// M.Kunze, Heidelberg University, 2018

#include <iostream>
#include <vector>
#include "Graph.h"
using namespace std;
  
Graph::Graph(int V)
: V(V), adj(V)
{
}
 
// A recursive function to print DFS starting from v
void Graph::DFSUtil(int v, std::vector<bool> &visited)
{
    // Mark the current node as visited and print it
    visited[v] = true;
    cout << v << " ";
 
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            DFSUtil(*i, visited);
}
 
Graph Graph::getTranspose()
{
    Graph g(V);
    for (int v = 0; v < V; v++)
    {
        for (const int w : adj[v])
        {
            g.adj[w].push_back(v);
        }
    }
    return g;
}
 
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
}
 
void Graph::fillOrder(int v, std::vector<bool> &visited, stack<int> &Stack)
{
    // Mark the current node as visited and print it
    visited[v] = true;
 
    // Recur for all the vertices adjacent to this vertex
    for (const int w : adj[v])
        if (!visited[w])
            fillOrder(w, visited, Stack);
 
    // All vertices reachable from v are processed by now, push v 
    Stack.push(v);
}
 
// The main function that finds and prints all strongly connected 
// components
int Graph::printSCCs()
{
    stack<int> Stack;
    int n = 0;
 
    // Mark all the vertices as not visited (For first DFS)
    std::vector<bool> visited(V, false);
 
    // Fill vertices in stack according to their finishing times
    for(int i = 0; i < V; i++)
        if(visited[i] == false)
            fillOrder(i, visited, Stack);
 
    // Create a reversed graph
    Graph gr = getTranspose();
 
    // Mark all the vertices as not visited (For second DFS)
    std::fill(visited.begin(), visited.end(), false);
 
    // Now process all vertices in order defined by Stack
    while (Stack.empty() == false)
    {
        // Pop a vertex from stack
        int v = Stack.top();
        Stack.pop();
 
        // Print Strongly connected component of the popped vertex
        if (visited[v] == false)
        {
            gr.DFSUtil(v, visited);
            cout << endl;
            n++;
        }
    }
    return n;
}
 
// Driver program to test above functions
/*
 int main()
{
    // Create a graph given in the above diagram
    Graph g(5);
    g.addEdge(1, 0);
    g.addEdge(0, 2);
    g.addEdge(2, 1);
    g.addEdge(0, 3);
    g.addEdge(3, 4);
 
    cout << "Following are strongly connected components in "
            "given graph \n";
    g.printSCCs();
 
    return 0;
}
*/
