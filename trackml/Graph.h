#ifndef _GRAPHh_H_
#define _GRAPHh_H_
// Implements weighted directed graph
// M.Kunze, Heidelberg University, 2018

#include <stdio.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stack>

template <typename T>
class Graph
{
private:
    std::set<T>                   fNodes;
    std::map<T, std::map<T, int>> fEdges;
public:
    
    void add(T n) // node
    {
        //if (fEdges.find(n)==fEdges.end()) return; // The node exists
        fNodes.insert(n);
        (void)fEdges[n];
    }
    
    void add(const T& n1, const T& n2, int d = 0) // edge
    {
        add(n1);
        add(n2);
        auto& adj = fEdges[n1];
        auto  n   = adj.find(n2);
        if (n != adj.end()) {
            int& d1 = n->second;
            if (d < d1) d1 = d;
        } else {
            adj[n2] = d;
        }
    }
    
    const std::set<T>& nodes() const
    {
        return fNodes;
    }
    
    const std::map<T, int>& edges(const T& n) const
    {
        static const std::map<T, int> null;
        if (fEdges.find(n)==fEdges.end()) return null; // The node does not exist
        return fEdges.at(n);
    }
    
    bool areConnected(const T& n1, const T& n2, int& d) const
    {
        if (fEdges.find(n1)==fEdges.end()) return false; // The node does not exist
        auto c = fEdges.at(n1);
        auto q = c.find(n2);
        if (q != c.end()) {
            d = q->second;
            return true;
        } else {
            return false;
        }
    }
    
    bool areConnected(const T& n1, const T& n2) const
    {
        int d;
        return areConnected(n1, n2, d);
    }
};

template <typename T>
inline std::vector<T> serialize(const Graph<T>& G)
{
    typedef std::function<void(const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R)> Visitfun;
    Visitfun visit = [&visit](const Graph<T>& G, const T& N, std::set<T>& V, std::vector<T>& R) {
        if (V.find(N) == V.end()) {
            V.insert(N);
            for (const auto& e : G.edges(N)) {
                visit(G, e.first, V, R);
            }
            R.push_back(N);
        }
    };
    
    std::vector<T> R;
    std::set<T>    V;
    for (const T& N : G.nodes()) {
        visit(G, N, V, R);
    }
    return R;
}

// print on stream
template <typename T>
inline std::ostream& operator<<(std::ostream& file, const Graph<T>& g)
{
    std::string sep      = "";
    bool   nok = false;
    
    file << "Graph {";
    for (const T& n : g.nodes()) {
        nok    = true;
        bool cok = false;
        if (g.edges(n).size()==0) continue;
        for (const auto& c : g.edges(n)) {
            cok = true;
            if (c.second == 0) {
                file << sep << n << "->" << (c.first);
            } else {
                file << sep << n << '-' << c.second << "->" << (c.first);
            }
            sep = ", ";
        }
        if (!cok) {
            file << sep << n;
        }
        sep = ", ";
    }
    
    return file << "}";
}

#endif

