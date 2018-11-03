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
inline std::vector<std::vector<T> > serialize(const Graph<T>& G)
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
    
    std::vector<std::vector<T> > v;
    std::set<T>    V;
    for (const T& N : G.nodes()) {
        std::vector<T> R;
        visit(G, N, V, R);
        std::sort(R.begin(),R.end());
        if (R.size()>0) v.push_back(R);
    }
    return v;
}


template <typename T>
inline std::vector<T> serialize(const Graph<T>& G, const T& N)
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
    visit(G, N, V, R);
    return R;
}

template <typename T>
inline std::vector<std::vector<T> > parallelize(const Graph<T>& g)
{
    //-----------------------------------------------------------
    // Find the level of a node n -> {m1,m2,...} such that
    //        level(n -> {})            = 0
    //        level(n -> {m1,m2,...})    = 1 + max(level(mi))
    //-----------------------------------------------------------
    typedef std::function<int(const Graph<T>& g, const T& n1, std::map<T, int>&)> Levelfun;
    
    Levelfun level = [&level](const Graph<T>& g, const T& n1, std::map<T, int>& levelcache) -> int {
        auto p = levelcache.find(n1);
        if (p != levelcache.end()) {
            return p->second;
        } else {
            int l = -1;
            for (const auto& e : g.edges(n1)) {
                l = std::max(l, level(g, e.first, levelcache));
            }
            return levelcache[n1] = l + 1;
        }
    };
    
    std::map<T, int> levelcache;
    // compute the level of each node in the graph
    int l = -1;
    for (const T& n : g.nodes()) {
        l = std::max(l, level(g, n, levelcache));
    }
    // create a graph for each level and place
    // each node in the appropriate level
    std::vector<std::vector<T> > v;
    v.resize(l + 1);
    for (const T& n : g.nodes()) {
        v[levelcache[n]].push_back(n);
    }
    
    return v;
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

