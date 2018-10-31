/*******************************************************************************
********************************************************************************

    digraph : directed graph

    Created by Yann Orlarey on 31/01/2017.
    Copyright © 2017 Grame. All rights reserved.

 *******************************************************************************
 ******************************************************************************/

#ifndef digraph_hh
#define digraph_hh

#include <stdio.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stack>

//===========================================================
// digraph : a directed graph, a set of nodes f type N and a
// set of connections between these nodes. Connections have an
// associated value, by defaut 0. This value is used in Faust
// to represent the time dependency between computations.
//===========================================================

template <typename N>
class digraph
{
   protected:
    //--------------------------------------------------------------------------
    // Real/internal structure of a graph. A graph is a set of nodes
    // and a set of connections between theses nodes. These connections
    // have an integer value attached.
    class internalgraph
    {
       private:
        std::set<N>                   fNodes;        // {n1,n2,...}
        std::map<N, std::map<N, int>> fConnections;  // {(ni -d-> nj),...}

       public:
#if 0
        internalgraph() { std::cout << "create internalgraph " << this << std::endl; }
        ~internalgraph() { std::cout << "delete internalgraph " << this << std::endl; }
#endif
        // Add the node n to the graph
        void add(N n)
        {
            if (fConnections.find(n)==fConnections.end()) return; // The node exists
            fNodes.insert(n);
            (void)fConnections[n];  // make sure we have an empty set of connections for n
        }

        // Add the nodes n1 and n2 and the connection (n1 -d-> n2) to the graph.
        // If a connection (n1 -d'-> n2) already exists, the connection is updated
        // with the min(d,d')
        void add(const N& n1, const N& n2, int d = 0)
        {
            add(n1);
            add(n2);
            auto& adj = fConnections[n1];
            auto  cnx = adj.find(n2);
            if (cnx != adj.end()) {
                int& d1 = cnx->second;
                if (d < d1) d1 = d;
            } else {
                adj[n2] = d;
            }
        }

        // returns the set of nodes of the graph
        const std::set<N>& nodes() const
        {
            return fNodes;
        }

        // returns the connections of node n in the graph
        const std::map<N, int>& connections(const N& n) const
        {
            static const std::map<N, int> null;
            if (fConnections.find(n)==fConnections.end()) return null; // The node does not exist
            return fConnections.at(n);
        }

        // tests if two nodes are connected
        bool areConnected(const N& n1, const N& n2, int& d) const
        {
            if (fConnections.find(n1)==fConnections.end()) return false; // The node does not exist
            auto c = fConnections.at(n1);
            auto q = c.find(n2);
            if (q != c.end()) {
                d = q->second;
                return true;
            } else {
                return false;
            }
        }

        // tests if two nodes are connected
        bool areConnected(const N& n1, const N& n2) const
        {
            int d;
            return areConnected(n1, n2, d);
        }
    };

    std::shared_ptr<internalgraph> fContent;

   public:
    digraph() : fContent(new internalgraph)
    {
    }

    ~digraph()
    {
    }
    
    // build the graph

    digraph& add(N n)
    {
        fContent->add(n);
        return *this;
    }

    // Add a graph with all its connections
    digraph& add(const digraph& g)
    {
        for (auto& n : g.nodes()) {
            add(n);
            for (auto& c : g.connections(n)) {
                add(n, c.first, c.second);
            }
        }
        return *this;
    }

    digraph& add(const N& n1, const N& n2, int d = 0)
    {
        fContent->add(n1, n2, d);
        return *this;
    }

    // query the graph

    const std::set<N>& nodes() const
    {
        return fContent->nodes();
    }
    const std::map<N, int>& connections(const N& n) const
    {
        return fContent->connections(n);
    }

    bool areConnected(const N& n1, const N& n2, int& d) const
    {
        return fContent->areConnected(n1, n2, d);
    }
    bool areConnected(const N& n1, const N& n2) const
    {
        return fContent->areConnected(n1, n2);
    }

    // compare graphs for maps and other containers

    friend bool operator<(const digraph& p1, const digraph& p2)
    {
        return p1.fContent < p2.fContent;
    }
    friend bool operator==(const digraph& p1, const digraph& p2)
    {
        return p1.fContent == p2.fContent;
    }
};

//===========================================================
//===========================================================
// file << graph : print graph on a stream
//===========================================================
//===========================================================

template <typename N>
inline std::ostream& operator<<(std::ostream& file, const digraph<N>& g)
{
    std::string sep      = "";
    bool   hasnodes = false;
    
    file << "Graph {";
    for (const N& n : g.nodes()) {
        hasnodes    = true;
        bool hascnx = false;
        if (g.connections(n).size()==0) continue;
        for (const auto& c : g.connections(n)) {
            hascnx = true;
            if (c.second == 0) {
                file << sep << n << "->" << (c.first);
            } else {
                file << sep << n << '-' << c.second << "->" << (c.first);
            }
            sep = ", ";
        }
        if (!hascnx) {
            file << sep << n;
        }
        sep = ", ";
    }
    
    return file << "}";
}

#endif /* digraph_hpp */

