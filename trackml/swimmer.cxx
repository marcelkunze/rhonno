// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include <vector>
#include <map>

using namespace std;

// Reconstruct tracks with a swimmer
map<int,vector<int> > Tracker::swimmer() {
    long nhits = points.size();
    map<int,vector<int> > tracklet;
    map<int,vector<int> > shortpath;
    map<int,treePoint*> hitmap;
    
    // Transfer the results from the weighted directed graph into hits
    for (auto &ip : points) {
        if (ip.id()==0) continue; // skip index 0
        int a = ip.id();
        int n = 0;
        map<int,float> const &path = tracking.edges(a);
        if (path.size()==0) continue;
        
        for (auto &it : path ) {
            //if (n>NEIGHBOURS) break;
            int id = it.first;
            float recall = it.second;
            ip.setneighbour(id,recall);
            points[a].setneighbour(id,recall);
            n++;
        }
    }
    
    if (verbose) {
        for (int i=1;i<nhits;i++) {
            auto adj = points[i].neighbours();
            cout << points[i].id() << " {" ;
            for (auto it : adj) cout << it << "," ;
            cout << "}" << endl;
        }
        cout << endl;
    }
    
    // fill the hitmap
    for (auto &p : points) {
        if (p.id()==0) continue; // skip index 0
        int id = p.id();
        hitmap[id] = &p;
    }
    
    // Swimmer
    
    int ntrack(1), nshort(1);
    while (hitmap.size()>0) {
        int n = 0;
        int neighbour = -1;
        vector<treePoint> pvec;
        auto it = hitmap.begin();
        treePoint *p0 = it->second;
        pvec.push_back(*p0);
        hitmap.erase(it++);
        while (it != hitmap.end()) { // Follow the path until there are no more neighbours
            neighbour = p0->neighbour(n);
            if (verbose) cout << p0->id() << "->" << neighbour << endl;
            if (neighbour < 0 || neighbour > nhits) break;
            auto it = hitmap.find(neighbour);
            if (it==hitmap.end()) { // hit is already assigned
                if (verbose) cout <<  "->" << neighbour << endl;
                n++;  // try an alternative neighbour
                continue;
            }
            if (it==hitmap.end()) break;
            treePoint *p1 = it->second;  // copy the point into the tracklet vector
            pvec.push_back(*p1);
            hitmap.erase(it++);
            n = 0;
            p0 = p1;
        }
        
        sort(pvec.begin(), pvec.end(), Point::sortById); // Sort the hits acording to the Id
        vector<int> tmpvec;
        for (auto &ip : pvec) tmpvec.push_back(ip.id()); // Note the hit indices
        if (verbose) print(tmpvec);
        if (pvec.size() >= TRACKLET) {
            tracklet[ntrack++] = tmpvec;
            trackletstotal += pvec.size();
            trackletsok += checkLabels(tmpvec);
        }
        else
            shortpath[nshort++] = tmpvec;
        if (verbose) cout << "---" << endl;
    }
    
    return tracklet;
}

