// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include <vector>
#include <map>

using namespace std;

// Driver function to sort the vector elements
// by second element of pairs
bool sortbysec(const pair<int,float> &a,
               const pair<int,float> &b)
{
    return (a.second < b.second);
}

// Look for seeding points by hit pair combinations in the innnermost layers
long Tracker::findSeeds()
{
    const int n=4; // Seeding layer combinations
    const int start_list[4][2] = {{0,1}, {11,12}, {4,5}, {18,19}};
    
    static int ntrack(1);
    
    pairs.clear();

    for (int i = 0; i < n; i++) {
        
        int layer1 = start_list[i][0];
        for (auto start : modules[layer1]) { // all modules in first layer
            const auto edgelist = paths.edges(start);
            if (edgelist.size() == 0) continue;
            int l1 = start/MODULES;
            int m1 = start%MODULES;
            if (_verbose) cout << "Start layer " << l1 << " module " << m1 << endl;
            for (auto edge : edgelist) {
                int nextindex = edge.first;
                if (nextindex<0) break;
                int l2 = nextindex/MODULES;
                int m2 = nextindex%MODULES;
                if (_verbose) cout << "-> Layer " << l2 << " module " << m2 << endl;
                for (auto a : module[start]) { // all hits in module
                    
                    if (assignment[a] > 0) continue;
                    
                    int twin = points[a].twin();
                    if (twin>0) {
                        pairs.push_back(make_pair(a,twin));
                        tracking.add(a,twin,1.0);
                        assignment[a] = ntrack;
                        assignment[twin] = ntrack;
                        if (_verbose) cout << "-> Add twin " << a << " " << twin << endl;
                        continue;
                    }
                    
                    // Generate seeding points
                    auto b = module[nextindex]; // all hits in following modules
                    if (b.size() == 0) continue;
                    vector<pair<int,float> > seed;
                    for (auto it:b)
                    {
                        if (assignment[it] > 0) continue;
                        double recall = checkTracklet(a,it); // Search for hit pairs
                        if (recall < THRESHOLD2) {
                            if (_verbose) cout << a << " " << it << ": R2 NOK " << recall << endl;
                            continue;
                        }
                        float d = distance(a,it);
                        seed.push_back(make_pair(it,d));
                    }
                    if (seed.size()==0) continue;
                    
                    // Sort seeds according to their value
                    sort(seed.begin(),seed.end(),sortbysec);
                    
                    auto s = seed.begin(); // Take the seed with the shortest distance
                    int id = s->first;
                    float d = s->second;
                    float recall = checkTracklet(a,id);
                    
                    //if (!checkTheta(a,id)) continue;
                    //if (!checkPhi(a,id)) continue;
                    if (_verbose) cout << a << " " << id << ": R2 OK " << recall << " d " << d << endl;
                    pairs.push_back(make_pair(a,id));
                    assignment[id] = ntrack;
                    assignment[a] = ntrack++;
                    tracking.add(a,id,recall);
                }
            }
        }
        
    }
    
    // Print the tracking graph
    if (_verbose) tracking.print();
    
    sort(pairs.begin(),pairs.end());
         
    return pairs.size();
}

// Look for seeding points by hit pair combinations in the innnermost layers
long Tracker::findSeedsPhiTheta()
{
    const int n=5; // Seeding layer combinations
    const int start_list[6][2] = {{0,1}, {11,12}, {4,5}, {0,4}, {0,11}, {18,19}};
    
    static int ntrack(1);

    pairs.clear();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < PHIDIM; j++) {
            for (int k=0;k<THEDIM;k++) {
                int tube1 = start_list[i][0];
                
                for (auto a : tube[tube1][j][k]) {
                    if (assignment[a] > 0) continue;
                    int tube2 = start_list[i][1];
                    auto b = tube[tube2][j][k];
                    if (b.size() == 0) continue;
                    
                    vector<pair<int,float> > seed;
                    for (auto it:b)
                    {
                        if (assignment[it] > 0) continue;
                        double recall = checkTracklet(a,it); // Search for hit pairs
                        seed.push_back(make_pair(it,(float)recall));
                    }
                    if (seed.size()==0) continue;
                    
                    // Sort seeds according to their value
                    sort(seed.begin(),seed.end(),sortbysec);
                    
                    auto s = seed.begin(); // Take the seed with the highest value
                    int id = s->first;
                    float recall = s->second;
                    if (recall < THRESHOLD2) continue;
                    if (!checkTheta(a,id)) continue;
                    if (!checkPhi(a,id)) continue;
                    double d = distance(a,id);
                    if (_verbose) cout << a << " " << id << ": R2 OK " << recall << " d " << d << endl;
                    pairs.push_back(make_pair(a,id));
                    assignment[id] = ntrack;
                    assignment[a] = ntrack++;
                    tracking.add(a,id,recall);
                }
            }
        }
    }
    
    if (_verbose) {
        tracking.print();
    }
    
    return pairs.size();
}


//Find pairs using a neural network
long Tracker::findPairs() {
    
    const int n = 50;//How many pairs of layers to consider. Roughly proportional to run-time, and setting this to 30 gave practically the same score (less than 0.0002 reduction)
    const pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    static int ntrack(0);
    ntrack++;

    pairs.clear();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <PHIDIM; j++) {
            for (int k=0;k<THEDIM;k++) {
                for (auto a : tube[start_list[i].first][j][k]) {
                    if (assignment[a] > 0) continue;
                    for (auto b : tube[start_list[i].second][j][k]) {
                        if (assignment[b] != 0) continue;
                        double recall = checkTracklet(a,b);
                        if (recall > THRESHOLD2) {
                            assignment[a] = ntrack;
                            //assignment[b] = ntrack;
                            pairs.push_back(make_pair(a, b));
                            tracking.add(a,b,recall);
                            tracking.add(b,a,recall);
                        }
                    }
                }
            }
        }
    }
    return pairs.size();
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::findTriples(int p0, int p1, std::vector<int> &seed)
{
    triple t;
    t.x = p0;
    t.y = p1;
    
    //triples.clear();
    
    for (auto &it : seed)
    {
        //if (!checkTheta(p1,it)) continue;
        //if (!checkRadius(p1,it)) continue;
        //if (!checkDistance(p1,it)) continue;
        
        double recall = checkTracklet(p0,p1,it);
        if (recall > 0) {
            t.z = it;
            t.r = recall;
            triples.push_back(t);
            if (_verbose) cout << t.x << " " << t.y << " " << it << ": R3 OK" << recall << ", ";
        }
        else
            if (_verbose) cout << t.x << " " << t.y << " " << it << ": R3 NOK" << recall << ", ";
    }
    
    //if (_verbose) { cout << "triples " << p0.id() << ":"; for (auto &it:triples) cout << " (" << it.x << "," << it.y << "," << it.z << ":" << it.r << ")"; cout << endl; }
    
    return triples.size();
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::findTriples() {

    triples.clear();

    for (int i=0;i<MAXDIM;i++) assignment[i] = 0;
    
    static int ntrack = 1;
    for (auto &it : pairs) {
        if (assignment[it.first]>0) continue;
        assignment[it.first] = ntrack++;
        int l = points[it.second].layer();
        int m = points[it.second].module();
        int index = MODULES*l + m;
        if (_verbose) cout << endl << "findTripes " <<  "{" << index << ","  << l << "," << m << "/" << it.first << "," << it.second << "}" << endl;
        vector<triple> trip;
        addHits(it.first,it.second,index,trip);
        for (auto t: trip) {
            tracking.add(t.y,t.z,t.r);
            //  tracking.add(t.z,t.y,t.r);
            //assignment[t.z] = assignment[t.y]  = assignment[t.x];
        }
        triples.insert(triples.end(),trip.begin(),trip.end()); // append the candidates
    }
    if (_verbose) tracking.print();
    return triples.size();
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::addHits(int p0,int p1,int start,std::vector<triple> &triples)
{
    
    triple t;
    t.x = p0;
    t.y = p1;
    
    if (start<0 || start>MODULES*LAYERS) return 0;
    const auto edgelist = paths.edges(start); // Get the list of possible connected modules
    if (edgelist.size() == 0) return 0;
    int l = start/MODULES;
    int m = start%MODULES;
    if (_verbose) cout << "{" << start << ","  << l << "," << m << "} -> ";
    
    long found(0);
    for (auto edge : edgelist) {
        int nextindex = edge.first;
        int nextlayer = nextindex/MODULES;
        int nextmodule = nextindex%MODULES;
        if (nextmodule<0 || nextmodule>MODULES*LAYERS) break;
        if (_verbose) cout << "{" << nextindex << "," << nextlayer << "," << nextmodule << "},";
        auto &seed1 = module[nextindex];
        for (auto &it1 : seed1)
        {
            if (assignment[it1] > 0) continue; // Point has benn already used
            //if (!checkTheta(p1,it1)) continue;
            Point &a = points[p0];
            Point &b = points[p1];
            Point &c = points[it1];
            float d = Point::distance3(a,b,c); // distance of it from line p0-p1
            if (d>DISTANCE) {
                if (_verbose) cout << it1 << " -> distance " << d << endl;
                nd++; continue;
            }

            double recall = checkTracklet(p0,p1,it1); // Point is a candidate on the next layer
            //double recall = scoreTriple(p0,p1,it1); // Point is a candidate on the next layer
            
            if (recall > THRESHOLD3) {
                t.z = it1;
                t.r = recall;
                triples.push_back(t);
                assignment[it1] = assignment[p0];
                found +=  addHits(p1,it1,nextindex,triples);
                if (_verbose) cout << endl << t.x << " " << t.y << " " << it1 << ": R3 OK " << recall << ", ";
            }
            else
                if (_verbose) cout << endl << t.x << " " << t.y << " " << it1 << ": R3 NOK " << recall << ", ";
        }
    }
    
    return found;
}
