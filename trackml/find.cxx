// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include <vector>
#include <map>

using namespace std;

// Look for seeding points using a KNN search and a neural network to identify hit pairs
std::vector<pair<int,float> > Tracker::findSeeds(int s, std::vector<int> &neighbours)
{
    vector<pair<int,float> > seed;
    if (assignment[s] !=  0) return seed;
    
    treePoint &p0 = points[s];
    
    // Generate seeding points
    for (auto it:neighbours)
    {
        if (assignment[it] != 0) continue;
        //        if (!checkTheta(s,it)) continue;
        //        if (!checkRadius(s,it)) continue;
        //        if (!checkDistance(s,it)) continue;
        
        double recall = checkTracklet(s,it); // Search for hit pairs
        if (recall > 0) {
            if (_verbose) cout << s << " " << it << ": R2 OK " << recall << endl;
            seed.push_back(make_pair(it,(float)recall));
            paths.add(s,it,1000*recall);
            paths.add(it,s,1000*recall);
            // Add double hits
            int twin = p0.twin();
            if (twin > 0) {
                recall = checkTracklet(twin,it);
                seed.push_back(make_pair(twin,(float)recall));
                paths.add(s,twin,1000*recall);
                paths.add(twin,s,1000*recall);
                assignment[twin] = assignment[s];
                
                if (_verbose) cout << "findSeeds: Added double hit " << twin << endl;
            }
        }
        else
            if (_verbose) cout << s << " " << it << ": R2 NOK " << recall << endl;
    }
    
    long size = seed.size();
    seedstotal += size;
    //seedsok += checkLabels(seed);
    
    return seed;
}


// Look for seeding points by hit pair combinations in the innnermost layers
vector<pair<int, int> > Tracker::findSeeds()
{
    const int n=5; // Seeding layer combinations
    const int start_list[6][2] = {{0,1}, {11,12}, {4,5}, {0,4}, {0,11}, {18,19}};
    
    static int ntrack(1);
    
    vector<pair<int, int> > pairs;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < PHIDIM; j++) {
            for (int k=0;k<THEDIM;k++) {
                int tube1 = start_list[i][0];
                for (auto a : tube[tube1][j][k]) {
                    if (assignment[a] != 0) continue;
                    int tube2 = start_list[i][1];
                    auto b = tube[tube2][j][k];
                    //auto b = knn[a][j][k];
                    if (b.size() == 0) continue;
                    vector<pair<int,float> > seed = findSeeds(a,b);
                    long n = seed.size();
                    if (n > 0) assignment[a] = ntrack++;
                    /*
                     if (n>0&&_verbose) {
                     cout << n << " seeds from " << a << " (Tube: " << tube1 << "->" << tube2 << "):" << endl;
                     for (auto it: seed) cout << it.first << "(" << it.second << ") ";
                     cout << endl;
                     }
                     */
                    // Make pairs
                    for (auto &it : seed) pairs.push_back(make_pair(a, it.first));
                }
            }
        }
    }
    
    if (_verbose) {
        cout << paths << endl;
    }
    
    return pairs;
}


//Find pairs using a neural network
vector<pair<int, int> > Tracker::findPairs() {
    
    const int n = 5;//How many pairs of layers to consider. Roughly proportional to run-time, and setting this to 30 gave practically the same score (less than 0.0002 reduction)
    const pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};
    
    static int ntrack(0);
    ntrack++;
    
    vector<pair<int, int> > pairs;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <PHIDIM; j++) {
            for (int k=0;k<THEDIM;k++) {
                for (auto a : tube[start_list[i].first][j][k]) {
                    if (assignment[a] != 0) continue;
                    //tracking.add(a);
                    for (auto b : tube[start_list[i].second][j][k]) {
                        if (assignment[b] != 0) continue;
                        double recall = checkTracklet(a,b);
                        if (recall > THRESHOLD2) {
                            assignment[a] = ntrack;
                            //assignment[b] = ntrack;
                            pairs.push_back(make_pair(a, b));
                            paths.add(a,b,recall*1000);
                            paths.add(b,a,recall*1000);
                        }
                    }
                }
            }
        }
    }
    return pairs;
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::findTriples(int p0, int p1, std::vector<int> &seed,std::vector<triple> &triples)
{
    triple t;
    t.x = p0;
    t.y = p1;
    
    for (auto &it : seed)
    {
        if (!checkTheta(p1,it)) continue;
        if (!checkRadius(p1,it)) continue;
        if (!checkDistance(p1,it)) continue;
        
        double recall = checkTracklet(p0,p1,it);
        if (recall > 0) {
            t.z = it;
            t.r = recall;
            triples.push_back(t);
            if (_verbose) cout << t.x << " " << t.y << " " << it << ": R3 OK " << recall << endl;
        }
        else
            if (_verbose) cout << t.x << " " << t.y << " " << it << ": R3 NOK " << recall << endl;
    }
    
    //if (_verbose) { cout << "triples " << p0.id() << ":"; for (auto &it:triples) cout << " (" << it.x << "," << it.y << "," << it.z << ":" << it.r << ")"; cout << endl; }
    
    return triples.size();
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::findTriples(vector<pair<int,int> > seed, vector<triple> &triples) {
    
    for (auto &it : seed) {
        int phi  = PHIFACTOR*(M_PI+points[it.second].phi());
        int the  = THEFACTOR*(M_PI+points[it.second].theta());
        int l = points[it.second].layer();
        vector<triple> trip;
        addHits(it.first,it.second,l,phi,the,trip);
        //addHitsCached(it.first,it.second,phi,trip);
        for (auto t: trip) {
            paths.add(t.y,t.z,1000*t.r);
            paths.add(t.z,t.y,1000*t.r);
            //assignment[t.z] = assignment[t.y]  = assignment[t.x];
        }
        triples.insert(triples.end(),trip.begin(),trip.end()); // append the candidates
    }
    if (_verbose) cout << paths << endl;
    return triples.size();
}


// Generate tracklets of 3 points wrt. the first point in seed
long Tracker::addHits(int p0,int p1,int start,int phi,int the,std::vector<triple> &triples)
{
    static const std::map<int,vector<int> > layers{
        {0,{1,4,11}},
        {1,{2,4,11}},
        {2,{3,18}},
        {3,{18}},
        {4,{5}},
        {5,{6,24}},
        {6,{7}},
        {7,{8,26}},
        {8,{9}},
        {9,{10,30}},
        {10,{34}},
        {11,{12,18}},
        {12,{13}},
        {13,{14,16,38}},
        {14,{15,38}},
        {15,{16,40}},
        {16,{17}},
        {17,{46}},
        {18,{19,24,36}},
        {19,{20,24,36}},
        {20,{21,24,36}},
        {21,{22,25,27,37}},
        {22,{23,25,37}},
        {23,{-1}},
        {24,{26}},
        {25,{27}},
        {26,{28}},
        {27,{29}},
        {28,{30,31,33}},
        {29,{31}},
        {30,{32}},
        {31,{33}},
        {32,{34}},
        {33,{35}},
        {34,{-1}},
        {35,{-1}},
        {36,{38,41}},
        {37,{39}},
        {38,{40}},
        {39,{41}},
        {40,{42}},
        {41,{43}},
        {42,{44,47}},
        {43,{45}},
        {44,{46}},
        {45,{47}},
        {46,{-1}},
        {47,{-1}},
    };
    
    triple t;
    t.x = p0;
    t.y = p1;
    float r = radius(p1);
    
    if (start<0 || start>47) return 0;
    const auto laylist = layers.at(start);
    
    long found(0);
    for (int i=0; i<laylist.size(); i++) {
        int nextlayer = laylist[i];
        if (nextlayer<0 || nextlayer>47) break;
        //cout << "dephth " << dephth << " layer " << l << endl;
        auto &seed1 = tube[nextlayer][phi][the];
        for (auto &it1 : seed1)
        {
            if (assignment[it1] != 0) continue; // Point has benn already used
            //if (!checkTheta(p1,it)) continue;
            //if (!checkRadius(p1,it)) continue;
            float d = distance(p1,it1);
            float dr = distance(p1,it1)*r;
            if (dr>DISTANCE) { nd++; continue; }
            //if (!checkDistance(p1,it)) continue;
            
            double recall = checkTracklet(p0,p1,it1); // Point is a candidate on the next layer
            
            if (recall > 0) {
                t.z = it1;
                t.r = recall;
                triples.push_back(t);
                assignment[it1] = assignment[p0];
                // Add double hits
                int twin = points[t.z].twin();
                if (twin > 0) {
                    int s = t.z;
                    t.z = twin;
                    triples.push_back(t);
                    paths.add(s,twin,1000*t.r);
                    //paths.add(twin,s,1000*t.r);
                    assignment[twin] = assignment[s];
                    if (_verbose) cout << "addHits: Added double hit " << twin << endl;
                }
                found +=  addHits(p1,it1,nextlayer,phi,the,triples);
                if (_verbose) cout << t.x << " " << t.y << " " << it1 << ": R3 OK " << recall << endl;
            }
            else
                if (_verbose) cout << t.x << " " << t.y << " " << it1 << ": R3 NOK " << recall << endl;
        }
        
    }
    
    return found;
}
