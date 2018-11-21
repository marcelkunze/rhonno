// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
#include "Point.h"
#include "XMLP.h"
#include <limits>
#include <vector>
#include <map>

using namespace std;

bool
Tracker::intersection(int A, int B, int C, int D, Point& ip)
// http://mathworld.wolfram.com/Line-LineIntersection.html
// in 3d; will also work in 2d if z components are 0
{
    Point &a = points[A];
    Point &b = points[B];
    Point &c = points[C];
    Point &d = points[D];

    Point da = b - a;
    Point db = d - c;
    Point dc = c - a;
    
    if (Point::dot(dc, Point::cross(da,db)) != 0.0) // lines are not coplanar
        return false;
    
    double s = Point::dot(Point::cross(dc, db), Point::cross(da, db)) / Point::norm2(Point::cross(da, db));
    double t = Point::dot(Point::cross(dc, da), Point::cross(da, db)) / Point::norm2(Point::cross(da, db));
    if ((s >= 0 && s <= 1) && (t >= 0 && t <= 1))
    {
        ip = Point(a.x()+da.x()*s,a.y()+da.y()*s,a.z()+da.z()*s);
        return true;
    }
    
    return false;
}

Point Tracker::intersection2d(int a, int b, int c, int d)
{
    treePoint &A = Tracker::points[a];
    treePoint &B = Tracker::points[b];
    treePoint &C = Tracker::points[c];
    treePoint &D = Tracker::points[d];

    // Line AB represented as a1x + b1y = c1
    double a1 = B.y() - A.y();
    double b1 = A.x() - B.x();
    double c1 = a1*(A.x()) + b1*(A.y());
    
    // Line CD represented as a2x + b2y = c2
    double a2 = D.y() - C.y();
    double b2 = C.x() - D.x();
    double c2 = a2*(C.x())+ b2*(C.y());
    
    double determinant = a1*b2 - a2*b1;
    
    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return Point(1.E7, 1.E7, 1.E7);
    }
    else
    {
        double x = (b2*c1 - b1*c2)/determinant;
        double y = (a1*c2 - a2*c1)/determinant;
        double z1 = A.z() + (B.z()-A.z())/(x-A.x());
        double z2 = A.z() + (B.z()-A.z())/(y-A.y());
        double z = 0.5*(z1+z2);
        return Point(x,y,z);
    }
}

// Recall function for pairs
double Tracker::checkTracklet(int p0,int p1)
{
    double recall = recallPair(points[p0],points[p1])[0];
    bool ok = recall>THRESHOLD2;
    recall = ok ? recall : -recall;
    if (ok) n4++;
    n2++;
    return recall;
}


// Recall function for triple
double Tracker::checkTracklet(int p1,int p2,int p3)
{
    treePoint &point1 = points[p1];
    treePoint &point2 = points[p2];
    treePoint &point3 = points[p3];
    /*
     Point v1 = point2 - point1;
     Point v2 = point3 - point1;
     double angle = acos(Point::dot(v1,v2)); // Check angle between the last points
     if (angle<0.5*M_PI) {
     if (angle>DELTANN) { np++; return 0.0; } // Check 0 deg.
     }
     else {
     if ((M_PI-angle)>DELTANN) { np++; return 0.0; } // Check 180 deg.
     }
     */
    //double recall = recall3(point1,point2,point3)[0];
    double recall = recallTriple(point1,point2,point3)[0];
    bool ok = recall>THRESHOLD3;
    recall = ok ? recall : -recall;
    if (ok) n1++;
    n3++;
    return recall;
}

// Recall function for 2 points
double* Tracker::recall2(Point &p1, Point &p2)
{
    static XMLP net(NETFILE2);
    static float x[6]={0.,0.,0.,0.,0.,0.};
    
    x[0]     = p1.rz()*0.001;   // rz1 [m]
    x[1]     = p1.phi();        // phi1
    x[2]     = p1.z()*0.001;    // z1 [m]
    x[3]     = p2.rz()*0.001;   // rz2 [m]
    x[4]     = p2.phi();        // phi2
    x[5]     = p2.z()*0.001;    // z2 [m]
    
    return net.Recallstep(x);
}

// Recall function for 2 points
double* Tracker::recallPair(treePoint &p1, treePoint &p2)
{
    static XMLP net(NETFILE1);
    static float x[12]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, y[6]={0.,0.,0.,0.,0.,0.};
    static double null[1]={0.};
    
    if (!getFeatures3(p1, p2, y)) return null;

    x[0]    = p1.rz()*0.001;   // rz1 [m]
    x[1]    = p1.phi();        // phi1
    x[2]    = p1.z()*0.001;    // z1 [m]
    x[3]    = p2.rz()*0.001;   // rz2 [m]
    x[4]    = p2.phi();        // phi2
    x[5]    = p2.z()*0.001;    // z2 [m]
    x[6]    = y[0];            // cell's data of p1
    x[7]    = y[1];            // cell's data of p2
    x[8]    = y[2]*0.001;      // distances from origin r (wdist x,y)
    x[9]    = y[3]*0.001;      // zdist2
    x[10]   = y[4]*0.001;      // wdistr
    x[11]   = y[5]*0.001;      // wdist x,y,z
/*
    x[0]    = y[0];            // cell's data of p1
    x[1]    = y[1];            // cell's data of p2
    x[2]    = y[2]*0.001;      // distances from origin r (wdist x,y)
    x[3]    = y[3]*0.001;      // zdist2
    x[4]    = y[4]*0.001;      // wdistr
    x[5]    = y[5]*0.001;      // wdist x,y,z
*/
    return net.Recallstep(x);
}


// Recall function for 2 points
double* Tracker::recallTriple(treePoint &p1, treePoint &p2, treePoint &p3)
{
    static XMLP net(NETFILE4);
    static float x[13]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, y[4]={0.,0.,0.,0.};
    
    y[3]=scoreTripleLogRadius_and_HitDirPoints(p1,p2,p3,y);
    
    x[0]    = p1.rz()*0.001;   // rz1 [m]
    x[1]    = p1.phi();        // phi1
    x[2]    = p1.z()*0.001;    // z1 [m]
    x[3]    = p2.rz()*0.001;   // rz2 [m]
    x[4]    = p2.phi();        // phi2
    x[5]    = p2.z()*0.001;    // z2 [m]
    x[6]    = p3.rz()*0.001;   // rz3 [m]
    x[7]    = p3.phi();        // phi3
    x[8]    = p3.z()*0.001;    // z3 [m]
    x[9]    = y[0];            //Cell's data of p1
    x[10]   = y[1];            //Cell's data of p2
    x[11]   = y[2];            //Cell's data of p3
    x[12]   = y[3];            //log(ir)

    return net.Recallstep(x);
}


// Recall function for 3 points
double* Tracker::recall3(Point &p1, Point &p2, Point &p3)
{
    static XMLP net(NETFILE3);
    static float x[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
    x[0]     = p1.rz()*0.001;   // rz1 [m]
    x[1]     = p1.phi();        // phi1
    x[2]     = p1.z()*0.001;    // z1 [m]
    x[3]     = p2.rz()*0.001;   // rz2 [m]
    x[4]     = p2.phi();        // phi2
    x[5]     = p2.z()*0.001;    // z2 [m]
    x[6]     = p3.rz()*0.001;   // rz3 [m]
    x[7]     = p3.phi();        // phi3
    x[8]     = p3.z()*0.001;    // z3 [m]
    
    return net.Recallstep(x);
}

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
    const int n=3; // Seeding layer combinations
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
            if (verbose) cout << "Start layer " << l1 << " module " << m1 << endl;
            for (auto edge : edgelist) {
                int nextindex = edge.first;
                if (nextindex<0) break;
                int l2 = nextindex/MODULES;
                int m2 = nextindex%MODULES;
                if (verbose) cout << "-> Layer " << l2 << " module " << m2 << endl;
                for (auto a : module[start]) { // all hits in module
                    
                    if (assignment[a] > 0) continue;
                    
                    // Recalculate the vertex
                    long n = pairs.size();
                    if (n>2) {
                        int a = pairs[n-1].first;
                        int b = pairs[n-1].second;
                        int c = pairs[n-2].first;
                        int d = pairs[n-2].second;
                        vertex = Point::distBetweenLines(points[a], points[b], points[c], points[d]);
                        if (verbose) cout << "vertex(" << a << "," << b << "," << c << "," << d <<"): " << vertex.x() << " " << vertex.y() << " " << vertex.z() << " " << endl;
                    }
                    
                    int twin = points[a].twin();
                    if (twin>0) {
                        pairs.push_back(make_pair(a,twin));
                        tracking.add(a,twin,1.0);
                        assignment[a] = ntrack;
                        assignment[twin] = ntrack;
                        if (verbose) cout << "-> Add twin " << a << " " << twin << endl;
                        continue;
                    }
                    
                    // Generate seeding points
                    auto p = module[nextindex]; // all hits in following modules
                    if (p.size() == 0) continue;
                    vector<pair<int,float> > seed;
                    for (auto b:p)
                    {
                        if (assignment[b] > 0) continue;
                        
                        float d = Point::distance3(vertex,points[a],points[b]); // distance of it from line a-b
                        if (d>DISTANCE) {
                            if (verbose) cout << a << " " << b << " <- distance " << d << endl;
                            nd++; continue;
                        }

                        double recall = checkTracklet(a,b); // Search for hit pairs
                        if (recall < THRESHOLD2) {
                            if (verbose) cout << a << " " << b << ": R2 NOK " << recall << endl;
                            continue;
                        }
                        
                        seed.push_back(make_pair(b,d));
                        
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
                    if (verbose) cout << a << " " << id << ": R2 OK " << recall << " d " << d << endl;
                    pairs.push_back(make_pair(a,id));
                    tracking.add(a,id,d);
                    assignment[a] = ntrack++;
                    assignment[id] = ntrack;
                    
                }
            }
        }
        
    }
    
    // Print the tracking graph
    if (verbose) tracking.print();
    
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
                    if (verbose) cout << a << " " << id << ": R2 OK " << recall << " d " << d << endl;
                    pairs.push_back(make_pair(a,id));
                    assignment[id] = ntrack;
                    assignment[a] = ntrack++;
                    tracking.add(a,id,recall);
                }
            }
        }
    }
    
    if (verbose) {
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
long Tracker::findTriples(int p0, int p1, std::vector<int> &pairs)
{
    triple t;
    t.x = p0;
    t.y = p1;
    
    triples.clear();
    
    for (auto &it : pairs)
    {
        //if (!checkTheta(p1,it)) continue;
        //if (!checkRadius(p1,it)) continue;
        //if (!checkDistance(p1,it)) continue;
        
        double recall = checkTracklet(p0,p1,it);
        if (recall > 0) {
            t.z = it;
            t.r = recall;
            triples.push_back(t);
            if (verbose) cout << t.x << " " << t.y << " " << it << ": R3 OK" << recall << ", ";
        }
        else
            if (verbose) cout << t.x << " " << t.y << " " << it << ": R3 NOK" << recall << ", ";
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
        if (verbose) cout << endl << "findTripes " <<  "{" << index << ","  << l << "," << m << "/" << it.first << "," << it.second << "}" << endl;
        vector<triple> trip;
        addHits(it.first,it.second,index,trip);
        for (auto t: trip) {
            tracking.add(t.y,t.z,t.r);
            //  tracking.add(t.z,t.y,t.r);
            //assignment[t.z] = assignment[t.y]  = assignment[t.x];
        }
        triples.insert(triples.end(),trip.begin(),trip.end()); // append the candidates
    }
    if (verbose) tracking.print();
    return triples.size();
}

double scoreTriple(int ai, int bi, int ci);

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
    if (verbose) cout << "{" << start << ","  << l << "," << m << "} -> ";
    
    long found(0);
    for (auto edge : edgelist) {
        float value = edge.second;
        //if (value<10.0) continue;
        int nextindex = edge.first;
        int nextlayer = nextindex/MODULES;
        int nextmodule = nextindex%MODULES;
        if (nextmodule<0 || nextmodule>3192) break;
        if (verbose) cout << "{" << nextindex << "," << nextlayer << "," << nextmodule << "},";
        auto &seed1 = module[nextindex];
        
        vector<pair<double, int> > v;
        for (auto &it1 : seed1)
        {
            if (verbose) cout << it1 << " ";
            if (assignment[it1] > 0) continue; // Point has benn already used
            //if (!checkTheta(p1,it1)) continue;
            treePoint &a = points[p0];
            treePoint &b = points[p1];
            treePoint &c = points[it1];
            double x = scoreTriplePoints(a,b,c); // Check helix propagation

            double s = scoreTriple(p0,p1,it1); // Check helix propagation
            v.push_back(make_pair(s, it1));
        }

        if (verbose) {
            cout << "addHits " << p0 << " " << p1 <<": " << v.size() << " matches" << endl;
            for (auto it : v) cout << "{" << it.first << "," << it.second << "},";
            cout << endl;
        }
        
        //Take only best ones
        sort(v.begin(), v.end());
        
        v.resize(5); // Take only the best solutions
        
        vector<pair<double, int> > r;
        for (int i = 0; i < v.size(); i++) {
            int p2 = v[i].second;
            double recall = checkTracklet(p0,p1,p2); // Point is a candidate on the next layer
            if (recall>THRESHOLD3) { //Prune most triples
                r.push_back(make_pair(recall,p2));
            }
        }
        
        if (r.size()==0) continue;

        //Take only best one
        sort(r.begin(), r.end(), greater<pair<double,int> >()); // Descending order

        int pbest = r.begin()->second;
        double rbest = r.begin()->first;

        t.z = pbest;
        t.r = rbest;
        if (verbose) cout << "push_back {" << t.x << "," << t.y << "," << t.z << "," << t.r << "},";
        triples.push_back(t);
        found +=  addHits(p1,pbest,nextindex,triples);
    }
    
    return found;
}

//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double Tracker::scoreTriplePoints(treePoint &ai, treePoint &bi, treePoint &ci) {
    Point center;
    double radius;
    Point::circle(ai, bi, ci, center, radius);
    if (radius == 0.0) return 1.e9;
    
    Point cb = ci-bi;
    Point ba = bi-ai;
    double ang_cb = asin(dist(cb.x(), cb.y())*.5/radius)*2;
    double ang_ba = asin(dist(ba.x(), ba.y())*.5/radius)*2;
    if (radius != radius || fabs(radius) > 1e50) {
        ang_cb = dist(cb.x(), cb.y());
        ang_ba = dist(ba.x(), ba.y());
    }
    if (ba.z()*cb.z() < 0) ang_ba *= -1;
    
    //if (dist(cb.x, cb.y)*.5/radius > M_PI/2 || dist(ba.x, ba.y)*.5/radius > M_PI/2) return 1e9;
    //-radius*2e-5+
    double x = ba.z() ? (fabs(cb.z()*ang_ba/ba.z()-ang_cb))*radius : 1e9;
    double y = ang_cb ? (fabs(cb.z()*ang_ba/ang_cb-ba.z())) : 1e9;
    double score = min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
    /*
     cout << endl;
     cout << truth_mom[bi]*part_q[truth_part[bi]] << endl;
     point rr = hits[bi]-center;
     if (cb.x*rr.y-cb.y*rr.x > 0) ang_cb *= -1;
     cout << point(-rr.y, rr.x, cb.z/ang_cb) << endl;*/
    if (score>100.) return 100.0;
    return score;
}

int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target);
double evaluateScore(int ci, point&dp, point&xp, point&bap);
double getDensity3(point&dp, point&xp, double tt, int li);


//How many outliers do we expect to fit better than "ci" in the triple "ai", "bi", "ci"?
double Tracker::scoreTripleDensityPoints(treePoint &a, treePoint &b, treePoint &c) {
    int ai = a.id();
    int bi = b.id();
    int ci = c.id();
    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, metai[ci], d, dp, xp, bap, polar[ci])) return 1e9;
    double s = evaluateScore(ci, dp, xp, bap);
    s = getDensity3(dp, xp, s, metai[ci]);
    return s+1.e-6;
}

