// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

// Reading trackml event data

#include "Tracker.h"
#include "PolarModule.h"

using namespace std;

// functions to read trackml data (taken from topquark)

void Tracker::readBlacklist(string base_path,int filenum) {
    if (filenum < 1000) return;
    char file[1000];
    sprintf(file, "%s/event%09d-blacklist_particles.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open blacklist\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    long long particle_id;
    while (fscanf(fp, "%lld", &particle_id) == 1) {
        blacklist.insert(particle_id);
    }
    fclose(fp);
}


void Tracker::readTruth(string base_path,int filenum) {
    if (filenum < 1000) return;
    char file[1000];
    sprintf(file, "%s/event%09d-truth.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open truth\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    while (1) {
        int hit_id;
        long long particle_id;
        double tx, ty, tz, tpx, tpy, tpz, weight;
        if (fscanf(fp, "%d,%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &hit_id, &particle_id, &tx, &ty, &tz, &tpx, &tpy, &tpz, &weight) != 9) break;
        if (!particle_id || blacklist.count(particle_id)) continue;
        
        truth_tracks[particle_id].push_back(hit_id);
        truth_pos[hit_id] = point(tx, ty, tz);
        truth_mom[hit_id] = point(tpx, tpy, tpz)*Bfield;
        truth_weight[hit_id] = weight;
        part_weight[particle_id] += weight;
        truth_part[hit_id] = particle_id;
        
        map<long long, int>::iterator it = partIDmap.find(particle_id);
        if( it==partIDmap.end() ){
            cout<<"Particle ID not found in map!!!"<<endl;
            cout<<"ID= "<<hit_id<<" hit "<<hit_id<<" iterator at ID "<<it->first<<endl;
            exit(0);
        }
        
        int newID = it->second;
        if( newID < 0 || newID>= (int)particles.size() ){
            cout<<"Mapped particle ID is wrong!!!"<<endl;
            cout<<"ID= "<<hit_id<<" new ID "<<newID<<endl;
            exit(0);
        }
        
        Particle &p = particles[newID];
        p.hit.push_back(hit_id);
        
    }
    fclose(fp);
    
    cout << truth_tracks.size() << " particles with truth" << endl;
}

void Tracker::readParticles(string base_path,int filenum) {
    char file[1000];
    sprintf(file, "%s/event%09d-particles.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) { printf("couldn't open particles\n"); return; }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    while (1) {
        long long id;
        int type;
        point p, m;
        int hits;
        int q;
        if (fscanf(fp, "%lld,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d", &id, &type, &p.x, &p.y, &p.z, &m.x, &m.y, &m.z, &q, &hits) == -1) break;
        if (blacklist.count(id)) continue;
        start_pos[id] = p;
        start_mom[id] = m*Bfield;
        part_q[id] = -q;
        part_hits[id] = hits;
        
        Particle part;
        partIDmap[ (long long) id ] = (int)particles.size();
        part.id = id;
        part.type = type;
        part.x = p.x;
        part.y = p.y;
        part.z = p.z;
        part.r = 0;
        part.px = m.x;
        part.py = m.y;
        part.pz = m.z;
        part.q = q;
        part.hits = hits;
        particles.push_back(part);
    }
    fclose(fp);
    cout << start_pos.size() << " particles" << endl;
}


bool Tracker::track_cmp(int a, int b) {
    point&ma = truth_mom[a];
    point&mb = truth_mom[b];
    double va = ma*ma, vb = mb*mb;
    if (fabs((va-vb)/va) > 1e-5) return va > vb;
    return (truth_pos[b]-truth_pos[a])*ma > 0;
}


//Sort the hits in each track chronologically
void Tracker::sortTracks() {
    int fails = 0, goods = 0;
    for (auto &p : truth_tracks) {
        vector<int> v;
        for (int hit_id : p.second) {
            v.push_back(hit_id);
        }
        sort(v.begin(), v.end(), track_cmp);
        for (unsigned int i = 0; i < v.size(); i++)
            p.second[i] = v[i];
        int bad = 0;
        for (unsigned int i = 2; i < v.size(); i++) {
            point&a = truth_pos[p.second[i-2]];
            point&b = truth_pos[p.second[i-1]];
            point&c = truth_pos[p.second[i]];
            if ((c.z-b.z)*(b.z-a.z) < 0) {
                fails++;
                bad++;
            }
            else goods++;
        }
    }
}


void Tracker::initOrder() {
    int handpicked[LAYERS] = {};
    int c = 0;
    for (int i = 0; i < 4; i++) handpicked[c++] = 7+i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 7-1-i;
    for (int i = 0; i < 7; i++) handpicked[c++] = 11+i;
    
    for (int i = 0; i < 4; i++) handpicked[c++] = 24+i;
    for (int i = 0; i < 2; i++) handpicked[c++] = 40+i;
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 24-1-i;
        handpicked[c++] = 40-1-i;
    }
    
    for (int i = 0; i < 6; i++) {
        handpicked[c++] = 28+i;
        handpicked[c++] = 42+i;
    }
    for (int i = 0; i < LAYERS; i++) {
        topo[i] = handpicked[i];
        itopo[topo[i]] = i;
        //cout << itopo[i] << "," ;
    }
}


//init layer geometries
void Tracker::initLayers() {
    double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
        { 600, 700, 820, 960, 1100, 1300, 1500}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 7; i++) {
            layer[k*11+i].minr = 30;
            layer[k*11+i].maxr = 176.5;
            layer[k*11+i].avgz = avgz1[k][i];
            layer[k*11+i].type = Disc;
        }
    double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
        { 1220, 1500, 1800, 2150, 2550, 2950}};
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < 6; i++) {
            layer[k*10+i+18].minr = 240;
            layer[k*10+i+18].maxr = 701;
            layer[k*10+i+18].avgz = avgz2[k][i];
            layer[k*10+i+18].type = Disc;
            
            layer[k*8+i+34].minr = 755;
            layer[k*8+i+34].maxr = 1018;
            layer[k*8+i+34].avgz = avgz2[k][i];
            layer[k*8+i+34].type = Disc;
        }
    
    double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
    double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
    double avgr3[2] = {820.2, 1020.2};
    
    for (int i = 0; i < 4; i++) {
        layer[i+7].minz =-491;
        layer[i+7].maxz = 491;
        layer[i+7].avgr = avgr1[i];
        layer[i+7].type = Tube;
    }
    for (int i = 0; i < 4; i++) {
        layer[i+24].minz =-1084;
        layer[i+24].maxz = 1084;
        layer[i+24].avgr = avgr2[i];
        layer[i+24].type = Tube;
    }
    for (int i = 0; i < 2; i++) {
        layer[i+40].minz =-1084;
        layer[i+40].maxz = 1084;
        layer[i+40].avgr = avgr3[i];
        layer[i+40].type = Tube;
    }
    Layer layer2[LAYERS];
    for (int i = 0; i < LAYERS; i++) layer2[i] = layer[i];
    for (int i = 0; i < LAYERS; i++) layer[i] = layer2[topo[i]];
    
    layer[0].var0 = 1e-3;
    layer[1].var0 = 5e-4;
    for (int i = 2; i < 18; i++) layer[i].var0 = 3e-4;
    for (int i = 18; i < 22; i++) layer[i].var0 = 5e-2;
    for (int i = 22; i < LAYERS; i++) layer[i].var0 = i%2 || i == 22 ? 9 : 0.1;
    
    for (int i = 0; i < 4; i++) layer[i].var1 = 0.5;
    for (int i = 4; i < 18; i++) layer[i].var1 = 5;
    for (int i = 18; i < 24; i++) layer[i].var1 = 7;
    for (int i = 24; i < LAYERS; i++) layer[i].var1 = i%2 ? 19 : 11;
}


void Tracker::readHits(string base_path, int filenum) {
    initOrder();
    
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-hits.csv", base_path.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-hits.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open hits\n");
        exit(1);
    }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    
    // For one indexing
    hits.push_back(point(0,0,0));
    polar.push_back(point(0,0,0));
    meta.push_back(point(0,0,0));
    metai.push_back(0);
    
    int layers[9] = {7,4,7,6,4,6,6,2,6};
    int metai_list[9][7];
    int c = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < layers[i]; j++)
            metai_list[i][j] = c++;
    }
    cout << "Detectors: " << c << endl;
    
    for (int i = 0; i < LAYERS; i++) {
        layer[i].minr = layer[i].minz = 1e9;
        layer[i].maxr = layer[i].maxz =-1e9;
    }
    
    while (1) {
        long long hit_id;
        double tx, ty, tz;
        int volume_id, layer_id, module_id;
        if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
        if (hit_id != hits.size()) cout << "Hit id's not as expected" << endl;
        meta.push_back(point(volume_id, layer_id, module_id));
        if (volume_id <= 9)
            volume_id -= 7;
        else if (volume_id <= 14)
            volume_id -= 9;
        else
            volume_id -= 10;
        
        int mi = itopo[metai_list[volume_id][layer_id/2-1]];
        
        hits.push_back(point(tx, ty, tz));
        track_hits[hit_id] = point(tx, ty, tz);
        polar.push_back(point(sqrt(tx*tx+ty*ty), atan2(ty,tx), tz));
        metai.push_back(mi);
        
        double r = sqrt(tx*tx+ty*ty);
        Layer&l = layer[metai_list[volume_id][layer_id/2-1]];
        l.minr = min(l.minr, r);
        l.avgr += r;
        l.maxr = max(l.maxr, r);
        l.minz = min(l.minz, tz);
        l.avgz += tz;
        l.maxz = max(l.maxz, tz);
        l.count++;
        //cerr << tz << ' ' << r << endl;
    }
    fclose(fp);
    cout << hits.size() << " hits" << endl;
    
    initLayers();
    
    for (unsigned int hit_id = 1; hit_id < hits.size(); hit_id++) {
        metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
    }
    
    map<double, double> mir[LAYERS], mar[LAYERS];
    for (unsigned int i = 1; i < hits.size(); i++) {
        int mi = metai[i];
        if (layer[mi].type != Disc) continue;
        double &mir_ = mir[mi][polar[i].z];
        double &mar_ = mar[mi][polar[i].z];
        if (!mir_) mir_ = 1e9;
        mir_ = min(mir_, polar[i].x);
        mar_ = max(mar_, polar[i].x);
    }
    
    map<double, int> zi[LAYERS];
    for (int mi = 0; mi < LAYERS; mi++) {
        if (layer[mi].type != Disc) continue;
        int k = 0;
        for (auto p : mir[mi]) {
            double mir_ = mir[mi][p.first]-1e-5;
            double mar_ = mar[mi][p.first]+1e-5;
            z_minr[mi][k] = mir_;
            z_maxr[mi][k] = mar_;
            disc_z[mi][k] = p.first;
            zi[mi][p.first] = k++;
        }
    }
    
    metaz.resize(hits.size());
    metaz[0] = 0;
    for (unsigned int i = 1; i < hits.size(); i++) {
        int mi = metai[i];
        metaz[i] = meta[i].z;
        if (layer[mi].type == Disc)
            metaz[i] = zi[mi][hits[i].z];
    }
    
}

// Volumes: 7,8,9, 12,13,14, 16,17,18
// Volume 8 is innermost, it has modules 2,4,6,8 concentric cylinders outwards


void Tracker::readDetectors(string base_path) {
    char file[1000];
    sprintf(file, "%s/detectors.csv", base_path.c_str());
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open detectors\n");
        exit(1);
    }
    
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    
    Detector d;
    point c,rx,ry,rz;
    while (fscanf(fp, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &d.volume_id, &d.layer_id, &d.module_id, &c.x, &c.y, &c.z,
                  // &d.rx.x, &d.rx.y, &d.rx.z,
                  // &d.ry.x, &d.ry.y, &d.ry.z,
                  // &d.rz.x, &d.rz.y, &d.rz.z,
                  &rx.x, &ry.x, &rz.x,
                  &rx.y, &ry.y, &rz.y,
                  &rx.z, &ry.z, &rz.z,
                  &d.d, &d.minw, &d.maxw, &d.h, &d.cell_w, &d.cell_h) == 21) {
        d.c =  c;
        d.rx = rx;
        d.ry = ry;
        d.rz = rz;

        //if (d.module_id >= 10000 || d.layer_id >= 1000 || d.volume_id >= 100) cout << "What!?" << endl;
        detectors[d.volume_id*10000000+d.layer_id*10000+d.module_id] = d;
    }
    fclose(fp);
}


void Tracker::readCells(string base_path,int filenum) {
    char file[1000];
    if (filenum >= 1000)
        sprintf(file, "%s/event%09d-cells.csv", base_path.c_str(), filenum);
    else
        sprintf(file, "%s/test/event%09d-cells.csv", base_path.c_str(), filenum);
    FILE*fp = fopen(file, "r");
    if (!fp) {
        printf("couldn't open cells\n");
        exit(1);
    }
    char tmpstr[1000];
    fscanf(fp, "%s", tmpstr);
    cout << tmpstr << endl;
    int hit_id, ch0, ch1;
    double value;
    while (fscanf(fp, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
        hit_cells[hit_id].push_back(make_pair(make_pair(ch0, ch1), value));
        //cout << hit_id << ' ' << ch0 << ' ' << ch1 << ' ' << value << endl;
    }
    fclose(fp);
}



//Calculate direction of each hit with cell's data
void Tracker::initHitDir() {
    for (auto &p : points) {
        int hit_id = p.hitid();
        point m = meta[hit_id];
        Detector &d = detectors[int(m.x)*10000000+int(m.y)*10000+int(m.z)];
        
        //if (!hit_cells[hit_id].size()) cout << "Hit with zero cells" << endl;
        //if (metai[hit_id] < 18) continue;
        
        //Use linear regression for direction
        double mx = 0, my = 0, mw = 0;
        auto &cells = hit_cells[hit_id];
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w;
            double y = c.first.second*d.cell_h;
            mw += w;
            mx += x*w;
            my += y*w;
        }
        mx /= mw;
        my /= mw;
        double mxx = 0, mxy = 0, myy = 0;
        for (auto&c : cells) {
            double w = c.second;
            double x = c.first.first*d.cell_w-mx;
            double y = c.first.second*d.cell_h-my;
            mxx += x*x*w;
            myy += y*y*w;
            mxy += x*y*w;
        }
        //Find eigenvector with minimum eigenvalue
        double a = mxx-myy, b = 2*mxy;
        double x = a+b+sqrt(a*a+b*b);
        double y =-a+b+sqrt(a*a+b*b);
        if (0) {
            double lambda = (mxx+myy+sqrt(a*a+b*b))/2;
            cout << lambda << ' ' << (mxx*x+mxy*y)/x << ' ' << (mxy*x+myy*y)/y << endl;
        }
        
        //Analytical formula for z
        double z = 2*d.d*(fabs(x)/d.cell_w+fabs(y)/d.cell_h+1e-8);
        x *= (cells.size()*1.-1.3);//1.3 != 1 was adjusted empirically
        y *= (cells.size()*1.-1.3);
        point d1(x,y,z), d2(x,y,-z);
        d1 = d.rx*d1.x+d.ry*d1.y+d.rz*d1.z;
        d2 = d.rx*d2.x+d.ry*d2.y+d.rz*d2.z;
        hit_dir[hit_id][0] = normalize(d1);
        hit_dir[hit_id][1] = normalize(d2);
        
        for (int k = 0; k < 2; k++)
            if (hit_dir[hit_id][k]*hits[hit_id] < 0)
                hit_dir[hit_id][k] = hit_dir[hit_id][k]*-1;

        p.setcx(hit_dir[hit_id][1].x);
        p.setcy(hit_dir[hit_id][1].y);
        p.setcz(hit_dir[hit_id][1].z);

        //Write out missalignment to ground truth for plotting
        //if (metai[hit_id] == 0)
        //    cerr << acos(max(fabs(hit_dir[hit_id][0]*truth_mom[hit_id]),
        //                     fabs(hit_dir[hit_id][1]*truth_mom[hit_id]))/dist(truth_mom[hit_id])) << endl;
    }
}


// Logistic regression model

//Return features for logistic regression of a triple. L has cell's data angle errors, return logarithm of inverse radius of helix
double Tracker::scoreTripleLogRadius_and_HitDir(int ai,int bi,int ci,float *L) {

    //Find circle with center p, radius r, going through a, b, and c (in xy plane)
    
    Point a = points[ai], b = points[bi], c = points[ci];
    double ax = a.x()-c.x(), ay = a.y()-c.y(), bx = b.x()-c.x(), by = b.y()-c.y();
    if (ax*by-ay*bx == 0.0) return 0.0; // Make sure to have independent points

    double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
    double idet = .5/(ax*by-ay*bx);
    double x = (aa*by-bb*ay)*idet;
    double y = (ax*bb-bx*aa)*idet;
    //double z = 0.0;
    double r = dist(x,y), ir = 1./r;
    x += c.x();
    y += c.y();

    for (int k = 0; k < 3; k++) {
        int di = k ? k==2 ? ci : bi : ai;
        double rx = points[di].x()-x, ry = points[di].y()-y;
        
        Point ca = points[ci]-points[ai];
        double ang_ca = asin(dist(ca.x(), ca.y())*.5*ir)*2;
        double cross = rx*ca.y()-ry*ca.x();
        
        double dx,dy,dz;
        if (ir) {
            dx =-ry*ang_ca;
            dy = rx*ang_ca;
            dz = ca.z();
            if (cross < 0) dz *= -1;
        } else {
            dx = ca.x();
            dy = ca.y();
            dz = ca.z();
        }

        Point dir(dx,dy,dz);
        Point d(points[ai].cx(),points[ai].cy(),points[ai].cz());
        double dot = Point::dot(dir,d);
        L[k] = acos(fabs(dot));
        if (isnan(L[k])) {
            cout << "NAN " << k << endl;
            return 0.0;
        }

    }

    return log(ir);
}


//Decides which pairs to fit logistic model to
int Tracker::good_pair(int a, int b) {
    if (!samepart(a, b)) return 0;
    point s = start_pos[truth_part[a]];
    if (s.x*s.x+s.y*s.y < 100) return 2; // within circle of 1 cm
    auto &v = truth_tracks[truth_part[a]];
    int ai = metai[a], bi = metai[b];
    if (ai > bi) swap(ai, bi);
    for (int i : v)
        if (metai[i] < ai || (metai[i] > ai && metai[i] < bi)) return 2;
    return 1;
}


//Angle between line through hits ai-bi and cell's data direction of at hit id ai
double Tracker::dir_miss(int ai, int bi) {
    Point d = points[ai]-points[bi];
    Point dir(points[ai].cx(),points[ai].cy(),points[ai].cz());
    double dot = Point::dot(dir,d);
    return acos(fabs(dot)); // Direction between hit and cell
}


//Get some features for the logistic regression model
bool Tracker::getFeatures3(int ai, int bi, float *feature) {
    Point &a = points[ai], &b = points[bi];
    Point d = a-b;
    
    double r1 = dist(a.x(), a.y());
    double r2 = dist(b.x(), b.y());
    double dr = r2-r1;
    
    feature[0] = dir_miss(ai, bi);//Cell's data of ai
    feature[1] = dir_miss(bi, ai);//Cell's data of bi
    //Different distances from origin (like how far does the line through ai-bi pass from the origin)
    feature[2] = wdist(a, d, 0);
    feature[3] = zdist2(a, b);
    feature[4] = wdistr(r1, dr, a.z(), d.z(), 1);
    feature[5] = wdist(a, d, 1);
    
    for (int i=0;i<6;i++) {
        if (isnan(feature[i])) {
            cout << "NAN " << i << endl;
            return false;
        }
        if (feature[i]==0.0) {
            return false;
        }
    }
    
    return true;
}

//Different distances from origin (like how far does the line through ai-bi pass from the origin)
double Tracker::wdistr(double r1, double dr, double az, double dz, double w) {
    double pp = r1*r1+az*az*w;
    double pd = r1*dr+az*dz*w;
    double dd = dr*dr+dz*dz*w;
    double result = pp-pd*pd/dd;
    if (result<0) return 0.0;
    return sqrt(result);
}

double Tracker::wdist(Point &a, Point &d, double w) {
    double pp = a.x()*a.x()+a.y()*a.y()+a.z()*a.z()*w;
    double pd = a.x()*d.x()+a.y()*d.y()+a.z()*d.z()*w;
    double dd = d.x()*d.x()+d.y()*d.y()+d.z()*d.z()*w;
    if (fabs(dd)<1.E-7) dd = 1.E-7;
    double result = pp-pd*pd/dd;
    if (result<0) return 0.0;
    return sqrt(result);
}

double Tracker::zdist(Point &a, Point&b) {
    static Point origin(0.,0.,0.);
    Point p;
    double r;
    Point::circle(origin, a, b, p, r);
    double ang_ab = 2*asin(dist(a.x()-b.x(), a.y()-b.y())*.5/r);
    double ang_a = 2*asin(dist(a.x(), a.y())*.5/r);
    return fabs(a.z()-(b.z()-a.z())*ang_a/ang_ab);
}

double Tracker::zdist2(Point &a, Point &b) {
    static Point origin(0.,0.,0.);
    Point p;
    double r;
    Point::circle(origin, a, b, p, r);
    double ang_ab = 2*asin(dist(a.x()-b.x(), a.y()-b.y())*.5/r);
    double d = dist(a.x(), a.y())*0.5/r;
    if (d>1.0) d = 1.0;
    double ang_a = 2*asin(d);
    double zdist = fabs(b.z()-a.z()-a.z()*ang_ab/ang_a);
    return zdist;
}


// Scores and assignment

void Tracker::scorePairs(vector<pair<int, int> >&pairs) {
    set<long long> found;
    for (auto p : pairs) {
        if (samepart(p.first, p.second)) {
            found.insert(truth_part[p.first]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    cout << "Score: " << score << " from " << pairs.size() << " pairs" << endl;
}


void Tracker::scoreTriples(vector<triple>&triples) {
    set<long long> found;
    for (auto p : triples) {
        if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
            found.insert(truth_part[p.x]);
        }
    }
    double score = 0;
    for (long long p : found) {
        score += part_weight[p];
    }
    cout << "Score: " << score << " from " << triples.size() << " triples" << endl;
}


void Tracker::scorePaths(map<int,vector<int> > &paths) {
    int total_length = 0;
    map<long long, double> score_part;
    for (auto &path : paths) {
        
        total_length += path.second.size();
        map<long long, int> count;
        for (int i : path.second)
            count[truth_part[i]]++;
        long long part=0;
        int ma = 0;
        for (auto p : count) {
            if (p.second >= ma) {
                ma = p.second;
                part = p.first;
            }
        }
        double w = 0;
        for (int j : path.second) {
            if (truth_part[j] == part) w += truth_weight[j];
        }
        score_part[part] = max(score_part[part], w);
    }
    double score = 0;
    for (auto p : score_part) if (p.first) score += p.second;
    cout << "Score: " << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}


void Tracker::scoreAssignment(map<int, int>& assignment) {
    map<int, int> track_length;
    for (auto p : assignment)
        track_length[p.second]++;
    
    double score = 0, score2 = 0;
    for (auto it : truth_tracks) {
        map<int, int> c;
        for (int i : it.second)
            if (assignment.count(i))
                c[assignment[i]]++;
        int pick = -1, maxlen = -1;
        for (auto p : c)
            if (p.second*2 > maxlen &&
                p.second*2 > track_length[p.first]) { pick = p.first; maxlen = p.second*2; }
        if (pick == -1) continue;
        
        for (int i : it.second)
            if (assignment[i] == pick) {
                if (maxlen > it.second.size()) score += truth_weight[i];
                else score2 += truth_weight[i];
            }
    }
    //"Final score: " should be the same as the official score (except for blacklisted electrons)
    cout << "Final score:  " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
    cout << "Short score:  " << score+score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
}

//Investigate some stats on the assignment, to try to see of there are any more points to gather
void Tracker::investigateAssignment(map<int,vector<int> > &solution_paths) {
    cout << solution_paths.size() << endl;
    int good = 0;
    map<long long, int> goodc;
    double stolen = 0, missing = 0, between = 0, duplicate = 0, wrong = 0;
    for (int i = 1; i < solution_paths.size(); i++) {
        map<long long, int> c;
        long long bestp = -1;
        int mi = solution_paths[i].size();
        int missed = 0;
        for (int j : solution_paths[i]) {
            if (j <= 0) {missed++;continue;}
            c[truth_part[j]]++;
            if (c[truth_part[j]]*2 > mi) bestp = truth_part[j];
        }
        if (bestp != -1) {
            goodc[bestp]++;
            if (goodc[bestp] == 1) {
                set<int> found, m;
                for (int j : truth_tracks[bestp]) found.insert(j);
                int mi = 1e9, ma = -1;
                for (int j : solution_paths[i]) {
                    int k = abs(j);
                    if (truth_part[k] == bestp && j <= 0) stolen += truth_weight[k];
                    if (truth_part[k] == bestp) {
                        found.erase(k);
                        mi = min(mi, metai[k]);
                        ma = max(ma, metai[k]);
                        m.insert(metai[k]);
                    }
                }
                for (int j : found) {
                    missing += truth_weight[j];
                    if (metai[j] > mi && metai[j] < ma) between += truth_weight[j];
                    if (m.count(metai[j])) duplicate += truth_weight[j];
                }
                good++;
            }
            //log(-scorepathDensity(solution_paths[i])) << endl;
        } else {
            int bad = 0;
            for (int k = 1; k < solution_paths[i].size()-1; k++) {
                auto&v = solution_paths[i];
                if (v[k-1] > 0 && v[k] < 0 && v[k+1] > 0) bad = 1;
            }
            //cerr << bad*1./solution_paths[i].size() << endl;
            for (int j : solution_paths[i])
                if (j > 0)
                    wrong += truth_weight[j];
        }
    }
    cout << good << endl;
    cout << "Stolen: " << stolen << endl;
    cout << "Missing: " << missing << endl;
    cout << "Between: " << between << endl;
    cout << "Duplicate: " << duplicate << endl;
    cout << "Wrong: " << wrong << endl;
    double outside = 0;
    for (auto&p : part_weight)
        if (!goodc.count(p.first)) outside += p.second;
    cout << "Outside: " << outside << endl;
    
    
    map<int, int> map_assignment;
    for (int i = 1; i < solution_paths.size(); i++)
        for (int j : solution_paths[i])
            if (j > 0)
                map_assignment[j] = i;
    
    double completely = 0;
    for (auto&t : truth_tracks) {
        map<int, int> c;
        for (int i : t.second)
            if (map_assignment[i])
                c[map_assignment[i]]++;
        int ma = 0;
        for (auto&p : c) ma = max(ma, p.second);
        if (ma <= 2) completely += part_weight[t.first];
    }
    cout << "Completely lost: " << completely << endl;
}


//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double Tracker::scoreTriple(int ai, int bi, int ci) {
    point center;
    double radius;
    circle(hits[ai], hits[bi], hits[ci], center, radius);
    
    point cb = hits[ci]-hits[bi];
    point ba = hits[bi]-hits[ai];
    double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
    double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
    if (radius != radius || fabs(radius) > 1e50) {
        ang_cb = dist(cb.x, cb.y);
        ang_ba = dist(ba.x, ba.y);
    }
    if (ba.z*cb.z < 0) ang_ba *= -1;
    
    //if (dist(cb.x, cb.y)*.5/radius > M_PI/2 || dist(ba.x, ba.y)*.5/radius > M_PI/2) return 1e9;
    //-radius*2e-5+
    double x = ba.z ? (fabs(cb.z*ang_ba/ba.z-ang_cb))*radius : 1e9;
    double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e9;
    double score = min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
    /*
     cout << endl;
     cout << truth_mom[bi]*part_q[truth_part[bi]] << endl;
     point rr = hits[bi]-center;
     if (cb.x*rr.y-cb.y*rr.x > 0) ang_cb *= -1;
     cout << point(-rr.y, rr.x, cb.z/ang_cb) << endl;*/
    return score;
}

int adj_thres = 1000; //TODO: tweak
//Decides how curved "approximately straight" helices are supposed to be
const double stretch = 0.02;
//List of coordinates for each layer, used to estimate outlier densities
vector<double> sorted_hits[48];
//Acceleration look-up table for faster lookup in sorted_hits
const int crude_steps = 1<<10;
int crudeIndex[48][crude_steps];
//Scaling to get sorted_hits in range [0,crude_steps)
pair<double, double> crudeIndex_a[48];
//Accumulated polynomial coefficients of second order polynomial for O(1) density look-up. Second dimension (20000) must be bigger than number of hits in the most populated layer
double poly[48][20000][3];
//reused global array for storing matching hits from PolarModule->getNear function
int match[200000];
//Number of directly adjacent hits in layers
int next_layer[48][48];

inline int getIndex(int&li, double x) {
    int ci = x*crudeIndex_a[li].first+crudeIndex_a[li].second;
    ci = min(crude_steps-1, max(0, ci));
    int i = crudeIndex[li][ci];
    //cout << x << ' ' << sorted_hits[li][i] << ' ' << crudeIndex_a[li].first << endl;
    //if (i < 1) i = 1;
    //if (i > sorted_hits[li].size()-2) i = sorted_hits[li].size()-2;
    //static int c = 0, d = 0;
    //while (i+1 < sorted_hits[li].size() && x >= sorted_hits[li][i]) i++, c++;
    //while (i && x < sorted_hits[li][i-1]) i--, c++;
    
    //Might segfault sometimes :)
    while (x >= sorted_hits[li][i]) i++;
    while (x < sorted_hits[li][i-1]) i--;
    
    //d++;
    //if (d%1000000 == 0) cout << c*1./d << ' ' << c << ' ' << d << endl;
    /*{
     int j = upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();
     if (j != i) {
     cout << i << ' ' << j << endl;
     cout << x << endl;
     cout << sorted_hits[li][i-1] << ' ' << sorted_hits[li][i] << ' ' << sorted_hits[li][i+1] << endl;
     cout << sorted_hits[li][j-1] << ' ' << sorted_hits[li][j] << ' ' << sorted_hits[li][j+1] << endl;
     cout << endl;
     
     }
     }*/
    return i;//max(0,min(i,int(sorted_hits[li].size())-1));
}

//Init everything needed for fast density calculations, includes most global variables above
void initDensity3() {
    vector<int>*tube = new vector<int>[48]();
    
    for (int i = 1; i < Tracker::hits.size(); i++) {
        if (!Tracker::assignment[i])
            tube[Tracker::metai[i]].push_back(i);
    }
    
    for (int li = 0; li < 48; li++) {
        sorted_hits[li].clear();
        if (Tracker::layer[li].type == Tracker::Tube) {
            for (int i : tube[li])
                sorted_hits[li].push_back(Tracker::hits[i].z);
        } else {
            for (int i : tube[li])
                sorted_hits[li].push_back(Tracker::polar[i].x);
        }
        sorted_hits[li].push_back(-1e50);
        sorted_hits[li].push_back(1e50);
        sort(sorted_hits[li].begin(), sorted_hits[li].end());
        
        double minx = *next(sorted_hits[li].begin())-1e-8;
        double maxx = *next(sorted_hits[li].rbegin())+1e-8;
        //cout << maxx << ' ' << minx << endl;
        double f = crude_steps/(maxx-minx);
        crudeIndex_a[li] = make_pair(f, -minx*f);
        
        for (int i = 0; i < crude_steps; i++) {
            double x = (i+.5)/f+minx;
            crudeIndex[li][i] = upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();
        }
        double acc[3] = {};
        for (int i = 1; i < sorted_hits[li].size(); i++) {
            for (int j = 0; j < 3; j++) poly[li][i][j] = acc[j];
            double x = sorted_hits[li][i];
            for (int j = 0; j < 3; j++)
                acc[j] += pow(x, j);
        }
    }
    delete[]tube;
}

//Get expected number of hits on layer "li" in area (in polar/cylindrical coordnates) spanned by (p-dp)^2-dot(p-dp, xp)^2 < tt
double getDensity3(point&dp, point&xp, double tt, int li) {
    Layer&l = Tracker::layer[li];
    
    double x0, dx, dy;
    if (l.type == Tracker::Tube) {
        x0 = dp.z;
        dx = xp.z;
        dy = xp.y;
    } else {
        x0 = dp.x;
        dx = xp.x;
        dy = xp.y;
    }
    double b = tt*(1-dy*dy)/(1-dx*dx-dy*dy);
    double a = sqrt(1-dx*dx-dy*dy)/((1-dy*dy)*M_PI*dp.x);
    double rx = sqrt(b);
    
    /*int ai = getIndex(li, x0-rx);
     int bi = getIndex(li, x0);
     int ci = getIndex(li, x0+rx);
     
     double ret =
     ((poly[li][bi][0]-poly[li][ai][0])*(rx-x0)+
     (poly[li][bi][1]-poly[li][ai][1])*(1)+
     (poly[li][ci][0]-poly[li][bi][0])*(rx+x0)+
     (poly[li][ci][1]-poly[li][bi][1])*(-1))*a;
     if (ret < 0) {
     cout << x0-rx << ' ' << x0 << ' ' << x0+rx << endl;
     cout << sorted_hits[li][ai] << ' ' << sorted_hits[li][bi] << ' ' << sorted_hits[li][ci] << endl;
     cout << ret << ' ' << getDensity2(dp, xp, tt, li) << endl;
     cout << ai << ' ' << bi << ' ' << ci << endl;
     cout << poly[li][ai][0] << ' ' << poly[li][ai][1] << endl;
     cout << poly[li][bi][0] << ' ' << poly[li][bi][1] << endl;
     cout << poly[li][ci][0] << ' ' << poly[li][ci][1] << endl;
     cout << endl;
     ret = 1e-8;
     }
     return ret;*/
    
    int ai = getIndex(li, x0-rx);
    int bi = getIndex(li, x0+rx);
    if (bi-ai > 10) {//Approximate integration by 2. order polynomial approximation to half disc
        //cout << ai << ' ' << bi << endl;
        const double A = 21*M_PI/64., B = -15*M_PI/64.;
        double ib = 1./b;
        double c0 = A+B*x0*x0*ib, c1 = -2*B*ib*x0, c2 = B*ib;
        double ret =
        ((poly[li][bi][0]-poly[li][ai][0])*c0+
         (poly[li][bi][1]-poly[li][ai][1])*c1+
         (poly[li][bi][2]-poly[li][ai][2])*c2)*a*rx;
        return max(ret,0.);
    } else { //Exact integration, uses half disc
        double density = 0;
        for(int i = ai; i < bi; i++) {
            double x = sorted_hits[li][i]-x0;
            double h = a*sqrt(b-x*x);/// *it;
            //cout << h << endl;
            density += h;
        }
        return density;
    }
}

//Find density by binary search
//This means we want to find (and return) tt such that getDensity3(dp, xp, tt, li) = target
double findDensity(point &dp, point &xp, double target, int li) {
    double Ad = 0, A = 0, B = 1, Bd;
    while (1) {
        Bd = getDensity3(dp, xp, B, li);
        //cout << B << ' ' << Bd << endl;
        if (B > 1e20) {
            cout << "No density?" << endl;
            //cout << dp << ' ' << xp << ' ' << li << endl;
            exit(0);
        }
        if (Bd > target) break;
        B *= 10;
        if (target/Bd < 1e8) B = max(B, target/Bd);
    }
    double mid = B/2;
    int cc = 0;
    while (1) {
        double density = getDensity3(dp, xp, mid, li);
        if (density > target) B = mid, Bd = density;
        else A = mid, Ad = density;
        
        //cout << A << ' ' << mid << ' ' << B << ' ' << density << endl;
        if ((B-A) < A*1e-3 || (density > target*0.9 && density < target*1.1) || cc >= 100) break;
        mid = max(A*0.9+B*0.1, min(B*0.9+A*0.1, (target-Ad)*(B-A)/(Bd-Ad)+A));
        if (++cc == 100) { //Should never happen
            cout << "Warning: Infinite loop in findDensity" << endl;
            /*cout << dp << endl;
             cout << xp << endl;
             cout << mid << endl;
             cout << target << endl;
             cout << li << endl;
             exit(0);*/
        }
    }
    return mid;
}

//Prepare ellipse equation of collision between line extrapolated through hits with id "ai" and "bi" and layer "li". Return collision coordinate "d", in polar coordinates "dp", ellipse stretching "xp", and direction of hit in polar coordnates "bap". "target" describes the layer, possibly corrected for a single point we are evaluating a helix quadruple
int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target) {
    const double slack = 1.00; //No slack
    
    Layer&l = Tracker::layer[li];
    point&a = Tracker::hits[ai], &b = Tracker::hits[bi];
    
    point ba = b-a;
    if (l.type == Tracker::Tube) {
        double vv = ba.x*ba.x+ba.y*ba.y;
        double pv = ba.x*a.x+ba.y*a.y;
        double pp = a.x*a.x+a.y*a.y;
        double RR = target.x*target.x;
        double sq = pv*pv-vv*(pp-RR);
        if (sq < 0) return -1;
        
        double t = (-pv+sqrt(sq))/vv;
        if (t < 0 || !vv) return -1;
        d.x = a.x+ba.x*t;
        d.y = a.y+ba.y*t;
        d.z = a.z+ba.z*t;
        
        if (d.z < l.minz*slack || d.z > l.maxz*slack) return -1;
        
        dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
        
        xp = point(0, -dp.x*(ba.x*ba.x+ba.y*ba.y), ba.z);
        bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
        
        bap = bap*(1./bap.x);
    } else if (l.type == Tracker::Disc) {
        double t = (target.z-a.z)/ba.z;
        if (t < 0 || !ba.z) return -1;
        d.x = a.x+ba.x*t;
        d.y = a.y+ba.y*t;
        d.z = a.z+ba.z*t;
        
        dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);
        
        if (dp.x < l.minr*(1./slack) || dp.x > l.maxr*slack) return -1;
        
        xp = point(ba.x*d.y-ba.y*d.x, d.x*ba.x+d.y*ba.y, 0);
        bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);
        
        bap = bap*(1./bap.z);
    }
    double xp2 = xp.x*xp.x+xp.y*xp.y+xp.z*xp.z;
    if (xp2)
        xp = xp*sqrt((1-stretch)/xp2);
    return 0;
}

//Default is using average position in the layer
int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap) {
    Layer&l = Tracker::layer[li];
    point target(l.avgr, 0, l.avgz);
    return prepareTripleScore(ai, bi, li, d, dp, xp, bap, target);
}

//Use the prepared "dp", "xp", "bap" and return the area that is closer to the collision line (taking into account xp for elliptic behaviour) compared to the hit with id "ci"
double evaluateScore(int ci, point&dp, point&xp, point&bap) {
    point&r = Tracker::polar[ci];
    point err = r-dp;
    if (err.y > M_PI) err.y -= M_PI*2;
    if (err.y <-M_PI) err.y += M_PI*2;
    err.y *= dp.x;
    
    err = err-bap*(Tracker::layer[Tracker::metai[ci]].type == Tracker::Disc ? err.z : err.x);
    double r2 = err*err-pow(err*xp, 2);
    return r2;
}


//Return this if no points were found, somewhat tunable parameter
const double density_eps = 1e-6;

//map<pair<pair<int, int>, int >, double> scoreTripleDensity_mem;

//How many outliers do we expect to fit better than "ci" in the triple "ai", "bi", "ci"?
double scoreTripleDensity(int ai, int bi, int ci) {
    /*auto memi = make_pair(make_pair(ai, bi), ci);
     double&memo = scoreTripleDensity_mem[memi];
     if (!memo) {
     */
    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, Tracker::metai[ci], d, dp, xp, bap, Tracker::polar[ci])) return 1e9;
    double s = evaluateScore(ci, dp, xp, bap);
    s = getDensity3(dp, xp, s, Tracker::metai[ci]);
    return s+density_eps;
    /*
     memo = s+density_eps;
     }
     return memo;*/
}

//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double scoreTriple(int ai, int bi, int ci) {
    point center;
    double radius;
    circle(Tracker::hits[ai], Tracker::hits[bi], Tracker::hits[ci], center, radius);
    
    point cb = Tracker::hits[ci]-Tracker::hits[bi];
    point ba = Tracker::hits[bi]-Tracker::hits[ai];
    double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
    double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
    if (radius != radius || fabs(radius) > 1e50) {
        ang_cb = dist(cb.x, cb.y);
        ang_ba = dist(ba.x, ba.y);
    }
    if (ba.z*cb.z < 0) ang_ba *= -1;
    
    //if (dist(cb.x, cb.y)*.5/radius > M_PI/2 || dist(ba.x, ba.y)*.5/radius > M_PI/2) return 1e9;
    //-radius*2e-5+
    double x = ba.z ? (fabs(cb.z*ang_ba/ba.z-ang_cb))*radius : 1e9;
    double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e9;
    double score = min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
    /*
     cout << endl;
     cout << truth_mom[bi]*part_q[truth_part[bi]] << endl;
     point rr = hits[bi]-center;
     if (cb.x*rr.y-cb.y*rr.x > 0) ang_cb *= -1;
     cout << point(-rr.y, rr.x, cb.z/ang_cb) << endl;*/
    return score;
}

//Return features for logistic regression of a triple. L has cell's data angle errors, return logarithm of inverse radius of helix
double scoreTripleLogRadius_and_HitDir(int ai, int bi, int ci, double*L) {
    point a = Tracker::hits[ai], b = Tracker::hits[bi], c = Tracker::hits[ci];
    //Find circle with center p, radius r, going through a, b, and c (in xy plane)
    double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
    double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
    double idet = .5/(ax*by-ay*bx);
    point p;
    p.x = (aa*by-bb*ay)*idet;
    p.y = (ax*bb-bx*aa)*idet;
    p.z = 0;
    double r = dist(p.x, p.y), ir = 1./r;
    p.x += c.x;
    p.y += c.y;
    
    for (int k = 0; k < 3; k++) {
        int di = k ? k==2 ? ci : bi : ai;
        double rx = Tracker::hits[di].x-p.x, ry = Tracker::hits[di].y-p.y;
        
        point ca = Tracker::hits[ci]-Tracker::hits[ai];
        double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;
        double cross = rx*ca.y-ry*ca.x;
        
        point dir;
        if (ir) {
            dir.x =-ry*ang_ca;
            dir.y = rx*ang_ca;
            dir.z = ca.z;
            if (cross < 0) dir.z *= -1;
        } else {
            dir = ca;
        }
        L[k] = acos(max(fabs(Tracker::hit_dir[di][0]*dir),
                        fabs(Tracker::hit_dir[di][1]*dir))/dist(dir));
    }
    return log(ir);
}

//1 for accepted, scores the triple based on a logistic regression model
int acceptTriple(const triple&t) {
    double A = log(scoreTriple(t.x, t.y, t.z)+1e-8);
    double B = log(scoreTripleDensity(t.x, t.y, t.z));
    double C = log(scoreTripleDensity(t.z, t.y, t.x));
    double L[3];
    double D = scoreTripleLogRadius_and_HitDir(t.x,t.y,t.z,L);
    double w[121] = {-13.215482291516638, -0.519368174205195, -0.6019168737814719, -0.400773825827796, -3.0689189279504614, -8.21987444638849, -1.773083608093787, -3.7271459966647913, -0.18753136282696767, -0.1700350202416788, -0.13325020734065293, -0.0712787103124509, 2.2365502305889295, -0.38264699613950004, -1.5996361946235698, 0.02602607302842127, -0.04090477074387659, -0.12800207114108786, -2.0616314706262706, 0.9350417490331662, -0.6313327964001432, 0.00830034532077729, -0.1021716887039019, 0.3719980432444666, 0.43818671158350325, 0.0338130884608543, 0.19225191422472998, -0.33226136371623366, -1.0631299954951279, -1.3290857109832128, 8.50340840417112, 4.489339976769724, -3.6866359902477703, -1.530677908510526, -0.3660375991432235, -0.2832850515900752, -0.003067393550643885, -0.06185860378584967, -0.004472073355177509, -0.034047188478205974, 0.056232208303619684, -0.09251101374546467, -0.3186456107148592, -0.011497406815609599, 0.0040898730087192275, 0.04166475101451824, 0.5313081554181062, 0.05691563704023761, 0.004054188315119864, 0.009440976068230187, 0.015452389083207108, 0.02857025202533131, -0.01788978369714811, -0.014820867725080077, -0.0032975179225221054, -0.2739810756530984, -0.209895536224461, -0.05555070596049059, -3.8393234795148445, -0.39189992715019867, 0.5302884318217037, -1.0560724732243318, 0.5808249742500916, 0.2085127159157602, -0.002879796716268462, -0.008289453513497825, -0.013327308424637882, 0.034516052559319284, 0.05612738574267425, -0.04698958101602463, 0.0007407605230615924, -0.015547995524776616, 0.06280040184070336, -0.056422842974113374, -0.02553695075115984, -0.030162351232030156, -0.216209409546151, 0.03852063554031595, -0.0693834129966963, -1.0570960495204662, 0.6811799884934827, 0.3386224510850844, -0.10244400357684635, -0.17437169642288441, 0.527447777429105, -0.0009197072806774356, -0.004512546596535816, -0.026048615023175962, -0.016165328699534447, -0.007957908851240184, -0.01677913671380496, 0.00448514125782629, -0.0164129789525374, -0.04792927265651915, 0.3459064488723725, 0.08305188504334206, -0.4177214300084773, -0.09227473501037928, 0.04508615512899353, -0.03988016215006392, 0.029600325479286028, -0.2533468783999991, -0.1438693183062194, -0.17942937900359165, -1.277174888294048, -0.12050721012197445, -1.306910361564254, -0.056617003726385146, -1.1681337555296898, 0.06259298498866638, -6.501290522349262, -10.841239719611016, 2.156866020752887, 1.3871764445557901, 4.945464722802966, -4.26463890108575, -1.510051189434741, -3.140021739172429, -5.693045331329942, 1.1610084032964856, 2.2204604560570425};
    double x[8] = {1,A,B,C,D,L[0],L[1],L[2]};
    int c = 0;
    double score = 0;
    for (int i = 0; i < 8; i++) {
        for (int j = i; j < 8; j++) {
            double a = x[i]*x[j];
            for (int k = j; k < 8; k++)
                score += a*x[k]*w[c++];
        }
    }
    //cout << score << ' ' << w[120] << endl;
    return score > w[120];
}

//Extend triples based on hits "ai" and "bi" to layer "li", add then to "triples", do this using an approximate straight line going through a and b (though taking into account some slightly curved helices by looking into an elliptic region given by "dp", "xp" and "bap". Again postfix "p" for polar coordinates)
double extendTripleLine(vector<triple>&triples, int ai, int bi, point&a, point&b, int li, PolarModule&mod, int rev = 0) {
    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, li, d, dp, xp, bap)) return 0;
    
    //xp = normalize(xp)*0.999;
    //xp = point(0,0,0.5);
    const double target0 = 0.1, target = 10;//0.5;
    
    //double mid0 = -findDensity(dp, xp, target0, li);
    double mid = findDensity(dp, xp, target, li);
    
    int matches = mod.getNear(dp, xp, bap, mid, match);
    
    int mini;
    double best = target;
    vector<pair<double, int> > v;
    for (int i = 0; i < matches; i++) {
        int ci = match[i];
        //if (ci == ai || ci == bi) continue;
        double s = scoreTriple(ai, bi, ci);//evaluateScore(ci, dp, xp, bap);//
        v.push_back(make_pair(s, ci));
    }
    //Take only best ones
    sort(v.begin(), v.end());
    for (int i = 0; i < v.size(); i++) {
        if (i >= target) break;
        triple t(ai, bi, v[i].second);
        if (rev) swap(t.x,t.z);
        if (acceptTriple(t)) //Prune most triples
            triples.push_back(t);
    }
    return 1;
    //cout << added << endl;
}

//Approximate magnetic field strengh as a function of z coordinate, decays drastically near the ends
double field(double z) {
    z *= 1./2750;
    double z2 = z*z;
    return 1.002-z*3e-2-z2*(0.55-0.3*(1-z2));
}

//Similar to prepareTripleScore, but now we extend the helix that passes through hits with id "ai", "bi", "ci". Assumes li > metai[ci] if sign = 1, and metai[bi] < li < metai[ci] if sign = -1
int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target, double sign = 1) {
    Layer&l = Tracker::layer[li];
    
    point p;
    double r, ir;
    
    point c = Tracker::hits[ci];
    point cb = Tracker::hits[ci]-Tracker::hits[bi];//c-b;
    double ang_cb;
    
    if (0) {
        point a = Tracker::hits[ai], b = Tracker::hits[bi], c = Tracker::hits[ci];
        
        //TODO: test if has bieffects
        if (l.type == Tracker::Disc && (c.z-b.z)*(target.z-c.z) < 0) return -1;
        
        //Find circle with center p, radius r, going through a, b, and c (in xy plane)
        double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
        double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
        double idet = .5/(ax*by-ay*bx);
        point p;
        p.x = (aa*by-bb*ay)*idet;
        p.y = (ax*bb-bx*aa)*idet;
        p.z = 0;
        double r = dist(p.x, p.y), ir = 1./r;
        p.x += c.x;
        p.y += c.y;
        
        ang_cb = asin(dist(cb.x, cb.y)*.5*ir)*2;
    }
    
    if (1) { //Take into account varying magnetic field strength
        point a = Tracker::hits[ai], b = Tracker::hits[bi], c = Tracker::hits[ci];
        if (l.type == Tracker::Disc && (c.z-b.z)*(target.z-c.z)*sign < 0) return -1;
        double B1 = field((a.z+b.z)*.5), B2 = field((b.z+c.z)*.5), B3;
        if (l.type == Tracker::Disc) B3 = field((c.z+target.z)*.5);
        else B3 = field(c.z);
        //B1 = B2 = B3 = 1;
        double ax = b.x-a.x, ay = b.y-a.y, bx = c.x-b.x, by = c.y-b.y;
        double aa = ax*ax+ay*ay, dot = ax*bx+ay*by, cross = ax*by-ay*bx;
        double alpha = B2/(2*B3), beta = (-B1*aa-B2*dot)/(2*cross*B3);
        //alpha *= -1;
        double rx = alpha*bx-beta*by, ry = alpha*by+beta*bx;
        p = point(c.x-rx, c.y-ry, 0);
        r = dist(rx, ry);
        ir = 1./r;
        /*cout << endl;
         //cout << (rx*ax+ry*ay)/(ax*ax+ay*ay) << endl;
         //cout << (rx*bx+ry*by)/(bx*bx+by*by) << endl;
         
         point p(c.x-rx, c.y-ry, 0);
         cout << dist(a.x-p.x, a.y-p.y) << ' ' << dist(b.x-p.x, b.y-p.y) << ' ' << dist(c.x-p.x, c.y-p.y) << endl;
         cout << r << ' ' << dist(rx, ry) << endl;*/
        ang_cb = B3/B2*asin(dist(cb.x, cb.y)*.5*ir*B2/B3)*2;
    }
    
    //exit(0);
    
    const double slack = 1.00;
    
    xp = point(0,0,0); //Circle for now
    
    if (l.type == Tracker::Tube) {
        double RR = target.x*target.x, pp = dist2(p.x, p.y);
        double s = .5+(RR-r*r)/(2*pp);
        double sq = RR/pp-s*s;
        
        /*
         if (truth_part[ai] == 4506417142702082LL) {
         cout << sq << endl;
         cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << endl;
         }
         */
        
        if (sq < 0) return -1;
        
        double t = sqrt(sq);
        if (p.y*c.x-p.x*c.y < 0) t *= -1;
        d.x = p.x*s+p.y*t;
        d.y = p.y*s-p.x*t;
        
        point dc = d-c;
        double A = dist(dc.x, dc.y);
        double B = A*.5*ir;
        double ang_dc = asin(B)*2;
        if (dc.x*cb.x+dc.y*cb.y < 0) ang_cb *= -1;
        
        d.z = c.z+cb.z*ang_dc/ang_cb;
        
        if (!(d.z > l.minz*slack && d.z < l.maxz*slack)) return -1;
        
        point dir;
        double s_ = target.x/pp, t_ = s_*(1-s)/t;
        dir.x = p.x*s_+p.y*t_;
        dir.y = p.y*s_-p.x*t_;
        dir.z = (dc.x*dir.x+dc.y*dir.y)*ir*cb.z/(ang_cb*A*sqrt(1-B*B));
        
        dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
        dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
        //cout << dirp << endl; //dirp.x = l.avgr
        
        dirp = dirp*(1./dirp.x);
        
    } else if (l.type == Tracker::Disc) {
        d.z = target.z;
        double fac = ang_cb/cb.z;
        double ang_dc = (d.z-c.z)*fac;
        
        double sa = sin(ang_dc), ca = cos(ang_dc);
        
        double rx = c.x-p.x, ry = c.y-p.y;
        double cross = rx*cb.y-ry*cb.x;
        if (cross < 0) sa *= -1;
        
        d.x = ca*rx-sa*ry+p.x;
        d.y = sa*rx+ca*ry+p.y;
        
        
        point dir;
        dir.x =-fac*(rx*sa+ry*ca);
        dir.y = fac*(rx*ca-ry*sa);
        dir.z = cross < 0 ? -1 : 1;
        
        
        dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
        
        if (!(dp.x > l.minr*(1./slack) && dp.x < l.maxr*slack)) return -1;
        
        dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
        //cout << dirp << endl; //dirp.x = l.avgr
        
        dirp = dirp*(1./dirp.z);
    }
    return 0;
}

//Default target is average position of layer
int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign = 1) {
    Layer&l = Tracker::layer[li];
    point target(l.avgr, 0, l.avgz);
    return prepareQuadrupleScore(ai, bi, ci, li, d, dp, xp, bap, target, sign);
}

//Extend triples using origin, similar to extendTripleLine. Except that it uses the assumption that the particle start in the origin instead of the straight line assumption. This is not used, as the triples starting from the origin are easy anyway, and therefore captured by extendTripleLine.
double extendTripleOrigin(vector<triple>&triples, int ai, int bi, point&a, point&b, int li, PolarModule&mod, double target = 0.5, int rev = 0) {
    Tracker::hits[0] = Tracker::polar[0] = point(0,0,0);
    
    point d, dp, xp, bap;
    if (prepareQuadrupleScore(0, ai, bi, li, d, dp, xp, bap)) return 0;
    
    //const double target0 = 0.1, target = 0.5;
    
    //double mid0 = findDensity(dp, xp, target0, li);
    double mid = findDensity(dp, xp, target, li);
    
    int matches = mod.getNear(dp, xp, bap, mid, match);
    
    int mini;
    double best = target;
    vector<pair<double, int> > v;
    for (int i = 0; i < matches; i++) {
        int ci = match[i];
        //double s = scoreTriple(ai, bi, ci);
        double s = evaluateScore(ci, dp, xp, bap);
        v.push_back(make_pair(s, ci));
    }
    if (!matches) return 0;
    sort(v.begin(), v.end());
    double thres = v[0].first*4;
    for (int i = 0; i < v.size(); i++) {
        if (i >= 2 || v[i].first > thres) break;// && v[i].first > mid0) break; //i >= 1 and added < 1 gives better than old
        if (rev)
            triples.push_back(triple(v[i].second, bi, ai));
        else
            triples.push_back(triple(ai, bi, v[i].second));
    }
    return 1;
    //cout << added << endl;
}

//Expand all pairs into triples (possibly many triples per pair). method 1 uses origin assumption, method 0 uses straight line assumption.
vector<triple> Tracker::findTriples(vector<pair<int, int> >&pairs, PolarModule* mod, int method = 0, double target = 0.5) {
    vector<triple> triples;
    vector<double> v[48];
    for (auto p : pairs) {
        int ai = p.first, bi = p.second;
        //if (!samepart(ai, bi)) continue;
        point&a = Tracker::hits[ai], &b = Tracker::hits[bi];
        
        int added = 0;
        for (int li = Tracker::metai[bi]+1; added < 100*method+1 && li < 48; li++) {
            if (next_layer[Tracker::metai[bi]][li] < adj_thres) continue;
            if (method == 1 && extendTripleLine(triples, ai, bi, a, b, li, mod[li])) added++;
            if (method == 0 && extendTripleOrigin(triples, ai, bi, a, b, li, mod[li], target)) added++;
        }
        //extendTriple(triples, ai, bi, a, b, metai[ai], mod[metai[ai]]);
        //extendTriple(triples, ai, bi, a, b, metai[bi], mod[metai[bi]]);
        /*for (int li = metai[ai]-1; li >= 0; li--) {
         if (next_layer[li][metai[ai]] < adj_thres) continue;
         if (method == 1 && extendTripleLine(triples, bi, ai, a, b, li, mod[li], 1)) break;
         //if (extendTripleOrigin(triples, bi, ai, a, b, li, mod[li], 1)) break;
         }*/
    }
    /*
     for (int i = 0; i < 48; i++) {
     sort(v[i].begin(), v[i].end());
     double t = v[i].size() ? v[i][v[i].size()*99/100]+0.1 : 0.;
     cout << t << ", ";
     }
     cout << endl;
     */
    
    return triples;
}




