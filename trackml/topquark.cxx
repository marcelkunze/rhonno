// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

// Reading trackml event data

#include "Tracker.h"

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

