// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include <vector>
#include <map>
#include "Tracker.h"

using namespace std;

// Data initialization

float* Tracker::_x;
float* Tracker::_y;
float* Tracker::_z;
float* Tracker::_cx;
float* Tracker::_cy;
float* Tracker::_cz;
int* Tracker::_volume;
int* Tracker::_layer;
int* Tracker::_module;
int* Tracker::_hitid;
long long *Tracker::_trackid;
bool Tracker::verbose(false);
int Tracker::maxpairs;
int Tracker::maxparticles;
Point Tracker::vertex;
vector<treePoint> Tracker::points; // hit Points
vector<pair<int, int> > Tracker::pairs; // hit pair combinations
vector<pair<int, int> > Tracker::truepairs; // true hit pair combinations
vector<triple> Tracker::triples; // hit triple combinations
Graph<int> Tracker::paths, Tracker::tracking; // graph to represent particle paths and tracking information
long Tracker::seedsok(0),Tracker::seedstotal(0);
long Tracker::trackletsok(0),Tracker::trackletstotal(0);
unsigned long Tracker::nd(0),Tracker::np(0),Tracker::nt(0);
unsigned long Tracker::n1(0),Tracker::n2(0),Tracker::n3(0),Tracker::n4(0),Tracker::ntwins(0);
vector<point> Tracker::hits; //hit position
vector<Particle> Tracker::particles; //true tracks
map<long long,int> Tracker::partIDmap; // create particle ID->index map
//vector<int> Tracker::knn[MAXDIM][MODULES];
vector<int> Tracker::module[LAYERS*MODULES]; // List of hits in each layer
set<int> Tracker::modules[LAYERS]; // List of modules in each layer
vector<int> Tracker::tube[LAYERS][PHIDIM][THEDIM]; // List of hits in each layer
map<long long, vector<int> > Tracker::truth_tracks; //truth hit ids in each track
map<long long, point> Tracker::track_hits; // Find points in hits
int Tracker::assignment[MAXDIM];
point Tracker::truth_pos[MAXDIM], Tracker::truth_mom[MAXDIM]; //truth position and momentum
double Tracker::truth_weight[MAXDIM]; //weighting of each hit
long long Tracker::truth_part[MAXDIM]; //particle this hit belongs to (particle id)
int Tracker::truth_assignment[MAXDIM]; //particle this hit belongs to (track number)
set<long long> Tracker::blacklist;
map<long long, double> Tracker::part_weight; //weighting of each particle
map<long long, map<int, double> > Tracker::metai_weight; //weighting of each particle hit, also adding duplicates
map<long long, point> Tracker::start_pos; //start position
map<long long, point> Tracker::start_mom; //start momentum
map<long long, int> Tracker::part_q; //start charge
map<long long, int> Tracker::part_hits; // = truth_tracks[particle_id].size()
int Tracker::topo[LAYERS], Tracker::itopo[LAYERS]; //reordering of layers for approximate sorting
vector<point> Tracker::polar; //hit position in polar / cylindrical coordinates
vector<point> Tracker::meta; //volume_id / layer_id / module_id
vector<int> Tracker::metai, Tracker::metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
double Tracker::disc_z[LAYERS][4];
Layer Tracker::layer[LAYERS];
double Tracker::z_minr[LAYERS][4], Tracker::z_maxr[LAYERS][4];
map<int, Detector> Tracker::detectors;
vector<std::pair<pair<int, int>, double> > Tracker::hit_cells[MAXDIM]; //pair<pair<ch0, ch1>, value>
point Tracker::hit_dir[MAXDIM][2]; //The two possible directions of the hit according to the cell's data for each hit

