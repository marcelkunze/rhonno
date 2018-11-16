#ifndef _POLARMOD_H_
#define _POLARMOD_H_

//Levels of detail
const int lods = 8;

struct PolarModuleInternal;

struct PolarModule {
    PolarModuleInternal *internal;
    int internals;
    PolarModule() {internal = NULL;internals = 0;}
    PolarModule(int li);
    int getNear(point&dp, point&xp, point&bap, double tt, int*match);
    ~PolarModule() {
        if (internal) {
            //delete[]internal;
            internal = NULL;
        }
    }
};

#endif

