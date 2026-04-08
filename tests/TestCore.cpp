#include "TDataServe.h"
#include "Graph.h"
#include "TMLP.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

int main()
{
    // Test TDataServe init and mean calculation
    TDataServe ds("test","test",2,1);
    float in1[2] = {1.0f, 2.0f};
    float out1[1] = {0.0f};
    float in2[2] = {2.0f, 3.0f};
    float out2[1] = {1.0f};

    ds.Putvec(in1, out1);
    ds.Putvec(in2, out2);
    ds.Init(1);

    assert(ds.GetNumTrnvecs() == 1u);
    assert(ds.GetNumTstvecs() == 1u);
    assert(ds.GetNumvecs() == 2u);

    float* inMean = ds.GetInputMean();
    assert(std::fabs(inMean[0] - 1.5f) < 1e-6f);
    assert(std::fabs(inMean[1] - 2.5f) < 1e-6f);

    float* outMean = ds.GetOutputMean();
    assert(std::fabs(outMean[0] - 0.5f) < 1e-6f);

    // Test VNeuralNet persistence via TMLP
    const std::string netFile = "test_net.txt";
    {
        TMLP net(0.1, 0.1, 2, 2, 1, 1.0, netFile);
        assert(net.GetNetID() == "XMLP");
        assert(net.GetParameters().fInNodes == 2);
        assert(net.GetParameters().fOutNodes == 1);
        net.Save();
    }
    {
        TMLP loaded(netFile);
        assert(loaded.GetFilename() == netFile);
        assert(loaded.GetNetID() == "XMLP");
        assert(loaded.GetParameters().fInNodes == 2);
        assert(loaded.GetParameters().fOutNodes == 1);
    }
    std::remove(netFile.c_str());

    // Test Graph SCC detection
    Graph g(5);
    g.addEdge(1, 0);
    g.addEdge(0, 2);
    g.addEdge(2, 1);
    g.addEdge(0, 3);
    g.addEdge(3, 4);

    int sccs = g.printSCCs();
    assert(sccs == 3);

    return 0;
}
