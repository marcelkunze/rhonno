// test_xmlpnpu
//
// Smoke test for TXMLPNPU: small synthetic regression, float-CPU vs NPU-Recall.
// Reports MSE of both paths, checks parity, and asserts PASS when errors are
// within tolerance.

#include "TXMLPNPU.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main()
{
    // XOR-ish regression on [0,1]^3 → target y
    const int N = 256;
    std::vector<std::array<double,3>> X(N);
    std::vector<double> Y(N);
    srand(42);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) X[i][k] = (rand() % 1000) / 1000.0;
        double a = X[i][0], b = X[i][1], c = X[i][2];
        Y[i] = std::max(a, std::max(b, c)); // smooth-ish target for ReLU MLP
    }

    // 3 → 8 → 8 → 1, input range 1.0 (inputs already in [0,1])
    TXMLPNPU net(3, /*inputRange*/1.0, /*netFile*/"",
                 /*innodes*/3, /*nodes*/{8, 8, 1},
                 /*learnSteps*/{0.05, 0.01, 0.005},
                 /*reluOutput*/false);

    const int epochs = 200;
    for (int e = 0; e < epochs; ++e) {
        // reshuffle per epoch
        for (int i = N - 1; i > 0; --i) {
            int j = rand() % (i + 1);
            std::swap(X[i], X[j]);
            std::swap(Y[i], Y[j]);
        }
        double S = 0.0;
        for (int i = 0; i < N; ++i) {
            NNO_INTYPE in[3] = { (NNO_INTYPE)X[i][0], (NNO_INTYPE)X[i][1], (NNO_INTYPE)X[i][2] };
            NNO_OUTTYPE target = (NNO_OUTTYPE)Y[i];
            S += net.Train(in, &target);
        }
        if (e % 50 == 49 || e == 0)
            printf("epoch %3d  sum_err=%.6g (mse=%.3g)\n", e + 1, S, S / N);
    }

    // Reference MSE via float forward (CPU, no NPU)
    net.SetDisableNPU(true);
    double mse_cpu = 0.0;
    for (int i = 0; i < N; ++i) {
        NNO_INTYPE in[3] = { (NNO_INTYPE)X[i][0], (NNO_INTYPE)X[i][1], (NNO_INTYPE)X[i][2] };
        double* y = net.Recall(in, nullptr);
        double d = Y[i] - y[0];
        mse_cpu += d * d;
    }
    mse_cpu /= N;
    printf("MSE CPU-EP:  %.6g\n", mse_cpu);
    printf("IsNPU() after CPU recall: %d (expected 0)\n", (int)net.IsNPU());

    // NPU path
    net.SetDisableNPU(false);
    double mse_npu = 0.0;
    for (int i = 0; i < N; ++i) {
        NNO_INTYPE in[3] = { (NNO_INTYPE)X[i][0], (NNO_INTYPE)X[i][1], (NNO_INTYPE)X[i][2] };
        double* y = net.Recall(in, nullptr);
        double d = Y[i] - y[0];
        mse_npu += d * d;
    }
    mse_npu /= N;
    printf("MSE NPU:     %.6g\n", mse_npu);
    printf("IsNPU():     %d\n", (int)net.IsNPU());

    double rel = std::fabs(mse_npu - mse_cpu) / std::max(1e-12, mse_cpu);
    printf("rel |mse_npu - mse_cpu| / mse_cpu: %.3g\n", rel);

    bool ok = (mse_cpu < 1e-2) && (rel < 0.05);
    printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
