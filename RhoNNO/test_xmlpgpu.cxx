// test_xmlpgpu — smoke test: small regression, CPU-train + GPU-Recall parity.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "TXMLPGPU.h"

int main()
{
    const int N = 256;
    std::vector<std::array<double,3>> X(N);
    std::vector<double> Y(N);
    srand(42);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k) X[i][k] = (rand() % 1000) / 1000.0;
        double a = X[i][0], b = X[i][1], c = X[i][2];
        Y[i] = std::max(a, std::max(b, c));
    }

    TXMLPGPU net(3, /*inputRange*/1.0, /*netFile*/"",
                 /*innodes*/3, /*nodes*/{8, 8, 1},
                 /*learnSteps*/{0.05, 0.01, 0.005},
                 /*reluOutput*/false);

    const int epochs = 200;
    for (int e = 0; e < epochs; ++e) {
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

    // CPU forward
    net.SetDisableGPU(true);
    double mse_cpu = 0.0;
    for (int i = 0; i < N; ++i) {
        NNO_INTYPE in[3] = { (NNO_INTYPE)X[i][0], (NNO_INTYPE)X[i][1], (NNO_INTYPE)X[i][2] };
        double* y = net.Recall(in, nullptr);
        double d = Y[i] - y[0];
        mse_cpu += d * d;
    }
    mse_cpu /= N;
    printf("MSE CPU:  %.6g\n", mse_cpu);
    printf("IsGPU() after CPU recall: %d (expected 0)\n", (int)net.IsGPU());

    // GPU forward
    net.SetDisableGPU(false);
    double mse_gpu = 0.0;
    for (int i = 0; i < N; ++i) {
        NNO_INTYPE in[3] = { (NNO_INTYPE)X[i][0], (NNO_INTYPE)X[i][1], (NNO_INTYPE)X[i][2] };
        double* y = net.Recall(in, nullptr);
        double d = Y[i] - y[0];
        mse_gpu += d * d;
    }
    mse_gpu /= N;
    printf("MSE GPU:  %.6g\n", mse_gpu);
    printf("IsGPU():  %d\n", (int)net.IsGPU());

    double md = net.CompareCpuVsGpu((NNO_INTYPE[]){0.3f, 0.7f, 0.5f});
    printf("max |cpu - gpu| on one sample: %.4g\n", md);

    double rel = std::fabs(mse_gpu - mse_cpu) / std::max(1e-12, mse_cpu);
    printf("rel diff: %.3g\n", rel);

    bool ok = TXMLPGPU::HasGpuSupport() && net.IsGPU() && (mse_cpu < 1e-2) && (md < 1e-3);
    printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
