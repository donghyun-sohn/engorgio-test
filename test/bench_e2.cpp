// E2: Filter+Agg Benchmark — HAMMER vs Engorgio E2E comparison
// Runs Q1, Q6, Q12 with records_num=60175 (SF=0.01)
// Reports both actual wall-clock time and Engorgio's amortized time.

#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <string>
#include <thread>
#include <fstream>
#include <set>
#include "math/chebyshev.h"

using namespace openfhe;
using namespace lbcrypto;

// Pull in the Eval_E2E functions from relational_query_test.cpp.
// We rename its main() so it doesn't conflict with ours.
#define main _relational_query_test_main
#include "relational_query_test.cpp"
#undef main

int main(int argc, char *argv[])
{
    int records_num = 60175;
    if (argc >= 2)
    {
        try { records_num = std::stoi(argv[1]); }
        catch (...) { records_num = 60175; }
    }

    std::cout << "========================================" << std::endl;
    std::cout << "E2: Engorgio Filter+Agg Benchmark" << std::endl;
    std::cout << "records_num = " << records_num << std::endl;
    std::cout << "========================================" << std::endl;

    // --- CPU / thread info ---
    // Count physical cores from /proc/cpuinfo ("core id" unique per "physical id")
    unsigned physical_cores = 0;
    {
        std::ifstream cpuinfo("/proc/cpuinfo");
        std::string line, cpu_model;
        std::set<std::pair<int,int>> seen;
        int phys_id = -1, core_id = -1;
        while (std::getline(cpuinfo, line)) {
            if (cpu_model.empty() && line.find("model name") != std::string::npos)
                cpu_model = line;
            if (line.find("physical id") != std::string::npos)
                phys_id = std::stoi(line.substr(line.find(':') + 1));
            if (line.find("core id") != std::string::npos) {
                core_id = std::stoi(line.substr(line.find(':') + 1));
                if (phys_id >= 0 && core_id >= 0)
                    seen.insert({phys_id, core_id});
                phys_id = core_id = -1;
            }
        }
        physical_cores = seen.empty() ? std::thread::hardware_concurrency() : seen.size();
        std::cout << std::endl << "--- Hardware ---" << std::endl;
        if (!cpu_model.empty()) std::cout << cpu_model << std::endl;
        std::cout << "physical cores = " << physical_cores << std::endl;
        std::cout << "hardware threads = " << std::thread::hardware_concurrency() << std::endl;
    }
#ifdef _OPENMP
    omp_set_num_threads(physical_cores);
    std::cout << "OpenMP threads = " << omp_get_max_threads() << " (bound to physical cores)" << std::endl;
#else
    std::cout << "OpenMP: disabled (WITH_OPENMP=OFF)" << std::endl;
#endif

    // --- FHE param summary (per-query params are set inside Eval_E2E_*) ---
    std::cout << std::endl << "--- FHE Parameters (E2 queries) ---" << std::endl;
#if NATIVEINT == 128
    std::cout << "NATIVEINT = 128" << std::endl;
    std::cout << "scalingModSize = 78, firstModSize = 89" << std::endl;
#else
    std::cout << "NATIVEINT = 64" << std::endl;
    std::cout << "scalingModSize = 50, firstModSize = 60" << std::endl;
#endif
    std::cout << "securityLevel = HEStd_128_classic" << std::endl;
    std::cout << "secretKeyDist = SPARSE_TERNARY" << std::endl;
    std::cout << "polyDegree = 119" << std::endl;
    std::cout << "Q1  multDepth = 27" << std::endl;
    std::cout << "Q6  multDepth = 30" << std::endl;
    std::cout << "Q12 multDepth = 30" << std::endl;
    std::cout << "========================================" << std::endl;

    // --- Q1 ---
    std::cout << std::endl << ">>> Running Q1 <<<" << std::endl;
    double q1_ms = Eval_E2E_q1(records_num);
    std::cout << "Q1_e2e_ms: " << q1_ms << std::endl;

    // --- Q6 ---
    std::cout << std::endl << ">>> Running Q6 <<<" << std::endl;
    double q6_ms = Eval_E2E_q6(records_num);
    std::cout << "Q6_e2e_ms: " << q6_ms << std::endl;

    // --- Q12 ---
    std::cout << std::endl << ">>> Running Q12 <<<" << std::endl;
    double q12_ms = Eval_E2E_q12(records_num);
    std::cout << "Q12_e2e_ms: " << q12_ms << std::endl;

    std::cout << std::endl << "========================================" << std::endl;
    std::cout << "E2 complete." << std::endl;
    return 0;
}
