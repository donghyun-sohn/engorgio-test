// E3: Sort + Sync Benchmark — HomSort vs HAMMER MPC sort comparison
// Sweeps n = {64, 256, 1024, 4096, 8192}
// For each n:
//   sort_time  = bitonic_sort_modular(8, 2, n)  — sorts n elements (2-block 16-bit key)
//   sync_time  = sync_test_small(8, n)           — syncs 1 payload column
//   total(n)   = sort_time + sync_time
// Fair comparison: HAMMER sorts 2 columns (1 key + 1 payload), so Engorgio
// does 1 sort (covers the key) + 1 sync (covers the payload).

#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <thread>
#include <set>
#include <omp.h>

using namespace openfhe;
using namespace lbcrypto;

int main(int argc, char *argv[])
{
    std::vector<int> ns = {64, 256, 1024, 4096, 8192};

    std::cout << "========================================" << std::endl;
    std::cout << "E3: Engorgio HomSort + Sync Benchmark" << std::endl;
    std::cout << "Sort: bitonic_sort_modular(8, 2, n) — 2-block sort key" << std::endl;
    std::cout << "Sync: sync_test_small(8, n) x1 — 1 payload column" << std::endl;
    std::cout << "n sweep: ";
    for (int n : ns) std::cout << n << " ";
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;

    // --- CPU / thread info ---
    unsigned physical_cores = 0;
    {
        std::ifstream cpuinfo_file("/proc/cpuinfo");
        std::string line, cpu_model;
        std::set<std::pair<int,int>> seen;
        int phys_id = -1, core_id = -1;
        while (std::getline(cpuinfo_file, line)) {
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

    // --- FHE param summary ---
    std::cout << std::endl << "--- FHE Parameters (E3) ---" << std::endl;
#if NATIVEINT == 128
    std::cout << "NATIVEINT = 128" << std::endl;
    std::cout << "Sort: scalingModSize=78, firstModSize=89" << std::endl;
    std::cout << "Sync: scalingModSize=78, firstModSize=89" << std::endl;
#else
    std::cout << "NATIVEINT = 64" << std::endl;
    std::cout << "Sort: scalingModSize=59, firstModSize=60" << std::endl;
    std::cout << "Sync: scalingModSize=55, firstModSize=60" << std::endl;
#endif
    std::cout << "securityLevel = HEStd_128_classic" << std::endl;
    std::cout << "secretKeyDist = SPARSE_TERNARY" << std::endl;
    std::cout << "Sort: polyDegree=247, levelsAfterBoot=25, bootstrap levelBudget={4,4}, numIter=2" << std::endl;
    std::cout << "Sync: bootstrap if n>=4096, levelBudget={4,4}, numIter=1" << std::endl;
    std::cout << "Sort multDepth = 25 + GetBootstrapDepth({4,4})" << std::endl;
    std::cout << "Sync multDepth: n<4096 -> (log2(n)+1)*log2(n)/2+1; else 28+GetBootstrapDepth({4,4})" << std::endl;
    std::cout << "========================================" << std::endl;

    // CSV output
    std::ofstream csv("engorgio_e3_results.csv");
    csv << "n,sort_ms,sync_ms,total_ms" << std::endl;

    // Header for stdout
    std::cout << std::endl;
    std::cout << std::setw(8) << "n"
              << std::setw(16) << "sort_ms"
              << std::setw(16) << "sync_ms"
              << std::setw(16) << "total_ms"
              << std::endl;
    std::cout << std::string(56, '-') << std::endl;

    for (int n : ns)
    {
        std::cout << std::endl << ">>> n=" << n << " <<<" << std::endl;

        // --- Sort ---
        // bitonic_sort_modular returns a vector of cumulative times per stage.
        // The last element is the total raw time (ms) for sorting n elements.
        std::cout << "  Running bitonic_sort_modular(8, 2, " << n << ")..." << std::endl;
        std::vector<double> sort_times = bitonic_sort_modular(8, 2, n);
        double sort_ms = sort_times.back();
        std::cout << "  sort raw time: " << sort_ms << " ms" << std::endl;

        // --- Sync (1 payload column) ---
        // sync_test_small returns amortized time: time_total * num_slots / length.
        // Called once — actual measurement, not estimated.
        std::cout << "  Running sync_test_small(8, " << n << ")..." << std::endl;
        double sync_ms = sync_test_small(8, n);
        std::cout << "  sync time: " << sync_ms << " ms" << std::endl;

        double total_ms = sort_ms + sync_ms;

        // Print result row
        std::cout << std::setw(8) << n
                  << std::setw(16) << std::fixed << std::setprecision(2) << sort_ms
                  << std::setw(16) << sync_ms
                  << std::setw(16) << total_ms
                  << std::endl;

        // Write CSV row
        csv << n << ","
            << std::fixed << std::setprecision(2)
            << sort_ms << ","
            << sync_ms << ","
            << total_ms << std::endl;
    }

    csv.close();
    std::cout << std::endl << "========================================" << std::endl;
    std::cout << "E3 complete. Results saved to engorgio_e3_results.csv" << std::endl;
    return 0;
}
