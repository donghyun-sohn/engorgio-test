#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"
using namespace openfhe;
using namespace lbcrypto;

void Eval_hybrid_query_1(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    std::cout << " Eval HQ1 with " << num << "records" << std::endl;
    // filter

    double filter_time = 0.;
    // 8*num comp
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    filter_time *= num;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;

    double agg_time = Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    double sort_time = bitonic_sort_query(plain_bits, num);
    std::cout << num << " sort_time time: " << sort_time << "ms" << std::endl;
    double sync_time = sync_test_small(plain_bits, num);
    std::cout << num << " sync_time time: " << sync_time << "ms" << std::endl;
    double total_time = filter_time + agg_time + sort_time + sync_time;
    std::cout << num << "HQ1 total_time time: " << total_time << "ms" << std::endl
              << std::endl
              << std::endl;
}
void Eval_hybrid_query_2(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << " Eval HQ2 with " << num << "records" << std::endl;
    double filter_time = 0.;
    // 8*num comp
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    filter_time *= num;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;

    double agg_time = Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    double topk_time = topk_sort_test(plain_bits, num, 4);
    std::cout << num << " topk_time time: " << topk_time << "ms" << std::endl;
    double sync_time = sync_test_small(plain_bits, num);
    std::cout << num << " sync_time time: " << sync_time << "ms" << std::endl;
    double total_time = filter_time + agg_time + topk_time + sync_time;
    std::cout << num << "HQ2 total_time time: " << total_time << "ms" << std::endl
              << std::endl
              << std::endl;
}
int main(int argc, char *argv[])
{
    std::cout << "--------------------- Eval HQ1 --------------------- " << std::endl;
    for (int i = 64; i < 32769; i *= 2)
    {
        Eval_hybrid_query_1(8, i, 2);
    }
    std::cout << std::endl
              << std::endl;
    std::cout << "--------------------- Eval HQ2 --------------------- " << std::endl;
    for (int i = 64; i < 32769; i *= 2)
    {
        Eval_hybrid_query_2(8, i, 2);
    }
    return 0;
}