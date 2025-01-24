#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"
using namespace openfhe;
using namespace lbcrypto;

void Eval_vec_query1(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    std::cout << num << "Eval VQ1 with " << num << " records" << std::endl;
    // filter
    int slots = 65536;
    double dis_time = 0.;
    dis_time = Euclid_distance(plain_bits) * block * num / 65536;
    // 8*num comp

    double agg_time = Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << " ms" << std::endl;
    double topk_time = topk_sort_test(plain_bits, num, 4);
    std::cout << num << " topk_time time: " << topk_time << " ms" << std::endl;
    double sync_time = 0.;
    for (int i = 0; i < 4; i++)
        sync_time += sync_test_small(plain_bits, num);
    std::cout << num << " sync time: " << sync_time << " ms" << std::endl;
    double total_time = sync_time + agg_time + topk_time;
    std::cout << num << "VQ1 with" << num << " records total_time time: " << total_time << "ms" << std::endl;
}

void Eval_vec_query2(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    std::cout << num << "Eval VQ2 with " << num << " records" << std::endl;
    // filter
    int slots = 65536;
    double dis_time = 0.;
    dis_time = Euclid_distance(plain_bits) * block * num / 65536;
    // 8*num comp

    double agg_time = Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << " ms" << std::endl;
    double topk_time = topk_sort_test(plain_bits, num, 16);
    std::cout << num << " topk_time time: " << topk_time << " ms" << std::endl;
    double sync_time = 0.;
    for (int i = 0; i < 10; i++)
        sync_time += sync_test_small(plain_bits, num);
    std::cout << num << " sync time: " << sync_time << " ms" << std::endl;
    double total_time = sync_time + agg_time + topk_time;
    std::cout << num << "VQ2 with" << num << " records total_time time: " << total_time << "ms" << std::endl;
}

int main(int argc, char *argv[])
{
    for (int i = 128; i < 1025; i *= 2)
    {
        Eval_vec_query1(8, i, 2);
    }
    for (int i = 128; i < 1025; i *= 2)
    {
        Eval_vec_query1(8, i, 2);
    }
    return 0;
}