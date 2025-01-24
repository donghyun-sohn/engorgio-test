#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/utils.h"
#include <chrono>
#include <omp.h>
#include "math/chebyshev.h"

using namespace openfhe;
using namespace lbcrypto;
// 1, 0.5, 0

void Eval_query_simple_q1(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter

    double filter_time = 0.;
    // 8*num comp
    // #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4; i++)
        filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    // std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    filter_time *= num;
    double agg_time = Eval_Agg_4(plain_bits, num);
    // std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;

    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q3(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q3" << std::endl;
    int comp_num = 5;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q4(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    // filter
    std::cout << "tpc-q4" << std::endl;
    int comp_num = 5;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}
// without order by
void Eval_query_simple_q5(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    // filter
    std::cout << "tpc-q5" << std::endl;
    int comp_num = 9;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q10(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q10" << std::endl;
    int comp_num = 6;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q12(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q12" << std::endl;
    int comp_num = 8;
    int agg_num = 2;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 10;
    for (int test = 0; test < num_test; test++)
        // #pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
        // #pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q17(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q17" << std::endl;
    int comp_num = 5;
    int agg_num = 3;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q18(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q18" << std::endl;
    int comp_num = 4;
    int agg_num = 2;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q19(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q19" << std::endl;
    int comp_num = 28;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q21(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q21" << std::endl;
    int comp_num = 11;
    int agg_num = 1;
    std::vector<double> filter_times(comp_num, 0);
    std::vector<double> agg_times(agg_num, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(comp_num)
        for (int i = 0; i < comp_num; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(agg_num)
        for (int i = 0; i < agg_num; i++)
        {
            agg_times[i] += Eval_Agg_1(plain_bits, num);
        }

    double agg_time = *max_element(agg_times.begin(), agg_times.end()) / num_test;
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

void Eval_query_simple_q6(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block)
{
    // filter
    std::cout << "tpc-q6 " << std::endl;
    std::vector<double> filter_times(5, 0);
    // 8*num comp
    int num_test = 1;
    for (int test = 0; test < num_test; test++)
#pragma omp parallel for num_threads(5)
        for (int i = 0; i < 5; i++)
        {
            filter_times[i] += Eval_modular_Strictly_greater_test(plain_bits, block);
            filter_times[i] *= num;
        }

    // std::cout << num << " Eval_modular_Strictly_greater_test finish " << std::endl;
    // for (int i = 0; i < 2; i++)
    //     filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
    // std::cout << num << " Eval_modular_Strictly_equal_test finish " << std::endl;
    double filter_time = *max_element(filter_times.begin(), filter_times.end()) / num_test;
    std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;
    double agg_time = 0.;
    for (int test = 0; test < num_test; test++)
        agg_time += Eval_Agg_1(plain_bits, num);
    std::cout << num << " agg_time time: " << agg_time << "ms" << std::endl;
    agg_time = agg_time / num_test;
    double total_time = filter_time + agg_time;
    std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
}

// void Eval_query_simple_q12(std::uint32_t plain_bits, std::uint32_t num, std::uint32_t block) {
//     //filter

//     double filter_time = 0.;
//     //8*num comp
//     for (int i = 0; i < 4; i++)
//         filter_time += Eval_modular_Strictly_greater_test(plain_bits, block);
//     for (int i = 0; i < 4; i++)
//         filter_time += Eval_modular_Strictly_equal_test(plain_bits, block);
//     filter_time *= num;
//     std::cout << num << " filter_time time: " << filter_time << "ms" << std::endl;

//     double total_time = filter_time;
//     std::cout << num << " total_time time: " << total_time << "ms" << std::endl;
// }

int main(int argc, char *argv[])
{
    for (int i = 512; i < 32769; i *= 2)
    {
        Eval_query_simple_q1(8, i, 2);
    }
    for (int i = 512; i < 32769; i *= 2)
    {
        Eval_query_simple_q12(8, i, 2);
    }
    return 0;
}
