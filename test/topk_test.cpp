#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <stdexcept>
using namespace openfhe;
using namespace lbcrypto;
void Eval_Topk(std::uint32_t plain_bits, std::uint32_t num)
{

    std::cout << num << "Eval Top-1/256" << std::endl;
    topk_sort(plain_bits, num, 1);
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << num << "Eval Top-4/256" << std::endl;
    topk_sort(plain_bits, num, 4);
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << num << "Eval Top-8/256" << std::endl;
    topk_sort(plain_bits, num, 8);
    std::cout << "---------------------------------------------" << std::endl
              << std::endl;
    std::cout << num << "Eval Top-16/256" << std::endl;
    topk_sort(plain_bits, num, 16);
}
int main(int argc, char *argv[])
{
    Eval_Topk(8, 256);
    return 0;
}