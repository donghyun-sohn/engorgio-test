
#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include <chrono>
using namespace openfhe;
using namespace lbcrypto;

void Eval_modular_Strictly_greater_test(std::uint32_t Bg)
{
    for (int i = 1; i < 9; i++)
    {
        Eval_modular_Strictly_greater_than(Bg, i);
    }
}

void Eval_modular_Strictly_equal_test(std::uint32_t Bg)
{
    for (int i = 1; i < 9; i++)
    {
        Eval_modular_Strictly_equal(Bg, i);
    }
}

int main(int argc, char *argv[])
{
    Eval_modular_Strictly_greater_test(8);
    Eval_modular_Strictly_equal_test(8);
    return 0;
}
