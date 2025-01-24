#include "ENGORGIO/comp.h"
#include "ENGORGIO/sort.h"
#include "ENGORGIO/quant.h"
#include "ENGORGIO/utils.h"
#include <stdexcept>
using namespace openfhe;
using namespace lbcrypto;

int main(int argc, char *argv[])
{
    bitonic_sort_modular(8, 2, 65536);
    // bitonic_sort_full_table_scan(8, 65536 * 2);
    return 0;
}