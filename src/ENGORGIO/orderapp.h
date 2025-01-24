#include "openfhe.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{
    void bitonic_swap_unbounded(int stage, int part, int slots, int length, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &compres,
                                std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void bitonic_swap_topk(int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &compres,
                           Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

    void bitonic_swap(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &compres,
                      Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey);

}