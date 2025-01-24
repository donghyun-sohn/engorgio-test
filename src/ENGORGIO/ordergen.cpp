#include "utils.h"
#include "comp.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{

    // length =2/4/8*slots
    void bitonic_comp_unbounded(int stage, int part, int slots, int length, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                                std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res, double precision, std::vector<double> &coefficients,
                                lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct[0]->GetCryptoContext();
        int num_ct = length / slots;
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<double> mask(length, 1), eliminate(length, 0), test_corssrot;
        // if (block < slots)
        for (int i = 0; i < length; i += 2 * Block)
        {
            for (int j = 0; j < Block; j++)
            {
                mask[i + j] = -1;
            }
        }

        // else

        for (int i = 0; i < length; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }

        std::vector<std::vector<double>> mask_split(num_ct), eliminate_split(num_ct);
        mask_split = split_slots(mask, slots);
        eliminate_split = split_slots(eliminate, slots);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_mask(num_ct);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_rot(num_ct);
        std::vector<Plaintext> plain_mask, plain_eliminate;
        for (int i = 0; i < num_ct; i++)
        {
            plain_mask.push_back(cc->MakeCKKSPackedPlaintext(mask_split[i]));
            plain_eliminate.push_back(cc->MakeCKKSPackedPlaintext(eliminate_split[i]));
        }
        for (int i = 0; i < num_ct; i++)
        {
            ct_mask[i] = cc->EvalMult(ct[i], plain_mask[i]);
            ct_mask[i] = cc->EvalMult(ct_mask[i], plain_eliminate[i]);
            Plaintext plaintextDec;
            cc->Decrypt(privateKey, ct_mask[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            std::vector<double> temp = extractRealParts(enc);
            test_corssrot.insert(test_corssrot.end(), temp.begin(), temp.end());
        }
        // std::cout << " gen_rotate " << std::endl;
        ct_rot = gen_rotate(slots, num_ct, block / 2, ct_mask, test_corssrot, privateKey);
        // std::cout << " gen_rotate finish" << std::endl;
        for (int i = 0; i < num_ct; i++)
        {

            ct_rot[i] = cc->EvalMult(ct_rot[i], plain_eliminate[i]);
            Ciphertext<lbcrypto::DCRTPoly> ploy_res;
            // std::cout << "comp_partial" << std::endl;
            comp_partial(ct_mask[i], ct_rot[i], precision, coefficients, res[i], privateKey);
            // std::cout << "comp_partial finish" << std::endl;
        }
    }

    void bitonic_comp_topk(int stage, int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &ct,
                           Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct->GetCryptoContext();
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k;
        int block = 1 << (int(std::log2(k)) + 1);
        int step = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;

        std::vector<double> mask(slots, 0), eliminate(slots, 0);
        // if (block < slots)
        for (int i = 0; i < slots; i += Block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                mask[i + j] = 1;
            }
            // std::cout << "mask i:" << i << "  block / 2: " << block / 2 << std::endl;
            // std::cout << "mask i:" << i << " i + 2 * block - 1: " << i + 2 * block - 1 << std::endl;
            for (int j = step; j < block / 2; j++)
            {
                mask[i + j] = 1;
            }
        }
        // else

        for (int i = 0; i < slots; i += Block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
            for (int j = step; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
            // std::cout << "eliminate i:" << i << "  i + 2 * block - 1: " << i + 2 * block - 1 << std::endl;
        }
        // for (int i = 0; i < 8; i++) {
        //     std::cout << "i:" << i << " mask: " << mask[i] << ", eliminate: " << eliminate[i] << std::endl;
        // }

        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        auto ct_mask = cc->EvalMult(ct, plain_mask);
        // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
        auto ct_rot = cc->EvalRotate(ct_mask, step);

        // std::cout << "\n comp_partial " << std::endl;
        ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

        ct_rot = cc->EvalMult(ct_rot, plain_eliminate);

        // std::cout << " before sign  " << std::endl;
        // output(ct_rot, 8, privateKey);
        // std::cout << std::endl;
        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        // comp_greater_than(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
        // comp_partial(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
        // std::cout << "\n Eval comp_partial: " << block / 2 << std::endl;
        comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
        res = cc->EvalMult(res, plain_eliminate);

        // std::cout << "compres  " << std::endl;
        // output(res, 8, privateKey);
        // std::cout << std::endl;

        // auto ct_rot = cc->EvalRotate(res, -Block / 2);
        // res         = cc->EvalAdd(res, ct_rot);
    }

    void bitonic_comp(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &ct,
                      Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct->GetCryptoContext();
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<double> mask(slots, 1), eliminate(slots, 0);
        // if (block < slots)
        for (int i = 0; i < slots; i += 2 * Block)
        {
            for (int j = 0; j < Block; j++)
            {
                mask[i + j] = -1;
            }
        }
        // else

        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        auto ct_mask = cc->EvalMult(ct, plain_mask);
        // std::cout << "\n EvalRotate: " << block / 2 << std::endl;
        auto ct_rot = cc->EvalRotate(ct_mask, block / 2);

        // std::cout << "\n comp_partial " << std::endl;
        ct_mask = cc->EvalMult(ct_mask, plain_eliminate);

        ct_rot = cc->EvalMult(ct_rot, plain_eliminate);
        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        comp_partial(ct_mask, ct_rot, precision, coefficients, res, privateKey);
    }

    void bitonic_comp_modular(int stage, int part, int slots, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct,
                              Ciphertext<lbcrypto::DCRTPoly> &res, double precision, std::vector<double> &coefficients,
                              lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        auto cc = ct[0]->GetCryptoContext();
        int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        int num_ct = ct.size();
        std::vector<double> mask(slots, 1), eliminate(slots, 0);
        for (int i = 0; i < slots; i += 2 * Block)
        {
            for (int j = 0; j < Block; j++)
            {
                mask[i + j] = -1;
            }
        }
        // else

        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ct_mask(num_ct), ct_rot(num_ct);
        Plaintext plain_mask = cc->MakeCKKSPackedPlaintext(mask);
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        for (int i = 0; i < num_ct; i++)
        {
            ct_mask[i] = cc->EvalMult(ct[i], plain_mask);
            ct_rot[i] = cc->EvalRotate(ct_mask[i], block / 2);
            ct_mask[i] = cc->EvalMult(ct_mask[i], plain_eliminate);
            ct_rot[i] = cc->EvalMult(ct_rot[i], plain_eliminate);
        }
        Ciphertext<lbcrypto::DCRTPoly> ploy_res;
        comp_partial_modular(ct_mask, ct_rot, precision, coefficients, ploy_res, privateKey);
        res = cc->EvalMult(ploy_res, plain_eliminate);
    }
}