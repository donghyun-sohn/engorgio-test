#include "utils.h"
#include "comp.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{

    void bitonic_swap_unbounded(int stage, int part, int slots, int length, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &compres,
                                std::vector<Ciphertext<lbcrypto::DCRTPoly>> &ct, std::vector<Ciphertext<lbcrypto::DCRTPoly>> &res, lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        int num_ct = length / slots;
        auto cc = ct[0]->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<std::vector<Ciphertext<lbcrypto::DCRTPoly>>> res_matrix(num_ct, std::vector<Ciphertext<lbcrypto::DCRTPoly>>(3));
        std::vector<double> eliminate(length, 0), eliminate_rot(length, 0), test_ct, test_compres, test_vec_lower;
        // 1level
        for (int i = 0; i < length; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }

        std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());
        std::vector<std::vector<double>> eliminate_split(num_ct);
        eliminate_split = split_slots(eliminate, slots);
        for (int i = 0; i < num_ct; i++)
        {
            Plaintext plaintextDec;
            cc->Decrypt(privateKey, ct[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            std::vector<double> temp = extractRealParts(enc);
            test_ct.insert(test_ct.end(), temp.begin(), temp.end());
        }
        for (int i = 0; i < num_ct; i++)
        {
            Plaintext plaintextDec;
            cc->Decrypt(privateKey, compres[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            std::vector<double> temp = extractRealParts(enc);
            test_compres.insert(test_compres.end(), temp.begin(), temp.end());
        }
        std::vector<Plaintext> plain_eliminate;
        for (int i = 0; i < num_ct; i++)
        {
            plain_eliminate.push_back(cc->MakeCKKSPackedPlaintext(eliminate_split[i]));
        }

        // std::cout << "ct_upper  gen_rotate" << std::endl;
        auto ct_upper = gen_rotate(slots, num_ct, length - block / 2, ct, test_ct, privateKey);
        // std::cout << "ct_lower  gen_rotate" << std::endl;
        auto ct_lower = gen_rotate(slots, num_ct, block / 2, ct, test_ct, privateKey);

        // auto compres_neg = cc->EvalNegate(compres);
        // auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

        // auto compres_rot = cc->EvalRotate(compres, -block / 2);
        // auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -block / 2);

        // auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);

        std::vector<Ciphertext<lbcrypto::DCRTPoly>> compres_neg(num_ct), compres_rot, swap_vec_lower(num_ct), swap_vec_upper, swap_vec_diagonal(num_ct);

        for (int i = 0; i < num_ct; i++)
        {
            compres_neg[i] = cc->EvalNegate(compres[i]);
            swap_vec_lower[i] = cc->EvalAdd(compres_neg[i], plain_eliminate[i]);
        }
        for (int i = 0; i < num_ct; i++)
        {
            Plaintext plaintextDec;
            cc->Decrypt(privateKey, swap_vec_lower[i], &plaintextDec);
            std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
            std::vector<double> temp = extractRealParts(enc);
            test_vec_lower.insert(test_vec_lower.end(), temp.begin(), temp.end());
        }
        // std::cout << "compres_rot  gen_rotate" << std::endl;

        compres_rot = gen_rotate(slots, num_ct, length - block / 2, compres, test_compres, privateKey);
        // std::cout << "swap_vec_upper  gen_rotate" << std::endl;

        swap_vec_upper = gen_rotate(slots, num_ct, length - block / 2, swap_vec_lower, test_vec_lower, privateKey);
        // std::cout << "swap_vec_diagonal  push_back" << std::endl;
        for (int i = 0; i < num_ct; i++)
            swap_vec_diagonal[i] = cc->EvalAdd(compres[i], compres_rot[i]);

        // std::cout << " res_martrix push_back" << std::endl;
        for (int i = 0; i < num_ct; i++)
        {
            res_matrix[i][0] = cc->EvalMult(ct_upper[i], swap_vec_upper[i]);
            res_matrix[i][1] = cc->EvalMult(ct[i], swap_vec_diagonal[i]);
            res_matrix[i][2] = cc->EvalMult(ct_lower[i], swap_vec_lower[i]);
            // }
        }
        if (res.size() < num_ct)
        {
            std::cerr << "Error: res vector size is insufficient." << std::endl;
            return;
        }
        // std::cout << " res_martrix EvalAdd" << std::endl;
        for (int i = 0; i < num_ct; i++)
        {
            // std::cout << "EvalAdd1" << std::endl;
            auto matrix_res = cc->EvalAdd(res_matrix[i][0], res_matrix[i][1]);
            // std::cout << "push_back EvalAdd" << std::endl;
            res[i] = cc->EvalAdd(matrix_res, res_matrix[i][2]);
            // std::cout << "push_back finish" << std::endl;
        }
    }

    void bitonic_swap_topk(int slots, int topk_i, int k, Ciphertext<lbcrypto::DCRTPoly> &compres,
                           Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                           lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        std::chrono::system_clock::time_point start, end;
        double time = 0.;
        start = std::chrono::system_clock::now();
        auto cc = ct->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int Block = (1 << int(std::log2(k) + 1)) + (pow(2, topk_i + 1) - 1) * 2 * k;
        int block = 1 << (int(std::log2(k)) + 1);
        int step = (1 << int(std::log2(k))) + (pow(2, topk_i + 1) - 1) * 2 * k;
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
        std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
        // 1level
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
        }
        std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());

        Plaintext plaintextDec;
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << "preprocess time: " << time << "ms" << std::endl;
        // Plaintext plain_eliminate_rot = cc->MakeCKKSPackedPlaintext(eliminate_rot);
        // if (block < slots) {
        // std::cout << " EvalRotate: " << std::endl;
        auto ct_upper = cc->EvalRotate(ct, -step);
        auto ct_lower = cc->EvalRotate(ct, step);

        // ct_upper = cc->EvalMult(ct_upper, plain_eliminate_rot);
        // ct_lower = cc->EvalMult(ct_lower, plain_eliminate);

        std::vector<std::complex<double>> enc;
        start = std::chrono::system_clock::now();
        auto compres_neg = cc->EvalNegate(compres);
        auto compres_rot = cc->EvalRotate(compres, -step);

        auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

        auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -step);

        auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << " EvalMult: " << std::endl;
        start = std::chrono::system_clock::now();
        res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
        res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
        res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
        // }
        res = cc->EvalAdd(res_martrix[0], res_martrix[1]);
        res = cc->EvalAdd(res, res_martrix[2]);
        // res = cc->EvalMult(res, plain_eliminate);
        // cc->Decrypt(privateKey, res, &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 16; i++) {
        //     std::cout << " swap_res: " << enc[i].real() << std::endl;
        // }
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    }

    void bitonic_swap(int stage, int part, int slots, Ciphertext<lbcrypto::DCRTPoly> &compres,
                      Ciphertext<lbcrypto::DCRTPoly> &ct, Ciphertext<lbcrypto::DCRTPoly> &res,
                      lbcrypto::PrivateKey<lbcrypto::DCRTPoly> &privateKey)
    {
        std::chrono::system_clock::time_point start, end;
        double time = 0.;
        start = std::chrono::system_clock::now();
        auto cc = ct->GetCryptoContext();
        // int Block = 1 << (stage + 1);
        int block = 1 << (stage + 1 - part);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
        std::vector<double> eliminate(slots, 0), eliminate_rot(slots, 0);
        // 1level
        for (int i = 0; i < slots; i += block)
        {
            for (int j = 0; j < block / 2; j++)
            {
                eliminate[i + j] = 1;
            }
        }
        std::rotate_copy(eliminate.begin(), eliminate.begin() + block / 2, eliminate.end(), eliminate_rot.begin());

        Plaintext plaintextDec;
        Plaintext plain_eliminate = cc->MakeCKKSPackedPlaintext(eliminate);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << "preprocess time: " << time << "ms" << std::endl;
        // Plaintext plain_eliminate_rot = cc->MakeCKKSPackedPlaintext(eliminate_rot);
        // if (block < slots) {

        auto ct_upper = cc->EvalRotate(ct, -block / 2);
        auto ct_lower = cc->EvalRotate(ct, block / 2);

        // ct_upper = cc->EvalMult(ct_upper, plain_eliminate_rot);
        // ct_lower = cc->EvalMult(ct_lower, plain_eliminate);

        std::vector<std::complex<double>> enc;
        start = std::chrono::system_clock::now();
        auto compres_neg = cc->EvalNegate(compres);
        auto compres_rot = cc->EvalRotate(compres, -block / 2);

        auto swap_vec_lower = cc->EvalAdd(compres_neg, plain_eliminate);

        auto swap_vec_upper = cc->EvalRotate(swap_vec_lower, -block / 2);

        // cc->Decrypt(privateKey, swap_vec_upper, &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 8; i++) {
        //     std::cout << " swap_vec_upper: " << enc[i].real() << std::endl;
        // }
        auto swap_vec_diagonal = cc->EvalAdd(compres, compres_rot);
        // cc->Decrypt(privateKey, swap_vec_diagonal, &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 8; i++) {
        //     std::cout << " swap_vec_diagonal: " << enc[i].real() << std::endl;
        // }
        // inverse(swap_vec_upper);
        // inverse(swap_vec_lower);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << "rotate time: " << time << "ms" << std::endl;

        start = std::chrono::system_clock::now();
        res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_upper));
        res_martrix.push_back(cc->EvalMult(ct, swap_vec_diagonal));
        res_martrix.push_back(cc->EvalMult(ct_lower, swap_vec_lower));
        // cc->Decrypt(privateKey, res_martrix[0], &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 8; i++) {
        //     std::cout << " res_martrix[0]: " << enc[i].real() << std::endl;
        // }
        // cc->Decrypt(privateKey, res_martrix[1], &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 8; i++) {
        //     std::cout << " res_martrix[1]: " << enc[i].real() << std::endl;
        // }
        // cc->Decrypt(privateKey, res_martrix[2], &plaintextDec);
        // plaintextDec->SetLength(slots);
        // enc = plaintextDec->GetCKKSPackedValue();
        // for (int i = 0; i < 8; i++) {
        //     std::cout << " res_martrix[2]: " << enc[i].real() << std::endl;
        // }
        res = cc->EvalAdd(res_martrix[0], res_martrix[1]);
        res = cc->EvalAdd(res, res_martrix[2]);
        end = std::chrono::system_clock::now();
        time = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        // std::cout << "mul and add time: " << time << "ms" << std::endl;
        // }
        // else {
        //     auto ct_partial_upper  = cc->EvalMult(res, plain_eliminate);
        //     auto ct_upper          = cc->EvalRotate(res, -block / 2);
        //     auto swap_vec_upper    = cc->EvalRotate(compres, -block / 2);·
        //     auto swap_vec_diagonal = cc->EvalAdd(compres, swap_vec_upper);
        //     res_martrix.push_back(cc->EvalMult(res, swap_vec_diagonal));
        //     inverse(swap_vec_diagonal);
        //     res_martrix.push_back(cc->EvalMult(ct_upper, swap_vec_diagonal));
        //     res = cc->EvalAddMany(res_martrix);
        // }
    }
}