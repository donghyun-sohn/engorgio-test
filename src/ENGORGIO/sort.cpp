#include "utils.h"
#include "orderapp.h"
#include "ordergen.h"
#include <chrono>
#include "math/chebyshev.h"
#include <stdexcept>
namespace openfhe
{
    using namespace lbcrypto;

    std::vector<std::vector<std::vector<double>>> bitonic_test_plain(std::vector<double> vec_unsort,
                                                                     std::vector<std::vector<double>> &plain_sort,
                                                                     std::uint32_t pack_slots)
    {
        // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<> dis(1, 10000);
        std::vector<std::vector<std::vector<double>>> martix;
        int len = vec_unsort.size();
        for (int stage = 0; pack_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                // printf("stage: %d part: %d\n", stage, part);
                std::vector<std::vector<double>> swap_martix; // dimension = 2 or 3, depends on stage and part
                int Block = 1 << (stage + 1);
                int block = 1 << (stage + 1 - part);
                // printf("block: %d\n", block);
                std::vector<double> swap_vector(len, 0);
                std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
                // 1level
                for (int i = 0; i < len; i += block)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        eliminate[i + j] = 1;
                    }
                }

                for (int i = 0; i < len / block; i++)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        int swap_bit;
                        if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                            swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        else
                            swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        swap_vector[i * block + j] = swap_bit;
                        swap_vector[i * block + block / 2 + j] = swap_bit;
                    }
                }

                swap_martix.push_back(swap_vector);
                std::vector<double> swap_vector_upper(len, 0);
                std::vector<double> swap_vector_below(len, 0);
                std::vector<double> swap_vector_tmp(len, 0);

                // if (block < len) {
                for (int i = 0; i < len / block; i++)
                    for (int j = 0; j < block / 2; j++)
                        swap_vector_tmp[i * block + j] = 1;
                for (int i = 0; i < len; i++)
                    swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
                rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                            swap_vector_upper.begin());
                swap_martix.push_back(swap_vector_below);
                swap_martix.push_back(swap_vector_upper);

                std::vector<double> vec_sort_temp(len, 0);

                // if (block < len) {
                std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                            vec_rot_upper.begin());

                for (int i = 0; i < len; i++)
                    vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                    swap_martix[2][i] * vec_rot_upper[i];

                plain_sort.push_back(vec_unsort);
                martix.push_back(swap_martix);
            }
        }
        return martix;
    }

    std::vector<std::vector<std::vector<double>>> bitonic_test_plain_full_table_scan(std::vector<double> vec_unsort,
                                                                                     std::vector<std::vector<double>> &plain_sort,
                                                                                     std::uint32_t pack_slots, int test_stage)
    {
        // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<> dis(1, 10000);
        std::vector<std::vector<std::vector<double>>> martix;
        int len = vec_unsort.size();
        for (int stage = test_stage; pack_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                // printf("stage: %d part: %d\n", stage, part);
                std::vector<std::vector<double>> swap_martix; // dimension = 2 or 3, depends on stage and part
                int Block = 1 << (stage + 1);
                int block = 1 << (stage + 1 - part);
                // printf("block: %d\n", block);
                std::vector<double> swap_vector(len, 0);
                std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
                // 1level
                for (int i = 0; i < len; i += block)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        eliminate[i + j] = 1;
                    }
                }

                for (int i = 0; i < len / block; i++)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        int swap_bit;
                        if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                            swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        else
                            swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        swap_vector[i * block + j] = swap_bit;
                        swap_vector[i * block + block / 2 + j] = swap_bit;
                    }
                }

                swap_martix.push_back(swap_vector);
                std::vector<double> swap_vector_upper(len, 0);
                std::vector<double> swap_vector_below(len, 0);
                std::vector<double> swap_vector_tmp(len, 0);

                // if (block < len) {
                for (int i = 0; i < len / block; i++)
                    for (int j = 0; j < block / 2; j++)
                        swap_vector_tmp[i * block + j] = 1;
                for (int i = 0; i < len; i++)
                    swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
                rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                            swap_vector_upper.begin());
                swap_martix.push_back(swap_vector_below);
                swap_martix.push_back(swap_vector_upper);

                std::vector<double> vec_sort_temp(len, 0);

                // if (block < len) {
                std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                            vec_rot_upper.begin());

                for (int i = 0; i < len; i++)
                    vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                    swap_martix[2][i] * vec_rot_upper[i];

                plain_sort.push_back(vec_unsort);
                martix.push_back(swap_martix);
            }
        }
        return martix;
    }

    void bitonic_sort(int plain_bits, int num_slots)
    {
        std::cout << "\nsort: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 495; //
        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "polyDegree " << polyDegree << std::endl;
        std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
        std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
        std::cout << "Modulus " << cc->GetModulus() << std::endl;
        std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;
        std::cout << "precision: " << precision << std::endl;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input(length, 0);
        // const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        // std::cout << "Moduli chain of pk: " << std::endl;
        // print_moduli_chain(ckks_pk[0]);
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }
        // input[0] = 80;
        // input[1] = 10;
        // input[2] = 40;
        // input[3] = 60;
        // input[4] = 5;
        // input[5] = 81;
        // input[6] = 189;
        // input[7] = 165;
        for (int i = 0; i < 8; i++)
        {
            std::cout << input[i] << std::endl;
        }
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
        //     std::cout << plain_sort_vec[0][i] << ",";
        // }
        std::cout << "ringDim: " << ringDim << std::endl;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
        start = std::chrono::system_clock::now();
        double time_total = 0;
        for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                std::cout << " stage " << stage << " part " << part << std::endl;
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                // std::cout << "bitonic_comp " << std::endl;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                std::cout << "bitonic_swap " << std::endl;
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();
                double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                                                      keyPair.secretKey, 2 * num_slots);
                std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
                          << std::endl;

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    std::cout << "boot " << std::endl;
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                    std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                              << std::endl;
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                    std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
                }

                double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                      keyPair.secretKey, 2 * num_slots);

                std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                          << std::log2(err_afterBoot) << std::endl;
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
                std::cout << "boot time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                          << "ms" << std::endl;
                std::cout << "swap time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                          << "ms" << std::endl;
                std::cout << "comp time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                          << "ms" << std::endl
                          << std::endl;
            }

            std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
            std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                      << "ms" << std::endl
                      << std::endl
                      << std::endl
                      << std::endl;
        }
        end = std::chrono::system_clock::now();
        Plaintext plaintextDec;
        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        std::cout << " sortres " << std::endl;
        for (int i = 0; i < int(num_slots); i++)
        {
            std::cout << enc[i].real() << std::endl;
        }
        std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
        std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
                  << std::endl;
    }

    double bitonic_sort_query(int plain_bits, int num_slots)
    {
        std::cout << "\nsort: " << num_slots << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 495; //
        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input(length, 0);
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }
        for (int i = 0; i < 8; i++)
        {
            std::cout << input[i] << std::endl;
        }
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        // std::cout << "ringDim: " << ringDim << std::endl;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        start = std::chrono::system_clock::now();
        double time_total = 0;
        for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                ;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                }

                // double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                //                                       keyPair.secretKey, 2 * num_slots);

                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
        }
        end = std::chrono::system_clock::now();
        Plaintext plaintextDec;
        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        return time_total / double(length) * num_slots;
    }

    void bitonic_sort_small(int plain_bits, int num_slots)
    {
        std::cout << "\nsort: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 119; //

        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "polyDegree " << polyDegree << std::endl;
        std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
        std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
        std::cout << "Modulus " << cc->GetModulus() << std::endl;
        std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;
        std::cout << "precision: " << precision << std::endl;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input(length, 0);
        // const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        // std::cout << "Moduli chain of pk: " << std::endl;
        // print_moduli_chain(ckks_pk[0]);
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }
        // input[0] = 80;
        // input[1] = 10;
        // input[2] = 40;
        // input[3] = 60;
        // input[4] = 5;
        // input[5] = 81;
        // input[6] = 189;
        // input[7] = 165;
        for (int i = 0; i < 8; i++)
        {
            std::cout << input[i] << std::endl;
        }
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        std::cout << "ringDim: " << ringDim << std::endl;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        std::cout << "number of levels fresh: " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
        start = std::chrono::system_clock::now();
        double time_total = 0;
        for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                std::cout << " stage " << stage << " part " << part << std::endl;
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                std::cout << "bitonic_swap " << std::endl;
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();
                double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                                                      keyPair.secretKey, 2 * num_slots);
                std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
                          << std::endl;

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    std::cout << "boot " << std::endl;
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                    std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort->GetLevel()
                              << std::endl;
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                    std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort->GetLevel() << std::endl;
                }

                double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                      keyPair.secretKey, 2 * num_slots);

                std::cout << "error afterBoot in stage " << stage << " part " << part << " is " << err_afterBoot << " ~ 2^"
                          << std::log2(err_afterBoot) << std::endl;
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
                std::cout << "boot time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                          << "ms" << std::endl;
                std::cout << "swap time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                          << "ms" << std::endl;
                std::cout << "comp time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                          << "ms" << std::endl
                          << std::endl;
            }

            std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
            std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
                      << "ms" << std::endl
                      << std::endl
                      << std::endl
                      << std::endl;
        }
        end = std::chrono::system_clock::now();
        Plaintext plaintextDec;
        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        std::cout << " sortres " << std::endl;
        for (int i = 0; i < int(num_slots); i++)
        {
            std::cout << enc[i].real() << std::endl;
        }
        std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
        std::cout << length << "slots sort amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
                  << std::endl;
    }

    std::vector<std::vector<std::vector<double>>> bitonictopk_test_plain(std::vector<double> vec_unsort,
                                                                         std::vector<std::vector<double>> &plain_sort,
                                                                         std::uint32_t pack_slots, std::uint32_t k)
    {
        // vector<int> vec_unsort = {5, 2, 1, 7, 3, 8, 6, 4};
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<> dis(1, 10000);
        std::vector<std::vector<std::vector<double>>> martix;
        int len = vec_unsort.size();

        for (int stage = 0; int(std::log2(k) + 1) >> stage > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                // printf("stage: %d part: %d\n", stage, part);
                std::vector<std::vector<double>> swap_martix; // dimension = 2 or 3, depends on stage and part
                int Block = 1 << (stage + 1);
                int block = 1 << (stage + 1 - part);
                // printf("block: %d\n", block);
                std::vector<double> swap_vector(len, 0);
                std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
                // 1level
                for (int i = 0; i < len; i += block)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        eliminate[i + j] = 1;
                    }
                }

                for (int i = 0; i < len / block; i++)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        int swap_bit;
                        if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                            swap_bit = vec_unsort[i * block + j] > vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        else
                            swap_bit = vec_unsort[i * block + j] < vec_unsort[i * block + block / 2 + j] ? 1 : 0;
                        swap_vector[i * block + j] = swap_bit;
                        swap_vector[i * block + block / 2 + j] = swap_bit;
                    }
                }

                swap_martix.push_back(swap_vector);
                std::vector<double> swap_vector_upper(len, 0);
                std::vector<double> swap_vector_below(len, 0);
                std::vector<double> swap_vector_tmp(len, 0);

                // if (block < len) {
                for (int i = 0; i < len / block; i++)
                    for (int j = 0; j < block / 2; j++)
                        swap_vector_tmp[i * block + j] = 1;
                for (int i = 0; i < len; i++)
                    swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
                rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - block / 2, swap_vector_below.end(),
                            swap_vector_upper.begin());
                swap_martix.push_back(swap_vector_below);
                swap_martix.push_back(swap_vector_upper);
                std::vector<double> vec_sort_temp(len, 0);

                // if (block < len) {
                std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + block / 2, vec_unsort.end(), vec_rot_below.begin());
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - block / 2, vec_unsort.end(),
                            vec_rot_upper.begin());

                for (int i = 0; i < len; i++)
                    vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                    swap_martix[2][i] * vec_rot_upper[i];
                plain_sort.push_back(vec_unsort);
                martix.push_back(swap_martix);
            }
        }
        std::cout << "partially sort:" << std::endl;
        for (int i = 0; i < len; i++)
        {
            std::cout << vec_unsort[i] << ",";
        }
        std::cout << std::endl;
        for (int i = 0; i < std::log2(pack_slots) - std::log2(k) - 1; i++)
        {
            std::cout << "times:" << i << std::endl;
            std::vector<std::vector<double>> swap_martix; // dimension = 2 or 3, depends on stage and part
            int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
            int block = 1 << (int(std::log2(k)) + 1);
            int step = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
            // printf("block: %d\n", block);
            std::vector<double> swap_vector(len, 0);
            std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
            // 1level
            for (int i = 0; i < len; i += block)
            {
                for (int j = 0; j < block / 2; j++)
                {
                    eliminate[i + j] = 1;
                }
            }

            for (int i = 0; i < len / Block; i++)
            {
                for (int j = 0; j < block / 2; j++)
                {
                    int swap_bit;
                    if (i % 2 == 0)
                        swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + step + j] ? 1 : 0;
                    else
                        swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + step + j] ? 1 : 0;
                    swap_vector[i * Block + j] = swap_bit;
                    swap_vector[i * Block + step + j] = swap_bit;
                    std::cout << "i:" << i * Block + j << ", j" << i * Block + step + j;
                    std::cout << ", comp1:" << vec_unsort[i * Block + j] << ", comp2" << vec_unsort[i * Block + step + j]
                              << std::endl;
                }
            }

            swap_martix.push_back(swap_vector);
            std::vector<double> swap_vector_upper(len, 0);
            std::vector<double> swap_vector_below(len, 0);
            std::vector<double> swap_vector_tmp(len, 0);

            // if (block < len) {
            for (int i = 0; i < len / Block; i++)
                for (int j = 0; j < block / 2; j++)
                    swap_vector_tmp[i * Block + j] = 1;
            for (int i = 0; i < len; i++)
                swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
            rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - step, swap_vector_below.end(),
                        swap_vector_upper.begin());
            swap_martix.push_back(swap_vector_below);
            swap_martix.push_back(swap_vector_upper);
            // }
            // else {
            //     for (int i = 0; i < len; i++)
            //         swap_vector_tmp[i] = 1;
            //     for (int i = 0; i < len; i++)
            //         swap_vector_upper[i] = swap_vector_tmp[i] - swap_martix[0][i];
            //     swap_martix.push_back(swap_vector_upper);
            //     swap_martix.push_back(swap_vector_upper);
            // }
            std::vector<double> vec_sort_temp(len, 0);

            // if (block < len) {
            std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + step, vec_unsort.end(), vec_rot_below.begin());
            rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - step, vec_unsort.end(), vec_rot_upper.begin());

            for (int i = 0; i < len; i++)
                vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                swap_martix[2][i] * vec_rot_upper[i];

            plain_sort.push_back(vec_unsort);
            martix.push_back(swap_martix);

            for (int part = 1; part < std::log2(k); part++)
            {
                // printf("stage: %d part: %d\n", stage, part);
                std::vector<std::vector<double>> swap_martix; // dimension = 2 or 3, depends on stage and part
                int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
                int block = 1 << (int(std::log2(k)) + 1);
                int step = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
                // printf("block: %d\n", block);
                std::vector<double> swap_vector(len, 0);
                std::vector<double> eliminate(len, 0), eliminate_rot(len, 0);
                // 1level
                for (int i = 0; i < len; i += block)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        eliminate[i + j] = 1;
                    }
                }

                for (int i = 0; i < len / Block; i++)
                {
                    for (int j = 0; j < block / 2; j++)
                    {
                        int swap_bit;
                        if ((part == 0 && i % 2 == 0) || (part > 0 && (i * block / Block) % 2 == 0))
                            swap_bit = vec_unsort[i * Block + j] > vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                        else
                            swap_bit = vec_unsort[i * Block + j] < vec_unsort[i * Block + block / 2 + j] ? 1 : 0;
                        swap_vector[i * Block + j] = swap_bit;
                        swap_vector[i * Block + block / 2 + j] = swap_bit;
                    }
                }

                swap_martix.push_back(swap_vector);
                std::vector<double> swap_vector_upper(len, 0);
                std::vector<double> swap_vector_below(len, 0);
                std::vector<double> swap_vector_tmp(len, 0);

                // if (block < len) {
                for (int i = 0; i < len / Block; i++)
                    for (int j = 0; j < block / 2; j++)
                        swap_vector_tmp[i * Block + j] = 1;
                for (int i = 0; i < len; i++)
                    swap_vector_below[i] = swap_vector_tmp[i] * (swap_vector_tmp[i] - swap_martix[0][i]);
                rotate_copy(swap_vector_below.begin(), swap_vector_below.begin() + len - step, swap_vector_below.end(),
                            swap_vector_upper.begin());
                swap_martix.push_back(swap_vector_below);
                swap_martix.push_back(swap_vector_upper);
                // }
                std::vector<double> vec_sort_temp(len, 0);

                // if (block < len) {
                std::vector<double> vec_rot_upper(len, 0), vec_rot_below(len, 0);
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + step, vec_unsort.end(), vec_rot_below.begin());
                rotate_copy(vec_unsort.begin(), vec_unsort.begin() + len - step, vec_unsort.end(), vec_rot_upper.begin());

                for (int i = 0; i < len; i++)
                    vec_unsort[i] = swap_martix[0][i] * vec_unsort[i] + swap_martix[1][i] * vec_rot_below[i] +
                                    swap_martix[2][i] * vec_rot_upper[i];

                plain_sort.push_back(vec_unsort);
                martix.push_back(swap_martix);
            }
            std::cout << "after topk" << std::endl;
            for (int i = 0; i < len; i++)
            {
                std::cout << vec_unsort[i] << ",";
            }
            std::cout << std::endl;
        }

        return martix;
    }

    void topk_sort(int plain_bits, int num_slots, int k)
    {
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 495;
        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input(length, 0);
        const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }
        std::cout << std::endl;
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }
        for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++)
        {
            rotstep.push_back((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k);
            rotstep.push_back(-((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k));
        }
        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        double time_total = 0;
        for (int stage = 0; int(std::log2(k)) * 2 >> stage > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                }

                double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                                                      keyPair.secretKey, 2 * num_slots);
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
        }
        Plaintext plaintextDec_temp;
        for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++)
        {

            cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec_temp);
            std::vector<std::complex<double>> enc = plaintextDec_temp->GetCKKSPackedValue();

            int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
            int block = 1 << (int(std::log2(k)) + 1);
            int step = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            Ciphertext<lbcrypto::DCRTPoly> sort_res;
            start_comp = std::chrono::system_clock::now();
            bitonic_comp_topk(std::log2(k), length, i, k, ciphertext_unsort, comp_res, precision, coefficients,
                              keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            cc->Decrypt(keyPair.secretKey, comp_res, &plaintextDec_temp);
            enc = plaintextDec_temp->GetCKKSPackedValue();
            start_swap = std::chrono::system_clock::now();
            bitonic_swap_topk(length, i, k, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap = std::chrono::system_clock::now();
            cc->Decrypt(keyPair.secretKey, sort_res, &plaintextDec_temp);
            enc = plaintextDec_temp->GetCKKSPackedValue();
            if ((multDepth - sort_res->GetLevel()) < 23)
            { // 23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                start_boot = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot = std::chrono::system_clock::now();
            }
            else
            {
                start_boot = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot = std::chrono::system_clock::now();
            }
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec_temp);
            enc = plaintextDec_temp->GetCKKSPackedValue();

            for (int part = 1; part < std::log2(k) + 1; part++)
            {
                int stage = std::log2(k);
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                }

                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
        }

        Plaintext plaintextDec;
        cc->Decrypt(keyPair.secretKey, ciphertext_unsort, &plaintextDec);
        std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
        std::cout << " topk res " << std::endl;
        for (int i = 0; i < k / 2; i++)
        {
            std::cout << "dec: " << enc[i].real() << std::endl;
            std::cout << "dec: " << enc[num_slots - i].real() << std::endl;
        }
        std::cout << length << "slots topk amortize time: " << time_total / double(length) * num_slots << "ms" << std::endl
                  << std::endl;
    }

    double topk_sort_test(int plain_bits, int num_slots, int k)
    {
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 495;
        std::uint32_t numIterations = 2;
        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input(length, 0);
        const std::vector<DCRTPoly> &ckks_pk = keyPair.publicKey->GetPublicElements();
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }
        std::cout << std::endl;
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }
        for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++)
        {
            rotstep.push_back((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k);
            rotstep.push_back(-((1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k));
        }
        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        double time_total = 0;
        for (int stage = 0; int(std::log2(k)) * 2 >> stage > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                Ciphertext<lbcrypto::DCRTPoly> sort_res;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();
                // double err_afterswap = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], sort_res,
                //                                       keyPair.secretKey, 2 * num_slots);

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                }

                // double err_afterBoot = error_estimate(plain_sort_vec[(1 + stage) * stage / 2 + part], ciphertext_unsort,
                //                                       keyPair.secretKey, 2 * num_slots);
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
        }
        Plaintext plaintextDec_temp;
        for (int i = 0; i < std::log2(num_slots) - std::log2(k) - 1; i++)
        {

            int Block = (1 << int(std::log2(k) + 1)) + (pow(2, i + 1) - 1) * 2 * k;
            int block = 1 << (int(std::log2(k)) + 1);
            int step = (1 << int(std::log2(k))) + (pow(2, i + 1) - 1) * 2 * k;
            Ciphertext<lbcrypto::DCRTPoly> comp_res;
            Ciphertext<lbcrypto::DCRTPoly> sort_res;

            start_comp = std::chrono::system_clock::now();
            bitonic_comp_topk(std::log2(k), length, i, k, ciphertext_unsort, comp_res, precision, coefficients,
                              keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();

            start_swap = std::chrono::system_clock::now();
            bitonic_swap_topk(length, i, k, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
            end_swap = std::chrono::system_clock::now();

            if ((multDepth - sort_res->GetLevel()) < 23)
            { // 23
                sort_res = cc->EvalMult(sort_res, 1 / precision);
                start_boot = std::chrono::system_clock::now();
                ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                end_boot = std::chrono::system_clock::now();
            }
            else
            {
                start_boot = std::chrono::system_clock::now();
                ciphertext_unsort = sort_res;
                end_boot = std::chrono::system_clock::now();
            }
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());

            for (int part = 1; part < std::log2(k) + 1; part++)
            {
                int stage = std::log2(k);
                start_comp = std::chrono::system_clock::now();
                bitonic_comp(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients, keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();
                bitonic_swap(stage, part, length, comp_res, ciphertext_unsort, sort_res, keyPair.secretKey);
                end_swap = std::chrono::system_clock::now();

                if ((multDepth - sort_res->GetLevel()) < 23)
                { // 23
                    sort_res = cc->EvalMult(sort_res, 1 / precision);
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = cc->EvalBootstrap(sort_res, numIterations, 8);
                    ciphertext_unsort = cc->EvalMult(ciphertext_unsort, precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort = sort_res;
                    end_boot = std::chrono::system_clock::now();
                }

                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
        }

        return time_total * num_slots / length;
    }

    void bitonic_sort_full_table_scan(int plain_bits, int length)
    {
        std::cout << "\nfull_table_scan: " << length << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 119; //
        std::uint32_t numIterations = 2;

        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
        std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
        std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int num_slots = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input_big;
        std::vector<std::vector<double>> input_vec;
        int num_ct = length / num_slots;
        int w = plain_bits / 8;
        for (int i = 0; i < length / 2; i++)
        {
            input_big.push_back(i / 512);
        }
        input_big.insert(input_big.end(), input_big.rbegin(), input_big.rend());
        // std::cout << "input_big 1 " << std::endl;
        // for (int i = 0; i < 256; i++)
        // {
        //     std::cout << input_big[i * 256] << ",";
        // }
        // std::cout << std::endl;
        // std::cout << "input_big 2 " << std::endl;
        // for (int i = 0; i < 256; i++)
        // {
        //     std::cout << input_big[i * 256 + length / 2] << ",";
        // }
        // std::cout << std::endl;
        input_vec = split_slots(input_big, num_slots);
        std::vector<std::vector<std::vector<double>>> quant_input(num_ct);
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain_full_table_scan(input_big, plain_sort_vec, length, 16);

        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_unsort;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, num_slots);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        std::vector<Plaintext> plain_big;
        // std::vector<Plaintext> plain_big;
        for (int i = 0; i < num_ct; i++)
        {
            Plaintext plain_0 = cc->MakeCKKSPackedPlaintext(input_vec[i]);
            ciphertext_unsort.push_back(cc->Encrypt(keyPair.publicKey, plain_0));
        }
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        start = std::chrono::system_clock::now();
        double time_total = 0;
        for (int stage = 16; length >> stage + 1 > 0; stage++)
        {
            // std::cout << " stage " << stage << " part " << 0 << std::endl;
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> comp_res_cross(num_ct);
            // std::cout << "bitonic_comp " << std::endl;
            start_comp = std::chrono::system_clock::now();
            // std::cout << " bitonic_comp_unbounded " << std::endl;
            bitonic_comp_unbounded(stage, 0, num_slots, length, ciphertext_unsort, comp_res_cross, precision, coefficients,
                                   keyPair.secretKey);
            end_comp = std::chrono::system_clock::now();
            // std::cout << "bitonic_comp_unbounded finish" << std::endl;
            // std::cout << "comp_res.size: " << comp_res.size() << std::endl;
            start_swap = std::chrono::system_clock::now();
            std::vector<Ciphertext<lbcrypto::DCRTPoly>> sort_res_cross(num_ct);
            // std::cout << "bitonic_swap_unbounded " << std::endl;
            bitonic_swap_unbounded(stage, 0, num_slots, length, comp_res_cross, ciphertext_unsort, sort_res_cross, keyPair.secretKey);
            // std::cout << "bitonic_swap finish" << std::endl;
            end_swap = std::chrono::system_clock::now();
            // double err_afterswap = error_estimate_unbounded(plain_sort_vec[part], sort_res,
            //                                                 keyPair.secretKey);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - sort_res[0]->GetLevel()
            //           << std::endl;
            // std::cout << "error after bitonic_swap is " << err_afterswap << " ~ 2^" << std::log2(err_afterswap)
            //           << std::endl;
            for (int i = 0; i < sort_res_cross.size(); i++)
            {
                if ((multDepth - sort_res_cross[i]->GetLevel()) < 23)
                { // 23

                    sort_res_cross[i] = cc->EvalMult(sort_res_cross[i], 1 / precision);
                    std::cout << "boot " << std::endl;
                    start_boot = std::chrono::system_clock::now();

                    ciphertext_unsort[i] = cc->EvalBootstrap(sort_res_cross[i], numIterations, 8);
                    ciphertext_unsort[i] = cc->EvalMult(ciphertext_unsort[i], precision);
                    end_boot = std::chrono::system_clock::now();
                }
                else
                {
                    start_boot = std::chrono::system_clock::now();
                    ciphertext_unsort[i] = sort_res_cross[i];
                    end_boot = std::chrono::system_clock::now();
                }
            }
            // double err_afterboot = error_estimate_unbounded(plain_sort_vec[part], ciphertext_unsort,
            // keyPair.secretKey);
            // std::cout << "number of levels remaining after bitonic_swap: " << multDepth - ciphertext_unsort->GetLevel()
            //           << std::endl;
            // std::cout << "error after boot is " << err_afterboot << " ~ 2^" << std::log2(err_afterboot)
            //           << std::endl;
            std::cout << "number of levels remaining: " << multDepth - ciphertext_unsort[0]->GetLevel() << std::endl;
            time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                          double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            std::cout << "boot time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                      << "ms" << std::endl;
            std::cout << "swap time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                      << "ms" << std::endl;
            std::cout << "comp time: "
                      << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                      << "ms" << std::endl
                      << std::endl;
            for (int part = 1; stage - part >= 0; part++)
            {
                std::cout << " stage " << stage << " part " << part << std::endl;
                for (int ct_i = 0; ct_i < num_ct; ct_i++)
                {

                    Ciphertext<lbcrypto::DCRTPoly> comp_res;
                    Ciphertext<lbcrypto::DCRTPoly> sort_res;
                    // std::cout << "bitonic_comp " << std::endl;
                    start_comp = std::chrono::system_clock::now();
                    // std::cout << " stage " << stage << " part " << part << std::endl;
                    bitonic_comp(stage - 1, part - 1, num_slots, ciphertext_unsort[ct_i], comp_res, precision, coefficients, keyPair.secretKey);
                    end_comp = std::chrono::system_clock::now();
                    // std::cout << "bitonic_swap " << std::endl;
                    start_swap = std::chrono::system_clock::now();
                    bitonic_swap(stage - 1, part - 1, num_slots, comp_res, ciphertext_unsort[ct_i], sort_res, keyPair.secretKey);
                    end_swap = std::chrono::system_clock::now();

                    if ((multDepth - sort_res->GetLevel()) < 23)
                    { // 23
                        sort_res = cc->EvalMult(sort_res, 1 / precision);
                        std::cout << "boot " << std::endl;
                        start_boot = std::chrono::system_clock::now();
                        ciphertext_unsort[ct_i] = cc->EvalBootstrap(sort_res, numIterations, 8);
                        ciphertext_unsort[ct_i] = cc->EvalMult(ciphertext_unsort[ct_i], precision);
                        end_boot = std::chrono::system_clock::now();
                        std::cout << "number of levels remaining  after boot: " << multDepth - ciphertext_unsort[ct_i]->GetLevel()
                                  << std::endl;
                    }
                    else
                    {
                        start_boot = std::chrono::system_clock::now();
                        ciphertext_unsort[ct_i] = sort_res;
                        end_boot = std::chrono::system_clock::now();
                        std::cout << "number of levels remaining : " << multDepth - ciphertext_unsort[ct_i]->GetLevel() << std::endl;
                    }

                    time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                                  double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                                  double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
                }
                std::cout << "boot time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count())
                          << "ms" << std::endl;
                std::cout << "swap time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count())
                          << "ms" << std::endl;
                std::cout << "comp time: "
                          << double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count())
                          << "ms" << std::endl
                          << std::endl;
            }
            std::vector<double> plain_sort_res;
            for (int i = 0; i < ciphertext_unsort.size(); i++)
            {
                Plaintext plaintextDec;
                cc->Decrypt(keyPair.secretKey, ciphertext_unsort[i], &plaintextDec);
                std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
                std::vector<double> dec_real = extractRealParts(enc);
                plain_sort_res.insert(plain_sort_res.end(), dec_real.begin(), dec_real.end());
            }
            std::cout << " sortres " << std::endl;
            for (int i = 0; i < pow(2, stage); i++)
            {
                std::cout << "sort res: " << plain_sort_res[i] << std::endl;
            }
            // std::cout << "stage " << stage << " sort total time: " << time_total << "ms" << std::endl;
            // std::cout << "stage " << stage << " sort amortize time: " << time_total / double(length) * (1 << (stage + 1))
            //           << "ms" << std::endl
            //           << std::endl
            //           << std::endl
            //           << std::endl;
        }
        end = std::chrono::system_clock::now();

        // std::cout << length << "slots sort total time: " << time_total << "ms" << std::endl;
        std::cout << length << "slots sort amortize time: " << time_total << "ms" << std::endl
                  << std::endl;
    }

    void bitonic_sort_modular(int plain_bits, int blocks, int num_slots)
    {
        std::cout << "\nmodular sort: " << plain_bits << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 59;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);
        std::vector<std::uint32_t> levelBudget = {4, 4};
        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        std::uint32_t polyDegree = 247; //

        std::uint32_t numIterations = 2;

        std::uint32_t levelsAvailableAfterBootstrap = 25;
        usint multDepth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        std::cout << "CyclotomicOrder " << cc->GetCyclotomicOrder() << std::endl;
        std::cout << "RingDimension " << cc->GetRingDimension() << std::endl;
        std::cout << "PlaintextModulus " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<double> input_big(length, 0);
        std::vector<std::vector<double>> input_vec;
        for (int j = 0; j < blocks; j++)
        {
            std::vector<double> input_0(length, 0);
            for (int i = 0; i < length; i++)
            {
                input_0[i] = message(engine);
            }
            input_vec.push_back(input_0);
        }

        for (int i = 0; i < length; i++)
        {
            for (int j = 0; j < blocks; j++)
            {
                input_big[i] += double(int(input_vec[j][i]) << (plain_bits * j));
            }
        }
        for (int i = 0; i < 16; i++)
        {
            std::cout << input_big[i] << std::endl;
        }
        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input_big, plain_sort_vec, length);

        std::vector<Ciphertext<lbcrypto::DCRTPoly>> ciphertext_unsort;
        // for (int i = 0; i < plain_sort_vec[0].size(); i++) {
        //     std::cout << plain_sort_vec[0][i] << ",";
        // }
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        cc->EvalBootstrapSetup(levelBudget);
        cc->EvalBootstrapKeyGen(keyPair.secretKey, length);

        std::vector<double> coefficients = EvalChebyshevCoefficients([](double x) -> double
                                                                     { return ((x > 0) ? 1 : 0); },
                                                                     lowerBound, upperBound, polyDegree);
        std::vector<Plaintext> plain_big;
        for (int i = 0; i < blocks; i++)
        {
            Plaintext plain_0 = cc->MakeCKKSPackedPlaintext(input_vec[i]);
            auto ciphertext_unsort_0 = cc->Encrypt(keyPair.publicKey, plain_0);
            ciphertext_unsort.push_back(ciphertext_unsort_0);
        }
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        double time_total = 0;
        for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                Ciphertext<lbcrypto::DCRTPoly> comp_res;
                std::vector<Ciphertext<lbcrypto::DCRTPoly>> sort_res_vec;
                // std::cout << "bitonic_comp " << std::endl;
                start_comp = std::chrono::system_clock::now();
                bitonic_comp_modular(stage, part, length, ciphertext_unsort, comp_res, precision, coefficients,
                                     keyPair.secretKey);
                end_comp = std::chrono::system_clock::now();
                start_swap = std::chrono::system_clock::now();

                for (int i = 0; i < ciphertext_unsort.size(); i++)
                {
                    Ciphertext<lbcrypto::DCRTPoly> sort_res;
                    bitonic_swap(stage, part, length, comp_res, ciphertext_unsort[i], sort_res, keyPair.secretKey);
                    sort_res_vec.push_back(sort_res);
                }
                end_swap = std::chrono::system_clock::now();
                for (int i = 0; i < ciphertext_unsort.size(); i++)
                {
                    if ((multDepth - sort_res_vec[i]->GetLevel()) < 23)
                    { // 23

                        sort_res_vec[i] = cc->EvalMult(sort_res_vec[i], 1 / precision);
                        start_boot = std::chrono::system_clock::now();

                        ciphertext_unsort[i] = cc->EvalBootstrap(sort_res_vec[i], numIterations, 8);
                        ciphertext_unsort[i] = cc->EvalMult(ciphertext_unsort[i], precision);
                        end_boot = std::chrono::system_clock::now();
                    }
                    else
                    {
                        start_boot = std::chrono::system_clock::now();
                        ciphertext_unsort[i] = sort_res_vec[i];
                        end_boot = std::chrono::system_clock::now();
                    }
                }

                double err_afterBoot = error_estimate_modular(plain_sort_vec[(1 + stage) * stage / 2 + part],
                                                              ciphertext_unsort, keyPair.secretKey, length);

                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end_boot - start_boot).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_swap - start_swap).count()) +
                              double(std::chrono::duration_cast<std::chrono::milliseconds>(end_comp - start_comp).count());
            }
            std::vector<double> plain_sort_res(num_slots, 0);
            for (int i = 0; i < ciphertext_unsort.size(); i++)
            {
                Plaintext plaintextDec;
                cc->Decrypt(keyPair.secretKey, ciphertext_unsort[i], &plaintextDec);
                std::vector<std::complex<double>> enc = plaintextDec->GetCKKSPackedValue();
                for (int j = 0; j < num_slots; j++)
                {
                    plain_sort_res[j] += enc[j].real() * pow(2, 8 * i);
                }
            }
            std::cout << " sortres " << std::endl;
            for (int i = 0; i < pow(2, stage + 1); i++)
            {
                std::cout << "plain:" << plain_sort_vec[(1 + stage) * stage / 2 + stage][i]
                          << " ,sort: " << plain_sort_res[i] << std::endl;
            }
            std::cout << pow(2, stage + 1) << " sort amortize time: " << time_total * (1 << (stage + 1)) / double(length)
                      << "ms" << std::endl
                      << "---------------------------------------------------------" << std::endl
                      << std::endl
                      << std::endl;
        }
    }

    double sync_test_big(int num_column, int plain_bits, int num_slots)
    {
        std::cout << "65536 slots " << num_column << " column sync test " << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
#if NATIVEINT == 128 && !defined(__EMSCRIPTEN__)
        ScalingTechnique rescaleTech = FIXEDAUTO;
        usint dcrtBits = 78;
        usint firstMod = 89;
#else
        usint dcrtBits = 40;
        usint firstMod = 60;
#endif
        parameters.SetScalingModSize(dcrtBits);
        parameters.SetFirstModSize(firstMod);

        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetRingDim(2 * num_slots);
        usint multDepth = (log2(num_slots) + 1) * log2(num_slots) / 2 + 1;
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<std::vector<double>> column;
        std::vector<double> input(length, 0);
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }

        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;

        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;

        double time_total = 0;
        for (int col = 0; col < num_column; col++)
        {
            Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
            auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);
            int part_count = 0;
            for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
            {
                for (int part = 0; stage - part >= 0; part++)
                {
                    int block = 1 << (stage + 1 - part);
                    int part_count = part + ((1 + stage) * stage) / 2;
                    std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;
                    for (int i = 0; i < matrix[part_count].size(); i++)
                    {
                        Plaintext plain_swap_vec = cc->MakeCKKSPackedPlaintext(matrix[part_count][i]);
                        Ciphertext<lbcrypto::DCRTPoly> swap_vec = cc->Encrypt(keyPair.publicKey, plain_swap_vec);
                        swap_matrix.push_back(swap_vec);
                    }
                    start = std::chrono::system_clock::now();
                    auto ct_upper = cc->EvalRotate(ciphertext_unsort, -block / 2);
                    auto ct_lower = cc->EvalRotate(ciphertext_unsort, block / 2);

                    std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
                    res_martrix.push_back(cc->EvalMult(ct_upper, swap_matrix[1]));
                    res_martrix.push_back(cc->EvalMult(ciphertext_unsort, swap_matrix[0]));
                    res_martrix.push_back(cc->EvalMult(ct_lower, swap_matrix[2]));

                    auto temp = cc->EvalAdd(res_martrix[0], res_martrix[1]);
                    ciphertext_unsort = cc->EvalAdd(temp, res_martrix[2]);
                    end = std::chrono::system_clock::now();
                    time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
                }
            }
        }

        std::cout << plain_bits << " bits " << num_column << " column sync time: " << time_total * num_slots / length
                  << "ms" << std::endl;
        return time_total * num_slots / length;
    }
    // 4096
    double sync_test_small(int plain_bits, int num_slots)
    {
        std::cout << num_slots << " records sync test " << std::endl;
        CCParams<CryptoContextCKKSRNS> parameters;

        parameters.SetSecurityLevel(HEStd_128_classic);
        std::uint32_t scaleModSize = 50;
        parameters.SetScalingModSize(scaleModSize);

        SecretKeyDist secretKeyDist = SPARSE_TERNARY;
        parameters.SetSecretKeyDist(secretKeyDist);
        // if (num_slots > 4096)
        //     parameters.SetRingDim(2 * num_slots);

        usint multDepth = (log2(num_slots) + 1) * log2(num_slots) / 2 + 1;
        parameters.SetMultiplicativeDepth(multDepth);
        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(FHE);
        cc->Enable(ADVANCEDSHE);
        // We need to enable Advanced SHE to use the Chebyshev approximation.

        double precision = (1 << (plain_bits - 1)) - 1;

        double lowerBound = -precision - 10;
        double upperBound = precision + 10;
        // double bound      = 3;
        auto keyPair = cc->KeyGen();
        usint ringDim = cc->GetRingDimension();
        int length = ringDim / 2;
        cc->EvalMultKeyGen(keyPair.secretKey);
        std::vector<int> rotstep;
        for (int i = 1; i < num_slots; i *= 2)
        {
            rotstep.push_back(i);
            rotstep.push_back(-i);
        }

        cc->EvalRotateKeyGen(keyPair.secretKey, rotstep);
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<int> message(0, precision);
        std::vector<std::vector<double>> column;
        std::vector<double> input(length, 0);
        for (int i = 0; i < length; i++)
        {
            input[i] = message(engine);
        }

        std::vector<std::vector<double>> plain_sort_vec;
        std::vector<std::vector<std::vector<double>>> matrix = bitonic_test_plain(input, plain_sort_vec, num_slots);
        std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;

        Plaintext plain = cc->MakeCKKSPackedPlaintext(input);
        auto ciphertext_unsort = cc->Encrypt(keyPair.publicKey, plain);
        int part_count = 0;
        std::chrono::system_clock::time_point start, end, start_comp, end_comp, start_swap, end_swap, start_boot, end_boot;
        double time_total = 0;

        for (int stage = 0; num_slots >> stage + 1 > 0; stage++)
        {
            for (int part = 0; stage - part >= 0; part++)
            {
                int block = 1 << (stage + 1 - part);
                int part_count = part + ((1 + stage) * stage) / 2;
                std::vector<Ciphertext<lbcrypto::DCRTPoly>> swap_matrix;
                for (int i = 0; i < matrix[part_count].size(); i++)
                {
                    Plaintext plain_swap_vec = cc->MakeCKKSPackedPlaintext(matrix[part_count][i]);
                    Ciphertext<lbcrypto::DCRTPoly> swap_vec = cc->Encrypt(keyPair.publicKey, plain_swap_vec);
                    swap_matrix.push_back(swap_vec);
                }
                start = std::chrono::system_clock::now();
                auto ct_upper = cc->EvalRotate(ciphertext_unsort, -block / 2);
                auto ct_lower = cc->EvalRotate(ciphertext_unsort, block / 2);

                std::vector<Ciphertext<lbcrypto::DCRTPoly>> res_martrix;
                res_martrix.push_back(cc->EvalMult(ct_upper, swap_matrix[1]));
                res_martrix.push_back(cc->EvalMult(ciphertext_unsort, swap_matrix[0]));
                res_martrix.push_back(cc->EvalMult(ct_lower, swap_matrix[2]));

                auto temp = cc->EvalAdd(res_martrix[0], res_martrix[1]);
                ciphertext_unsort = cc->EvalAdd(temp, res_martrix[2]);
                end = std::chrono::system_clock::now();
                time_total += double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
            }
        }

        std::cout << " sync time: " << time_total * num_slots / length << "ms" << std::endl;
        return time_total * num_slots / length;
    }
}