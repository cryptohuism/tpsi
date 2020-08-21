#include "basic_psi.h"
#include "full_psi.h"

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}


//EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key); //TODO: change the name and put in a .h
//EXPORT void tfhe_createLweBootstrappingKeyFFT(LweBootstrappingKeyFFT* bk, const LweKey* key_in, const TGswKey* rgsw_key);
//EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const LweSample* x);


#ifndef NDEBUG
extern const TLweKey *debug_accum_key;
extern const LweKey *debug_extract_key;
extern const LweKey *debug_in_key;
#endif

int32_t main(int32_t argc, char **argv) {
#ifndef NDEBUG
    cout << "DEBUG MODE!" << endl;
#endif
    const int32_t nb_re_samples = 1; // nb_re_samples is equal to s_y, the number of receiver's items.
    cout << "nb_re_samples is " << nb_re_samples << endl;
    const int32_t nb_se_samples = 200; // nb_re_samples is equal to s_x, the number of sender's items.
    //const int32_t nb_samples = 64;// nb_samples is equal to sx
    const int32_t nb_trials = 10;

    /** sigma_N is the number of substrings divied by 1 item
     * sigma_N = (\omega - log(m')) / (log(N)-1)
     * sigma_N = 13 when \omega = 128, m' = 1024, N = 1024 
     */
    const int32_t sigma_N = 14;

    // generate params 
    int32_t minimum_lambda = 100;
    TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *in_out_params = params->in_out_params;
    // generate the secret keyset
    TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);

    const int32_t basic_sigma_exp = params->tgsw_params->tlwe_params->N /2; // basic_sigma = 9 when N = 1024

    srand((unsigned)time(NULL));
    for (int32_t trial = 0; trial < nb_trials; ++trial) {
        // sample Y for receiver
        int32_t *receiver_item = new int32_t[sigma_N];
        sampleSingleItemforRe(receiver_item, sigma_N, basic_sigma_exp);

        // receiver encrypts Y and sends ciphertext reciever_se to sender
        TLweSample *reciever_se = new_TLweSample_array(sigma_N, params->tgsw_params->tlwe_params);
        // generate inputs (0-->sigma_N=14)
        for (int32_t i = 0; i < sigma_N; ++i) {
            bootsRlweEncrypt(reciever_se + i, *(receiver_item + i), keyset);
        }

        cout << "sizeof(reciever) is " << sizeof(reciever_se) << endl;

        // sample X for sender
        int32_t *sender_items = new int32_t[sigma_N * nb_se_samples];
        sampleItemsforSeOnebin(sender_items, sigma_N, nb_se_samples, basic_sigma_exp);


        /** test
         * let X[0] = Y[0], and receiver would decrypt reciever_re[0] to true.
         *  
         */
        if (1 == rand() % 2) {
            for (int32_t i = 0; i < sigma_N; ++i)
                sender_items[i] = receiver_item[i];
        }

        // evaluate the full PSI tree
        cout << "starting full PSI tree...trial " << trial << endl;
        clock_t begin = clock();

        LweSample *reciever_re = new_LweSample(in_out_params);
        bootsFullPSI(reciever_re, reciever_se, sender_items, sigma_N, nb_se_samples, keyset, params);
        
        clock_t end = clock();
        cout << "finished full PSI tree" << endl;
        cout << "time per full PSI (microsecs)... " << end - begin << endl;


        // verification
        bool mess = 0;
        for (int32_t i = 0; i < nb_se_samples; ++i) {
            // If the substrings of receiver hits sigma_N times, it is considered to be equal.
            int32_t mess2 = 0;
            for (int32_t j = 0; j < sigma_N; ++j) {
                if (*(receiver_item + j) == *(sender_items + i*sigma_N + j)) mess2++;
            }
            if (mess2 == sigma_N) {
                mess = 1;
                break;
            }
        }

        bool out = bootsSymDecrypt(reciever_re, keyset);
        if (out != mess) {
            cout << "ERROR!!! " << trial << "," << " - ";
            cout << t32tod(lwePhase(reciever_re, keyset->lwe_key)) << endl;
        }

    
        delete[] receiver_item;
        delete[] sender_items;
        delete_TLweSample_array(sigma_N, reciever_se);
        delete_LweSample(reciever_re);
    }

    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
