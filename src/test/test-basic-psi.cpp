#include "basic_psi.h"

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
    const int32_t nb_re_samples = 55; // nb_re_samples is equal to s_y, the number of receiver's items.
    const int32_t nb_se_samples = 1000; // nb_re_samples is equal to s_x, the number of sender's items.
    const int32_t nb_trials = 10;

    // generate params 
    int32_t minimum_lambda = 100;
    TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *in_out_params = params->in_out_params;
    // generate the secret keyset
    TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);

    const int32_t sigma_exp = params->tgsw_params->tlwe_params->N /2; // sigma = 9 when N = 1024

    srand((unsigned)time(NULL));
    for (int32_t trial = 0; trial < nb_trials; ++trial) {
        // sample Y for receiver
        int32_t *receiver_items = new int32_t[nb_re_samples];
        sampleItems(receiver_items, nb_re_samples, sigma_exp);

        // receiver encrypts Y and sends ciphertext reciever_se to sender
        TLweSample *reciever_se = new_TLweSample_array(nb_re_samples, params->tgsw_params->tlwe_params);
        // generate inputs (0-->nb_re_samples=63)
        for (int32_t i = 0; i < nb_re_samples; ++i) {
            bootsRlweEncrypt(reciever_se + i, *(receiver_items + i), keyset);
        }

        // sample X for sender
        int32_t *sender_items = new int32_t[nb_se_samples];
        sampleItems(sender_items, nb_se_samples, sigma_exp);


        /** test
         * let sender_items[i] = receiver_items[i] for some random items, and receiver
         * would decrypt reciever_re[i] to true for these i.
         *  
         */
        for (int32_t i = 0; i < nb_re_samples; ++i)
        if(rand() % 3 == 1) sender_items[i] = receiver_items[i];


        // evaluate the EXTRACT function
        cout << "starting extract function...trial " << trial << endl;
        clock_t begin = clock();

        LweSample *reciever_re = new_LweSample_array(nb_re_samples, &params->tgsw_params->tlwe_params->extracted_lweparams);
        for (int32_t i = 0; i < nb_re_samples; ++i) {
            bootsExtract(reciever_re + i, reciever_se + i, sender_items, nb_se_samples, in_out_params, params->tgsw_params->tlwe_params);
            
            /** In fact, we decrypt it using n-lwe secret key instead of N-lwe secret key now if
             * we use Keyswitch function. What's more, this code can also be simplified within
             * generating n-Lwe secret key and bootstrarpping keys.
             * 
             * LweSample *reciever_re = new_LweSample_array(nb_re_samples, in_out_params);
             * bootsKeySwitch(reciever_re + i, keyset->cloud.bk->ks, u + i);
            */
            
        }

        clock_t end = clock();
        cout << "finished extract function tree" << endl;
        cout << "time extract (microsecs)... " << end - begin;
        cout << " with (|X|=" << nb_se_samples << ", |Y|=" << nb_re_samples << ")" << endl;


        /** verification
         * The receiver encrypts messages to tlwe ciphertext while receives lwe ciphertext. 
         * Therefore, it's worth noting that the decryption lwe key 'Nkey' is extracted 
         * from the tlwe key, and the decryption key is a vector while the tlwe key is a
         * polynomial. What's more, they have the same information.  
        */
        LweKey *Nkey = new_LweKey(&params->tgsw_params->tlwe_params->extracted_lweparams);
        tLweExtractKey(Nkey, &keyset->tgsw_key->tlwe_key);

        for (int32_t i = 0; i < nb_re_samples; ++i) {
            bool mess = 0;
            for (int32_t j = 0; j < nb_se_samples; ++j) {
                if (*(receiver_items + i) == *(sender_items + j)) {
                    mess = 1;
                    break;
                }
            }
            
            bool out = bootsLweDecrypt(reciever_re + i, Nkey);

            if (out != mess) {
                cout << "ERROR!!! " << trial << "," << i << " - ";
                cout << t32tod(lwePhase(reciever_re + i, Nkey)) << endl;
            }

        }


        // test full psi begin
        LweSample *test_and = new_LweSample_array(nb_re_samples-1, in_out_params);
        LweSample *test_or = new_LweSample_array(nb_re_samples-1, in_out_params);
        LweSample *test_nand = new_LweSample_array(nb_re_samples-1, in_out_params);

        LweSample *reno = new_LweSample_array(nb_re_samples, in_out_params);
        for (int32_t i = 0; i < nb_re_samples; ++i) {
            bootsKeySwitch(reno + i, keyset->cloud.bk->ks, reciever_re + i);
            //cout << "mu1 is " << lwePhase(reno + i, keyset->lwe_key) << endl;
            renormalize(reno + i, in_out_params);
        }

        //for(int32_t i = 0; i < nb_re_samples; ++i) cout << "mu2 is " << lwePhase(reno + i, keyset->lwe_key) << endl; // after renormalize

        for (int32_t i = 0; i < nb_re_samples - 1; ++i) {
            bootsAND(test_and + i, reno + i, reno + i+1, &keyset->cloud);
            bootsOR(test_or + i, reno + i, reno + i+1, &keyset->cloud);
            bootsNAND(test_nand + i, reno + i, reno + i+1, &keyset->cloud);
        }
        for (int32_t i = 0; i < nb_re_samples - 1; ++i) {
            bool out_and = bootsSymDecrypt(test_and + i, keyset);
            bool out_or = bootsSymDecrypt(test_or + i, keyset);
            bool out_nand = bootsSymDecrypt(test_nand + i, keyset);
            
            if (out_and != (bootsSymDecrypt(reno + i, keyset) && bootsSymDecrypt(reno + i+1, keyset))) {
                cout << "ERROR!!! AND" << trial << "," << i << " - ";
                cout << t32tod(lwePhase(test_and + i, keyset->lwe_key)) << endl;
            }
            cout << "(out_and= " << out_and << ") ?= " << (bootsSymDecrypt(reno + i, keyset) && bootsSymDecrypt(reno + i+1, keyset)) << endl;
            if (out_or != (bootsSymDecrypt(reno + i, keyset) || bootsSymDecrypt(reno + i+1, keyset))) {
                cout << "ERROR!!! OR" << trial << "," << i << " - ";
                cout << t32tod(lwePhase(test_or + i, keyset->lwe_key)) << endl;
            }
            cout << "(out_or = " << out_or << ") ?= " << (bootsSymDecrypt(reno + i, keyset) || bootsSymDecrypt(reno + i+1, keyset)) << endl;
            if (out_nand != 1 - bootsSymDecrypt(reno + i, keyset) * bootsSymDecrypt(reno + i+1, keyset)) {
                cout << "ERROR!!! NAND" << trial << "," << i << " - ";
                cout << t32tod(lwePhase(test_nand + i, keyset->lwe_key)) << endl;
            }
            cout << "(out_nand = " << out_nand << ") ?= " << 1 - bootsSymDecrypt(reno + i, keyset) * bootsSymDecrypt(reno + i+1, keyset) << endl;
        }
        // test full psi end


        delete[] sender_items;
        delete[] receiver_items; 
        delete_TLweSample_array(nb_re_samples, reciever_se);
        delete_LweSample_array(nb_re_samples, reciever_re);
    }

    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
