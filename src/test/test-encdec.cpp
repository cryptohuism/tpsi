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

#include "basic_psi.h"

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
    //const int32_t nb_samples = 64;
    const int32_t nb_trials = 10;

    // generate params 
    int32_t minimum_lambda = 100;
    TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
    const LweParams *in_out_params = params->in_out_params;
    // generate the secret keyset
    TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);

    int32_t N = params->tgsw_params->tlwe_params->N; 

    for (int32_t trial = 0; trial < nb_trials; ++trial) {
        srand((unsigned)time(NULL));

        // **********************************************************************************
        // ********************      Test lwe encryption           **************************
        // **********************************************************************************
        /**
        cout << "starting lwe Encryption tree...trial " << trial << endl;
        LweSample *test_in = new_LweSample_array(nb_samples, in_out_params);
        int32_t *messages1 = new int32_t[nb_samples];
        bool *mess1 = new bool[nb_samples];

        // generate inputs (64-->127)
        for (int32_t i = 0; i < nb_samples; ++i) {
            *(messages1 + i) = rand() % 2;
            cout << *(messages1 + i);
            bootsSymEncrypt(test_in + i, *(messages1 + i), keyset);
            *(mess1 + i) = bootsSymDecrypt(test_in + i, keyset);
        }
        delete_LweSample_array(nb_samples, test_in);
        cout << "finished lwe Encryption tree" << endl;
        */

        // **********************************************************************************
        // ********************      Test rlwe encryption          **************************
        // **********************************************************************************
        /** test rlwe encryption, just a single encryption */
        cout << "starting tlwe Encryption tree...trial " << trial << endl;
        //int32_t message2 = 32;
        TLweSample *test_in2 = new_TLweSample(params->tgsw_params->tlwe_params);
        //bootsRlweEncrypt(test_in2, message2, keyset);
        TorusPolynomial *messages_poly = new_TorusPolynomial(params->tgsw_params->tlwe_params->N);
        //double alpha = keyset->params->in_out_params->alpha_min;
        double alpha = 0.0000298;

        
        TorusPolynomial *mess2 = new_TorusPolynomial(params->tgsw_params->tlwe_params->N);
        for (int32_t i = 0; i < messages_poly->N; ++i)
            messages_poly->coefsT[i] = 0;
        messages_poly->coefsT[0] = 536870912; messages_poly->coefsT[4] = -536870912;
        //836870912

        TorusPolynomial *phi = new_TorusPolynomial(params->tgsw_params->tlwe_params->N);
        tLweSymEncrypt(test_in2, messages_poly, alpha, &keyset->tgsw_key->tlwe_key);
        tLwePhase(phi, test_in2, &keyset->tgsw_key->tlwe_key); // PHASE COMUPTATION
        tLweSymDecrypt(mess2, test_in2, &keyset->tgsw_key->tlwe_key, 8);
        cout << "finished tlwe Encryption tree" << endl;

        
        /**   new lwe key, this key is N length which is used for TLWE origin */
        LweParams *Nparam = new_LweParams(N, 0.0000244, 0.012467000000000001);
        LweKey *Nkey = new_LweKey(Nparam);
        for (int32_t i = 0; i < N; ++i)
            Nkey->key[i] = keyset->tgsw_key->tlwe_key.key->coefs[i];
        


        // **********************************************************************************
        // ********************      test Extract function         **************************
        // **********************************************************************************
        //LweParams *extract_params = params->tgsw_params->tlwe_params->extracted_lweparams;
        LweSample *u = new_LweSample(&params->tgsw_params->tlwe_params->extracted_lweparams);
        LweSample *reciever_re = new_LweSample(in_out_params);
        
        tLweExtractLweSample(u, test_in2, &params->tgsw_params->tlwe_params->extracted_lweparams, params->tgsw_params->tlwe_params);
        
        // Nkey2
        const LweParams *extract_params = &params->tgsw_params->tlwe_params->extracted_lweparams;
        //LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
        const TLweKey *accum_key = &keyset->tgsw_key->tlwe_key;
        LweKey *Nkey2 = new_LweKey(extract_params);
        tLweExtractKey(Nkey2, accum_key);

        cout << "Nmu is " << lwePhase(u, Nkey) << endl;
        cout << "Nmu2 is " << lwePhase(u, Nkey2) << endl;
        lweKeySwitch(reciever_re, keyset->cloud.bk->ks, u);
        bool out = bootsSymDecrypt(reciever_re, keyset);

        cout << "mu is " << lwePhase(reciever_re, keyset->lwe_key) << endl;
        cout << "out is " << out << endl;


        // the second lwe ciphertext
        LweSample *reciever_re2 = new_LweSample(in_out_params);
        int32_t reciever_nu2 = rand() % 2;
        bootsSymEncrypt(reciever_re2, reciever_nu2, keyset);
        cout << "reciever_nu2 is " << reciever_nu2 << endl;



        LweSample *AND = new_LweSample(in_out_params);
        LweSample *OR = new_LweSample(in_out_params);
        LweSample *NAND = new_LweSample(in_out_params);

        bootsAND(AND, reciever_re, reciever_re2, &keyset->cloud);
        bootsOR(OR, reciever_re, reciever_re2, &keyset->cloud);
        bootsNAND(NAND, reciever_re, reciever_re2, &keyset->cloud);

        cout << "AND is " << bootsSymDecrypt(AND, keyset) << endl;
        cout << "OR is " << bootsSymDecrypt(OR, keyset) << endl;
        cout << "NAND is " << bootsSymDecrypt(NAND, keyset) << endl;

        /**
        cout << "mess2 is (";
        for (int32_t i = 0; i < params->tgsw_params->tlwe_params->N; ++i)
            cout << mess2->coefsT[i] << ",";
        cout << ")." << endl;
        */

        cout << "finished trial " << trial << endl;

        delete_TLweSample(test_in2);
        delete_TorusPolynomial(messages_poly);
        delete_TorusPolynomial(mess2);
        delete_TorusPolynomial(phi);
        delete_LweParams(Nparam);
        delete_LweKey(Nkey2);
        delete_LweSample(u);
        delete_LweSample(reciever_re);

    }

    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
