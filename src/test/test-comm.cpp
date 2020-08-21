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

#define nb_X 5

using namespace std;


struct LweSampleArry {
    LweSample send[nb_X];
};

struct TLweSampleArry {
    TLweSample receive[nb_X];
};

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
    const int32_t nb_re_samples = nb_X; // nb_re_samples is equal to s_y, the number of receiver's items.
    const int32_t nb_se_samples = 10; // nb_re_samples is equal to s_x, the number of sender's items.
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

        // evaluate the EXTRACT function
        cout << "starting extract function...trial " << trial << endl;

        LweSample *reciever_re = new_LweSample_array(nb_re_samples, &params->tgsw_params->tlwe_params->extracted_lweparams);
        for (int32_t i = 0; i < nb_re_samples; ++i) {
            bootsExtract(reciever_re + i, reciever_se + i, sender_items, nb_se_samples, in_out_params, params->tgsw_params->tlwe_params);     
        }
        cout << "finished extract function tree" << endl;

        // **********************************************************************************
        // **************************** communication cost **********************************
        // **********************************************************************************
        // receiver send an array of tlwe ciphertexts
        //struct LweSampleArray *lwe_arry_for_write = (struct LweSampleArray *)malloc(sizeof(struct LweSampleArray));
        //struct TLweSampleArray *tlwe_arry_for_write = (struct TLweSampleArray *)malloc(sizeof(struct TLweSampleArray));

        struct TLweSample *tlwe_for_write = (struct TLweSample *)malloc(sizeof(struct TLweSample));
        struct TLweSample *tlwe_for_read = (struct TLweSample *)malloc(sizeof(struct TLweSample));
        if (tlwe_for_write == NULL || tlwe_for_read == NULL) {
            perror("Failed to apply memory!");
            exit(EXIT_FAILURE);
        }

        FILE *fp;  
        // write in
        if ((fp = fopen("./send.txt", "w")) == NULL) {
            printf("Failed to open file!");
            exit(EXIT_FAILURE);
        }
        tlwe_for_write = reciever_se;
        fwrite(tlwe_for_write, sizeof(struct TLweSample), 1, fp);
        fclose(fp);
        // read
        if ((fp = fopen("./send.txt", "r")) == NULL) {
            perror("Failed to open file!");
            exit(EXIT_FAILURE);
        }   
        fread(tlwe_for_read, sizeof(struct TLweSample), 1, fp);
        cout << "The a[0] of tlwe is " << *tlwe_for_read->a[0].coefsT << endl;
        fclose(fp);
        
        // receiver receive an array of lwe ciphertexts
        struct LweSample *lwe_for_write = (struct LweSample *)malloc(sizeof(struct LweSample));
        struct LweSample *lwe_for_read = (struct LweSample *)malloc(sizeof(struct LweSample));
        if (lwe_for_write == NULL || lwe_for_read == NULL) {
            perror("Failed to apply memory!");
            exit(EXIT_FAILURE);
        }
 
        // write in
        if ((fp = fopen("./receive.txt", "w")) == NULL) {
            printf("Failed to open file!");
            exit(EXIT_FAILURE);
        }
        lwe_for_write = reciever_re;
        fwrite(lwe_for_write, sizeof(struct LweSample), 1, fp);
        fclose(fp);
        // read
        if ((fp = fopen("./receive.txt", "r")) == NULL) {
            perror("Failed to open file!");
            exit(EXIT_FAILURE);
        }   
        fread(lwe_for_read, sizeof(struct LweSample), 1, fp);
        cout << "The a[0] of lwe is " << lwe_for_read->a[0] << endl;
        fclose(fp);





        delete[] sender_items;
        delete[] receiver_items; 
        delete_TLweSample_array(nb_re_samples, reciever_se);
        delete_LweSample_array(nb_re_samples, reciever_re);
    }

    delete_gate_bootstrapping_secret_keyset(keyset);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
