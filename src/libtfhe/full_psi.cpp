#include "full_psi.h"

#include <iostream>
#include "basic_psi.h"
#include "time.h"
#include "tfhe.h"

using namespace std;


/** sample a single item for receiver, partition sigma_N substrings. */
EXPORT void sampleSingleItemforRe(int32_t *receiver_item, const int32_t sigma_N, const int32_t basic_sigma_exp) {
    sampleItems(receiver_item, sigma_N, basic_sigma_exp);
}

/** sample a array items for senders in one bin */
EXPORT void sampleItemsforSeOnebin(int32_t *sender_items, const int32_t sigma_N, const int32_t nb_items, const int32_t basic_sigma_exp) {
    sampleItems(sender_items, sigma_N * nb_items, basic_sigma_exp);
}

EXPORT void verify(const int32_t trial, int32_t *receiver_item, int32_t *sender_items, const TFheGateBootstrappingSecretKeySet *keyset) {
    //TODO
}

/** run full psi */
EXPORT void bootsFullPSI(LweSample *reciever_re,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        const int32_t sigma_N,
                        const int32_t nb_se_samples,
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params) {
    
    const int32_t N = params->tgsw_params->tlwe_params->N;
    int32_t lines = ceil(double(nb_se_samples)*20/N);
    // sample temporary lwes
    LweSample *lwes_temp = new_LweSample_array(lines * sigma_N, params->in_out_params);
    handleExtractAndKS(lwes_temp, reciever_se, sender_items, sigma_N, nb_se_samples, keyset, params);
    
    for (int32_t i = 0; i < lines - 1; ++i) {
        for (int32_t j = 0; j < sigma_N - 1; ++j) {
            bootsANDandEqual(lwes_temp + i*sigma_N + j + 1, lwes_temp + i*sigma_N + j, params->in_out_params, &keyset->cloud);
        }
    }
    for (int32_t i = 0; i < lines - 1; ++i) {
        // note: the vaild lwe ciphertext is in the location of sigma_N - 1 for every lines
        bootsORandEqual(lwes_temp + (i+1)*sigma_N + sigma_N - 1, lwes_temp + i*sigma_N + sigma_N - 1, params->in_out_params, &keyset->cloud);
    }

    lweCopy(reciever_re, lwes_temp + (lines-1)*sigma_N + sigma_N - 1, params->in_out_params);

    delete_LweSample_array(lines * sigma_N, lwes_temp);
}

/** handle extract and keyswitch */
EXPORT void handleExtractAndKS(LweSample *results,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        const int32_t sigma_N,
                        const int32_t nb_se_samples,
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params) {
    
    const int32_t N = params->tgsw_params->tlwe_params->N;
    int32_t lines = ceil(double(nb_se_samples)*20/N);
    // sample temporary lwes
    // LweSample *lwes_temporary = new_LweSample_array(lines * sigma_N, params->in_out_params);
    
    for (int32_t i = 0; i < lines - 1; ++i) {
        /** The lines expect the last one contains N/20 substrings, while the last
         * one cantains nb_se_samples%(N/20). Therefore, we handle it separately.
        */
        for (int32_t j = 0; j < sigma_N; ++j) {
            // sample u and keyswitch to obtain result(reciever_re)
            
            // partition, create virtual items which is one of blocks.
            int32_t *virtual_items = new int32_t[N/20];
            for (int32_t index = 0; index < N/20; ++index)
                virtual_items[index] = *(sender_items + i * sigma_N * (N/20) + index * sigma_N + j);

            handleSingleExtractAndKS(results + i * sigma_N + j, reciever_se + j, virtual_items, N/20, keyset, params);

            delete[] virtual_items;
        }
    }

    // handle last line
    int32_t i = lines - 1;
    for (int32_t j = 0; j < sigma_N; ++j) {
        //cout << "nb_se_samples%(N/20) = " << nb_se_samples%(N/20) << endl;
        int32_t *virtual_items_last = new int32_t[nb_se_samples%(N/20)];
        for (int32_t index = 0; index < nb_se_samples%(N/20); ++index)
            virtual_items_last[index] = *(sender_items + i * sigma_N * (N/20) + index * sigma_N + j);
        handleSingleExtractAndKS(results + i * sigma_N + j, reciever_se + j, virtual_items_last, nb_se_samples%(N/20), keyset, params);
        
        delete[] virtual_items_last;
    }

}

/** processe single Extract and keyswitch */
EXPORT void handleSingleExtractAndKS(LweSample *results,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        //const int32_t sigma_N,
                        const int32_t nb_se_substring, // the numbers of substrings in one block
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params) {
    
    LweSample *u = new_LweSample(&params->tgsw_params->tlwe_params->extracted_lweparams);

    bootsExtract(u, reciever_se, sender_items, nb_se_substring, params->in_out_params, params->tgsw_params->tlwe_params);
    bootsKeySwitch(results, keyset->cloud.bk->ks, u);

    //Torus32 mu = lwePhase(results, keyset->lwe_key);
    //cout << "mu = " << mu << endl;

    renormalize(results, params->in_out_params);

    //Torus32 mu2 = lwePhase(results, keyset->lwe_key);
    //cout << "mu2 = " << mu2 << endl;
    
    delete_LweSample(u);
}

/** gate AND and equal */
EXPORT void bootsANDandEqual(LweSample *results, const LweSample *lwe, const LweParams *params, const TFheGateBootstrappingCloudKeySet *bk) {
    LweSample *u = new_LweSample(params);
    bootsAND(u, results, lwe, bk);

    lweCopy(results, u, params);
    delete_LweSample(u);
}

/** gate OR and equal */
EXPORT void bootsORandEqual(LweSample *results, const LweSample *lwe, const LweParams *params, const TFheGateBootstrappingCloudKeySet *bk) {
    LweSample *u = new_LweSample(params);
    bootsOR(u, results, lwe, bk);

    lweCopy(results, u, params);
    delete_LweSample(u);
}