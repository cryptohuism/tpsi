#ifndef FULL_PSI_H
#define FULL_PSI_H

#include "tfhe_core.h"


#include "full_psi.h"

#include <iostream>
#include "basic_psi.h"
#include "time.h"
#include "tfhe.h"

using namespace std;


/** sample a single item for receiver, partition sigma_N substrings. */
EXPORT void sampleSingleItemforRe(int32_t *receiver_item, const int32_t sigma_N, const int32_t basic_sigma_exp);

/** sample a array items for senders in one bin */
EXPORT void sampleItemsforSeOnebin(int32_t *sender_items, const int32_t sigma_N, const int32_t nb_items, const int32_t basic_sigma_exp);

EXPORT void verify(const int32_t trial, int32_t *receiver_item, int32_t *sender_items, const TFheGateBootstrappingSecretKeySet *keyset);

/** run full psi */
EXPORT void bootsFullPSI(LweSample *reciever_re,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        const int32_t sigma_N,
                        const int32_t nb_se_samples,
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params);

/** handle extract and keyswitch */
EXPORT void handleExtractAndKS(LweSample *results,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        const int32_t sigma_N,
                        const int32_t nb_se_samples,
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params);

/** processe single Extract and keyswitch */
EXPORT void handleSingleExtractAndKS(LweSample *results,
                        TLweSample *reciever_se,
                        const int32_t *sender_items,
                        //const int32_t sigma_N,
                        const int32_t nb_se_substring, // the numbers of substrings in one block
                        const TFheGateBootstrappingSecretKeySet *keyset,
                        const TFheGateBootstrappingParameterSet *params);

/** gate AND and equal */
EXPORT void bootsANDandEqual(LweSample *results, const LweSample *lwe, const LweParams *params, const TFheGateBootstrappingCloudKeySet *bk);

/** gate OR and equal */
EXPORT void bootsORandEqual(LweSample *results, const LweSample *lwe, const LweParams *params, const TFheGateBootstrappingCloudKeySet *bk);


#endif
