#ifndef BASCI_PSI_H
#define BASCI_PSI_H

#include "tfhe_core.h"

#include "basic_psi.h"
#include "time.h"
#include "tfhe.h"
#include "polynomials_arithmetic.h"


/** encrypts a item using RLWE encrypt */
EXPORT void bootsRlweEncrypt(TLweSample *result, const int32_t message, const TFheGateBootstrappingSecretKeySet *key);

/** decrypts a boolean */
EXPORT int32_t bootsLweDecrypt(const LweSample *sample, const LweKey *Nkey);

/** sample a bunch of random items */
EXPORT void sampleItems(int32_t *items, const int32_t nb_samples, const int32_t sigma_exp);

/** run the Extract function */
EXPORT void bootsExtract(LweSample *result,
                         const TLweSample *x,
                         const int32_t *sender_items,
                         const int32_t nb_se_samples,
                         const LweParams *params,
                         const TLweParams *rparams);

/** run the KeySwitch function */
EXPORT void bootsKeySwitch(LweSample* result, const LweKeySwitchKey* ks, const LweSample *u);

/** establish a polynomial */
EXPORT void generatePolynomial(IntPolynomial* testvect,
                              const int32_t *sender_items,
                              const int32_t nb_se_samples,
                              const TLweParams *rparams);

/** multiply of tlwe ciphertext and a polynomial */
EXPORT void multiplyPolynomial(TLweSample *result, const TLweSample *tlwe, const IntPolynomial *testvect, const TLweParams *rparams);

/** convert a lwe ciphertext encrypted {0,1} to a lwe ciphertext encrypted {-1,1}, the noise increased by 2 times.
 * 2(a,b)-(0,1/8)
*/
EXPORT void renormalize(LweSample *results, const LweParams *params);

#endif