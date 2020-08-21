#include "basic_psi.h"
#include "time.h"
#include "tfhe.h"
#include "polynomials_arithmetic.h"


/** encrypts a item using RLWE encrypt */
EXPORT void bootsRlweEncrypt(TLweSample *result, const int32_t message, const TFheGateBootstrappingSecretKeySet *key) {
    double alpha = key->params->in_out_params->alpha_min;
    int32_t N = key->params->tgsw_params->tlwe_params->N;

    TorusPolynomial *message_poly = new_TorusPolynomial(N);
    Torus32 _1s8 = modSwitchToTorus32(1, 8);

    for(int32_t i = 0; i < N; ++i) message_poly->coefsT[i] = 0;
    message_poly->coefsT[N-message] = -_1s8;
    tLweSymEncrypt(result, message_poly, alpha, &key->tgsw_key->tlwe_key);
}

/** decrypts a boolean */
EXPORT int32_t bootsLweDecrypt(const LweSample *sample, const LweKey *Nkey) {
    Torus32 mu = lwePhase(sample, Nkey);
    /** It's worth noting that we use modSwitchToTorus32(1, 16) insead of 0 on accout of
     * the decrypt value is (0, 1/8) instead of (-1/8, 1/8). In fact, origin lwe decryp
     * tion use `return (mu > 0 ? 1 : 0);`. However, both of them meet safety requirements.
    */
    return (mu > modSwitchToTorus32(1, 16) ? 1 : 0);
}

/** sample a bunch of random items */
EXPORT void sampleItems(int32_t *items, const int32_t nb_samples, const int32_t sigma_exp) {
    for (int32_t i = 0; i <= nb_samples; i++) items[i] = rand() % sigma_exp;
}

/** run the Extract function */
EXPORT void bootsExtract(LweSample *result,
                         const TLweSample *x,
                         const int32_t *sender_items,
                         const int32_t nb_se_samples,
                         const LweParams *params,
                         const TLweParams *rparams) {
    // test polynomial 
    IntPolynomial *testvect = new_IntPolynomial(rparams->N);
    generatePolynomial(testvect, sender_items, nb_se_samples, rparams);

    // multiply test polynomial
    TLweSample *mulx = new_TLweSample(rparams);
    multiplyPolynomial(mulx, x, testvect, rparams);

    // extract a lwe ciphertext
    const LweParams *extract_params = &rparams->extracted_lweparams;
    tLweExtractLweSample(result, mulx, extract_params, rparams);
}

/** run the KeySwitch function */
EXPORT void bootsKeySwitch(LweSample* result, const LweKeySwitchKey* ks, const LweSample *u) {
    lweKeySwitch(result, ks, u);
}

/** establish a polynomial */
EXPORT void generatePolynomial(IntPolynomial* testvect,
                              const int32_t *sender_items,
                              const int32_t nb_se_samples,
                              const TLweParams *rparams) {
    
    for (int32_t i = 0; i < rparams->N; ++i) testvect->coefs[i] = 0; // initialization
    for (int32_t i = 0; i < nb_se_samples; ++i) {
        int32_t index = *(sender_items + i);
        testvect->coefs[index] = 1;
    }
}


/** multiply of tlwe ciphertext and a polynomial */
EXPORT void multiplyPolynomial(TLweSample *result, const TLweSample *tlwe, const IntPolynomial *testvect, const TLweParams *rparams) {
    /*
    TorusPolynomial *tlwe_a = new_TorusPolynomial(rparams->N);
    TorusPolynomial *tlwe_b = new_TorusPolynomial(rparams->N);
    torusPolynomialCopy(tlwe_a, result->a);
    torusPolynomialCopy(tlwe_b, result->b);
    */
    
    // results = (a,b)*testvect
    torusPolynomialMultKaratsuba(result->a, testvect, tlwe->a);
    torusPolynomialMultKaratsuba(result->b, testvect, tlwe->b);
}

/** convert a lwe ciphertext encrypted {0,1} to a lwe ciphertext encrypted {-1,1}, the noise increased by 2 times.
 * 2(a,b)-(0,1/8)
*/
EXPORT void renormalize(LweSample *results, const LweParams *params) {
    LweSample *lwe_temp = new_LweSample(params);
    Torus32 _1s8 = modSwitchToTorus32(1, 8);

    for (int32_t i = 0; i < params->n; ++i) lwe_temp->a[i] = 2*results->a[i];
    lwe_temp->b = 2 * results->b - _1s8;

    lweCopy(results, lwe_temp, params);
}