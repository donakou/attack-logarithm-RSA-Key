#ifndef __log_H__
#define __log_H__

void aff(mpz_t*, int, int);
int AlgoRS(mpz_t **, int, mpz_t);
int AlgoRFN(mpz_t **, int, mpz_t, mpz_t);
void gen(mpz_t, mpz_t *, int,mpz_t*);
void gen_groupe (unsigned int, mpz_t, mpz_t*);//sortie générateur sans factorisation
void gen_groupe2 (unsigned int k,mpz_t p, mpz_t* gamma,mpz_t *ListP,int *sizeP);//sortie générateur avec factorisation
int babystep_giantstep(mpz_t,mpz_t, mpz_t,mpz_t);
void verif(mpz_t,mpz_t,mpz_t,mpz_t);
#endif
