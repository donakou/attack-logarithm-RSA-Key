#ifndef __rho_H__
#define __rho_H__

void fct (mpz_t module , mpz_t gene, mpz_t ord, mpz_t alpha , mpz_t* x,mpz_t* a,mpz_t* b);
int rho( mpz_t module , mpz_t gene, mpz_t ord , mpz_t alpha,gmp_randstate_t rstate,mpz_t* exposant,int mode ); // alpha est le nb dont on veut le log; // le mode 1 pour la comparaison; 0 pour le mode travail
void Pohlig_Hellmann(mpz_t module,mpz_t gene,mpz_t alpha, gmp_randstate_t rstate,mpz_t * fact , int taille, mpz_t* exposant);// fact est la facto de p-1
int test1 (gmp_randstate_t);
void test_Rho_Pollard(gmp_randstate_t );
int test_Pohlig_Hellman(gmp_randstate_t );
void Euristique_Rho(gmp_randstate_t );
void Euristique_Pohlig(gmp_randstate_t );
#endif
