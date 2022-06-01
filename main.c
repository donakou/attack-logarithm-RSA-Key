#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "log.h"
#include "rho.h"

gmp_randstate_t rstate;

int main()
{
	int stop =0 ;
	
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,time(NULL));
	mpz_t p,alpha,beta,l,gamma,exposant,ver,ordre;
	int a,b;
	while (stop ==0){
	printf("********* Techniques simples de calcul de logarithme discret *********\n");
	printf("Quelle méthode souhaitez-vous utiliser ?\n");
	printf("[1] Baby-step giant-step\n");
	printf("[2] rho de Pollard\n");
	printf("[3] Pohlig-Hellmam\n");
	scanf("%d",&a);
	switch(a)
	{
		case 1:
			printf("Vous avez choisi la méthode Baby-step giant-step :\n");
			printf("[1] Saisir des paramètres \n[2] Générer des paramètres \n");
			scanf("%d",&b);
			switch(b)
			{
				case 1:
					printf("Saisie des paramètres :\n");
					
					mpz_inits(p,alpha,beta,NULL);
					printf("Un nombre premier p :\n");
					gmp_scanf("%Zd",p);// tester si premier ? ou confiance utilisateur?
					printf("Un générateur du groupe alpha :\n");
					gmp_scanf("%Zd",alpha);
					printf("L'éléments beta dont on veut le logarithme discret :\n");
					gmp_scanf("%Zd",beta);
				
					mpz_init(l);
					if(babystep_giantstep(p,alpha,beta,l)==0)
					{
						gmp_printf("Le log de beta est : %Zd \n",l);
						verif(alpha,l,beta,p);
					}
					else
						printf("Il y a eu une erreur");
					break;
				case 2:
					printf("Générer des paramètres :\n");
					int k;
					printf("Entrez le nombre de chiffres voulus pour p:\n");
					scanf("%d",&k);					
					mpz_inits(p,gamma,NULL);
					gen_groupe(k,p,&gamma);
					gmp_printf("p=%Zd \ngamma=%Zd \n",p,gamma);				
					mpz_init(beta);
					mpz_urandomm(beta,rstate,p);
					if(mpz_cmp_ui(beta,0)==0)					
						mpz_add_ui(beta,beta,1);
					gmp_printf("beta = %Zd \n",beta);
					mpz_init(l);
					if(babystep_giantstep(p,gamma,beta,l)==0)
					{
						gmp_printf("Le log de beta est : %Zd \n",l);
						verif(gamma,l,beta,p);
					}
					else
						printf("Il y a eu une erreur");
					break;
				default:
					printf("Vous n'avez pas rentré un nombre correct\n");
					break;
			}
			break;
		case 2:
			printf("Vous avez choisi la méthode rho de Pollard :\n ");
			printf("[1] pour utiliser des paramètres aléatoires :\n ");
			printf("[2] pour entrer vos paramètres :\n ");
			scanf("%d",&b);
			switch (b)
			{
				case 1 :
					printf("entrez la taille en base 10 du module p\n");
					scanf("%d",&b);
					if (b<0)
					{
						printf(" La taille est négative \n");
						break;
					}
					mpz_inits(p,gamma,beta,exposant,ver,ordre,NULL);
					gen_groupe(b,p,& gamma);
					gmp_printf("le nombre 1er est %Zd \n",p);
					mpz_sub_ui(ordre,p,1);
					mpz_urandomm(beta,rstate,p);
					rho(p,gamma,ordre,beta,rstate,& exposant,0);	//brent
					mpz_powm(ver,gamma,exposant,p);
					if (mpz_cmp(ver,beta)==0) 		// si le calcul est bon
						gmp_printf("%Zd est %Zd à la puissance %Zd \n mod %Zd \n",beta,gamma,exposant,p);
					break;
				case 2:
					printf("Entrez le nombre premier :(le module) \n");
					mpz_inits(p,gamma,beta,exposant,ordre,NULL);
					gmp_scanf("%Zd",&p);
					mpz_sub_ui(ordre,p,1);
					printf("\nEntrez la valeur du générateur de Z*p (assurez vous que s'en est bien un) \n");
					gmp_scanf("%Zd",& gamma);
					printf("\nEntrez le nombre dont vous voulez le log \n");
					gmp_scanf("%Zd",& beta);
					rho(p,gamma,ordre,beta,rstate,& exposant,0);
					gmp_printf("%Zd est %Zd à la puissance %Zd \n mod %Zd \n",beta,gamma,exposant,p);
					break;
				default:
					printf("Vous n'avez pas rentré un nombre correct\n");
				break;
			}
			break;
		case 3:
			printf("Vous avez choisi la méthode Pohlig-Hellmam :\n");	
			test_Pohlig_Hellman(rstate);
			break;
		default:
			printf("Vous n'avez pas rentré un nombre correct\n");
			break;
	}
	
	printf("\n \n");
	printf("Voulez-vous continuer[0] ou arrêter[1] ?\n");
	scanf("%d",&stop);
	}
	return 0;
}
