#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <gmp.h>
#include "rho.h"
void fct (mpz_t module , mpz_t gene, mpz_t ord, mpz_t alpha , mpz_t* x,mpz_t* a,mpz_t* b)
{
	mpz_t q;
	mpz_t result;
	mpz_init(q);
	mpz_init(result);
	
	mpz_cdiv_q_ui(q ,module,3); // on a int (p/3)
	if (mpz_cmp(*x,q)<=0)
	{
		mpz_mul(*x, alpha,*x);
		mpz_mod(*x,*x,module);
		mpz_add_ui(*a,*a,1);
		mpz_mod(*a,*a,ord); 
	}
		
	else
	{
		mpz_mul_ui(q,module,2);    // 2p
		mpz_cdiv_q_ui(q,q,3); // on a int (2p/3)
		if (mpz_cmp(*x,q)<=0)  // *x < ou = 2p/3		
		{
			mpz_mul(*x, *x,*x);
			mpz_mod(*x,*x,module);
			mpz_mul_ui(*a,*a,2);
			mpz_mod(*a,*a,ord);
			mpz_mul_ui(*b,*b,2);
			mpz_mod(*b,*b,ord);
			
		}
		else
		{
			mpz_mul(*x,gene,*x);
			mpz_mod(*x,*x,module);
			mpz_add_ui(*b,*b,1);
			mpz_mod(*b,*b,ord);
		}
	}

	mpz_clear(q);
	mpz_clear(result);
}


/******************************************************************************
								RHO
******************************************************************************/
int rho( mpz_t module , mpz_t gene, mpz_t ord , mpz_t alpha,gmp_randstate_t rstate,mpz_t* exposant,int mode )
{
	if(mpz_divisible_p(alpha,module)!=0)
		return 1;
	int cp=0;//cpt=0;u_1,u_2,tps=0;
	mpz_t x,y,ax,inx,bx,ay,by;
	mpz_inits(x,y,ax,ay,bx,by,inx,NULL);
	//FILE* floyd,* brent,* logd;
	mpz_set_ui(ax,0);
	mpz_set_ui(ay,0);
	mpz_urandomm(inx,rstate,ord); // on prendra plutôt un exposant de g aléatoire
	mpz_set(bx,inx);
	mpz_set(by,bx);
	mpz_powm(x,gene,bx,module);
	mpz_powm(y,gene,by,module);

	/****************  détermination de collision  *******************/

	/*********	floyd  	***********/
	if (1)//(mode==1) // si l'on est en mode comparaison
	{	
		
		do 
		{	

			fct(module,gene,ord,alpha,&x,&ax,&bx);
			fct(module,gene,ord,alpha,&y,&ay,&by);
			fct(module,gene,ord,alpha,&y,&ay,&by);
			//cpt++;
			//printf("j'y suis \n");
		}
		while(mpz_cmp(x,y)!=0);

		
	}
	
	else
	{
		/**********	Brent 	*******/

		// on reprend au début 
		mpz_set_ui(ax,0);
		mpz_set_ui(ay,0);
		mpz_set(bx,inx);
		mpz_set(by,bx);
		mpz_powm(x,gene,bx,module);
		mpz_powm(y,gene,by,module);
		int pow =1,lam=1;
		/*if (mode==1)
		{
			brent=fopen("tps_brent","a+");
			u_1=time(NULL);
		}*/

		fct(module,gene,ord,alpha,&y,&ay,&by);  // y=f(x)
		while (mpz_cmp(x,y) != 0)
		{
			if (pow==lam)
			{
				mpz_set(x,y);
				mpz_set(ax,ay);
				mpz_set(bx,by);
				pow=pow*2;
				lam =0;
			}
			fct(module,gene,ord,alpha,&y,&ay,&by);
			lam++;
			cp++;
		//printf(" %d",cp);

		}		
	}

	
	//if (mode==1) gmp_printf (" collision trouvée en \n   FLOYD : %d  appels de f   \n   BRENT : %d appels de f  \n",3*cpt,cp);


	/******************* détermination de l'exposant **************/
	mpz_t diffa,diffb,tmp,gcd,h,inv,m,sum,al;
	mpz_inits(diffa,diffb,tmp,NULL);
	mpz_inits(gcd,h,inv,m,sum,al,NULL);
	int d,k;
	mpz_sub(diffa,ax,ay);
	mpz_sub(diffb,by,bx);
	mpz_gcd(gcd,diffa , ord);
	if (1)	//(mpz_cmp_ui(gcd,1000000)<0)
	{
		d= mpz_get_ui(gcd);
		mpz_cdiv_q(h,ord,gcd);
		mpz_invert(inv,diffa,h);
		mpz_mul(m,diffb,inv);
		mpz_mod(m,m,h);
		for (k=0;k<d;k++)
		{
			mpz_mul_ui(tmp,h,k);
			mpz_add(sum,m,tmp) ;
			mpz_mod(sum,sum,ord);  // on prend l'exposant modulo odr
			mpz_powm(tmp,gene,sum,module);

			mpz_mod(al,alpha,module);
			if (mpz_cmp(tmp,al)==0)
			{
				mpz_set(*exposant,sum);
				break;
			}
		}
	}

	else 
	{	printf ("pgcd trop grand ");
	}
	/*
	mpz_clears(x,y,ax);
	mpz_clear(inx);
	mpz_clear(bx);
	mpz_clears(ay,by);
	mpz_clears(diffa,diffb,tmp);
	mpz_clear(gcd);
	mpz_clears(h,inv,m);
	mpz_clear(sum); 
	mpz_clear(al);*/

	return 0;
}



/******************************************************************************
*						Pohlig_Hellmann                                        *
*******************************************************************************/


void Pohlig_Hellmann(mpz_t module,mpz_t gene,mpz_t alpha, gmp_randstate_t rstate,mpz_t * fact , int taille, mpz_t* exposant)

{
	mpz_t* l;
	mpz_t aw,acc,m,tmp,mp,h,inv,odr,p,u,frac,cpt;
	int i,j,e;
	mpz_inits(aw,acc,m,tmp,mp,h,inv,odr,p,u,frac,cpt,NULL);
	mpz_set_ui(acc,1);
	mpz_set_ui(aw,1);
	mpz_sub_ui(mp,module,1);						// l'on a déjà sa factorisation dans fact
	l=(mpz_t*)malloc((taille/2)*sizeof(mpz_t)); 	// stoque les différentes valeurs auquels est congru le log mod les dv de l'ord

	//printf("lataille du tableau est %d\n",taille);
	for (i=0;i<taille/2;i++)						// dans fact on a le nb 1er et son nb d'occurence
	{
		mpz_init(l[i]);
		mpz_set_ui(l[i],0);
	}

	for (i=0;i<taille/2;i++)				// pour chaque diviseur premier p de moule-1
	{
		mpz_set(p,fact[2*i]);
		e= mpz_get_ui(fact[2*i+1]); 				// on réccupère p puis son exposant 

		mpz_cdiv_q(tmp,mp,p);						// on définit le géné du sg d'ordre p (h= gene^((mod-1)/p))
		mpz_powm(h,gene,tmp,module);

		mpz_set(odr,p);								// on met à jour l'ordre du new sgp
		mpz_set_ui(acc,1);

		
		mpz_set_ui(cpt,0);						
		//while(mpz_cmp(fact[2*i+1],cpt) > 0)		//for cpt de 0 à e-1 
		for(j=0;j<e;j++)
		{
			/***aw = (alpha*acc)^((mod-1)/p_i^(j+1))****/
			mpz_pow_ui(tmp,fact[2*i],j+1);
			mpz_cdiv_q(tmp,mp,tmp);
			mpz_mul(aw,alpha,acc);
			
			mpz_powm(aw,aw,tmp,module);
			rho(module,h,odr,aw,rstate,& m,0); 	// m = exposant de aw dans le sg d'ordre p 

			//gmp_printf("  dans le sg d'ordre %Zd \n",odr);
			
			/**** l[i]= l[i]+m.p^j*****/
			mpz_pow_ui(tmp,p,j);
			mpz_mul(tmp,m,tmp);
			mpz_add(l[i],l[i],tmp);

			// mise à jour de l'accumulateur acc = acc.gene^(-m*p^j)= acc*(gene^tmp)^-1
			mpz_powm(tmp,gene,tmp,module);
			mpz_invert(tmp,tmp,module);
			mpz_mul(acc,acc,tmp);
			mpz_mod(acc,acc,module);
			

		}
		//gmp_printf(" dans le sg d'ordre %Zd l'exposant est congru à %Zd\n",odr,l[i]);
	} // à la fin de la boucle log(alpha) = l[i]mod (P_i) ^(e_i) 

	//********************** CRT  ****************************************
	/* le log(alpha) = l1.(p2^(e2)p3^(e3)...pr^(er))*inv_mod_p1^(e1)[p2^(e2)p3^(e3)...pr^(er)]+
					   l2.(p1^(e1)p3^(e3)...pr^(er))*inv_mod_p2^(e2)[p1^(e1)p3^(e3)...pr^(er)]+...
					  +lr.(p1^(e1)p2^(e2)...(pr-1)^(er-1) * inv_mod_pr^(er)[p1^(e1)p2^(e2)...(pr-1)^(er-1)] mod (module-1)	*/

	
	mpz_set_ui(m,0);
	for (i=0;i<taille/2;i++)
	{
		// m = m + l[i]*[(module-1)/(p_i^(ei))]*inv_mod_pi^(ei)
		e= mpz_get_ui(fact[2*i+1]);
		mpz_pow_ui(u,fact[2*i],e);			// p_i^(ei)
		mpz_cdiv_q(frac,mp,u);				
		mpz_invert(inv,frac,u);
		mpz_mul(tmp,inv,frac);
		mpz_mul(tmp,tmp,l[i]);
		mpz_add(m,m,tmp);

	}
	mpz_mod(m,m,mp); // on prend l'exposant modulo module-1 (card du gp multiplicatif)
	mpz_set(*exposant,m);
	mpz_clears(aw,acc,m,tmp,mp,h,inv,odr,p,u,frac);
}




/****************************************************************************
*								Test 1										*
***************************************************************************/
int test1 (gmp_randstate_t rstate)
{

	mpz_t p,g,j,x,odr,exposant;
	//mpz_t fact[4];
	mpz_inits(p,g,j,x,odr,exposant,NULL);
	mpz_set_str(p,"689740383853",10);
	mpz_set_str(g,"394319546118",10);
	mpz_powm_ui(j,g,996565889,p);
	mpz_set_str(x,"123455454546",10);
	mpz_sub_ui(odr,p,1);
	rho(p,g,odr,j,rstate,& exposant,0); // le 1 pour le mode comparaison; 0 pour le mode travail
	gmp_printf("la valeur  %Zd correspond à g exposant %Zd \n ",j,exposant);

	/******* test pohlig-hellmann ********
	mpz_inits(fact[0],fact[1],fact[2],fact[3],NULL);
	mpz_set_ui(p,101);
	mpz_set_ui(fact[0],2);
	mpz_set_ui(fact[1],2);
	mpz_set_ui(fact[2],5);
	mpz_set_ui(fact[3],2);
	mpz_set_ui(g,2);
	mpz_powm_ui(j,g,198,p);
	Pohlig_Hellmann(p,g,j,rstate,fact,4,& exposant);
	gmp_printf("la valeur  %Zd correspond à g exposant %Zd \n ",j,exposant);
	*/
	return 0 ;
}


/************* Test Rho de polard ********************************
******************************************************************/
void test_Rho_Pollard(gmp_randstate_t rstate)
{
	
	int i;
	mpz_t pi[4],gi[4],b,exposant,ordre, verif;
	
	for (i=0;i<4;i++)
	{
		mpz_inits(pi[i],gi[i],NULL);
	}
	mpz_set_str(pi[0],"689740383853",10);
	mpz_set_str(pi[1],"12263445054821",10);
	mpz_set_str(pi[2],"7100573378083121",10);
	mpz_set_str(pi[3],"239108134213568687",10);

	mpz_set_str(gi[0],"394319546118",10);
	mpz_set_str(gi[1],"6959266850417",10);
	mpz_set_str(gi[2],"4195130626214895",10);
	mpz_set_str(gi[3],"11196532448230429",10);

	mpz_inits(b,exposant,ordre,verif,NULL);
	mpz_set_ui(b,2010);
	mpz_set_ui(exposant,0);

	printf("/**************  test Rho de polard ****************/ \n");
	
	for (i=0;i<4;i++)
	{
		mpz_sub_ui(ordre,pi[i],1); // on donne l'ordre du gp multiplicatif
		rho(pi[i],gi[i],ordre,b,rstate,& exposant,0);
		mpz_powm(verif,gi[i],exposant,pi[i]);
		if (mpz_cmp(verif,b)==0) 		// si le calcul est bon
			gmp_printf("%Zd est %Zd à la puissance %Zd \n mod %Zd \n",b,gi[i],exposant,pi[i]);
	}
	
}







/*************** Test Pohlig- Hellman *********************
***********************************************************/


int test_Pohlig_Hellman(gmp_randstate_t rstate)
{
	printf("/*************** Test Pohlig- Hellman *********************/ \n");
	mpz_t q,p_1;
	mpz_t fact[200];
	mpz_init(p_1);
	mpz_set_str(p_1,"93829103457701351047678150545874455874415200013145296892799092824465729890236733039644041875602059462560643180671763707936371771172424325634235538844763279109070226245563765394479104625421576563590785443467192831069583054627852233430288407576160818322777666044553172636695128063366168185951512357172769367900",10);
	gmp_printf("p-1 = %Zd \n",p_1);	
	int i;
	i=0;
	mpz_init(q);
	mpz_set_ui(q,2);
	mpz_init(fact[0]);
	mpz_init(fact[1]);
	mpz_set_ui(fact[0],1);
	mpz_set_ui(fact[1],0);
	// on cherche les facteurs 1er de p-1
	printf("Les facteurs de p-1 sont : \n");  
	while ((mpz_cmp(p_1,q)>=0) )//(mpz_cmp_ui(q,2000000)<0) &&)
	{
		while ((mpz_divisible_p(p_1,q)!=0) && (mpz_cmp(p_1,q)>=0)) // tant que p-1 est dv par p
		{	
			if (mpz_cmp(q,fact[2*i])==0)			// si on est sur le mê facteur
			{
				mpz_divexact(p_1,p_1,q);
				mpz_add_ui(fact[2*i+1],fact[2*i+1],1);
			}
			else									// on a un new fact
			{
				mpz_set(fact[2*i],q);
				mpz_set_ui(fact[2*i+1],1);
				mpz_divexact(p_1,p_1,q);

			}	
			// tester la taille du tableau
		}

		mpz_nextprime(q,q);
		if(mpz_divisible_p(p_1,q)!=0)		// q divise p_1
				{	      
					i++;
					mpz_init(fact[2*i]);
					mpz_init(fact[2*i+1]);
					mpz_set_ui(fact[2*i],1);
					mpz_set_ui(fact[2*i+1],0);
				}
			
	}


	/************  vérification ******************/

	// le tableau a 2*(i+1) éléments donc on met le i à jour
	if (mpz_cmp_ui(fact[0],1)!=0) i++; 		
	int j;
	//printf ("  p-1  a %d facteurs qui sont  : \n",i);
	mpz_t tmp, tmp3;
	int tmp2;
	mpz_init(tmp);
	mpz_init(tmp3);
	mpz_set_ui(tmp,1);

	for(j=0;j<i;j++)
	{
		gmp_printf(" %Zd avec multiplicité %Zd \n" , fact[2*j],fact[2*j+1]);
		tmp2 = mpz_get_ui(fact[2*j+1]);
		mpz_pow_ui(tmp3,fact[2*j],tmp2);		// on élève chq facteur à la puissance qui convient
		mpz_mul(tmp,tmp,tmp3);
	}
	
	mpz_set_str(p_1,"93829103457701351047678150545874455874415200013145296892799092824465729890236733039644041875602059462560643180671763707936371771172424325634235538844763279109070226245563765394479104625421576563590785443467192831069583054627852233430288407576160818322777666044553172636695128063366168185951512357172769367900",10);
	//if (mpz_cmp(tmp,p_1)==0) gmp_printf(" c'est  bon \n");



	/************ L'Histoire d'Alice et Bob ***********/
	
	mpz_t p,gene, y_Alice,y_Bob,ex_Alice;
	mpz_inits(p,gene,y_Alice,y_Bob,ex_Alice,NULL);
	mpz_set_str(p,"93829103457701351047678150545874455874415200013145296892799092824465729890236733039644041875602059462560643180671763707936371771172424325634235538844763279109070226245563765394479104625421576563590785443467192831069583054627852233430288407576160818322777666044553172636695128063366168185951512357172769367901",10);
	mpz_set_str(gene,"92191477025440111726923875183538607191632201592088010786683684660712544433485839734552882570141978547914319891956607297274759593396083685267217314624115920729471288043472055672572743613774043169714657806752596729139979180925465108475772533639928148511100860491285555856554901597826620477860782633975817418049",10);
	mpz_set_str(y_Alice,"38981678469915864463459170512652786604902813717618386306528422388688531343589071644026761072484487191523441968773332572940728118810803616661131890277161228582061817925693294540745748766603854391531221451896311310433829505177134753306116534950915000974209271489260904268077018384864727303802539174256609098090",10);
	mpz_set_str(y_Bob,"49719309207600533969024809590905088229745960359873961211001626233784769809008286790931420271863316318751970715966431812374282172185267651027944366680719952779982427965152011258747865347947162909915005549093542247327203894659997504116571213929538952549650839258280681814451584925120598894819564449439844998046",10);	

	
	// on détermine l'exposant privé d'Alice
	Pohlig_Hellmann(p,gene,y_Alice,rstate,fact,2*i,& ex_Alice);
	mpz_powm(tmp,gene,ex_Alice,p);
	if (mpz_cmp(tmp,y_Alice)!=0)
	{
		printf(" Erreur lors du processus \n");
		return 1;
	}

	mpz_powm(tmp,y_Bob,ex_Alice,p);		// secret : y_Bob^exp_Alice
	gmp_printf("\n\n\n***************************** \n\nLe secret partagé par Alice et Bob est\n %Zd\n", tmp);
	#if 0
mpz_t f[4];
	for (i=0;i<4;i++)
	{
		mpz_init(f[i]);
	}

	mpz_set_ui(f[0],2);
	mpz_set_ui(f[1],2);
	mpz_set_ui(f[2],3);
	mpz_set_ui(f[3],2);
	mpz_set_ui(p_1,36);
gmp_printf( "hhi  %Zd \n",tmp);
#endif
	return 0;
}





void Euristique_Rho(gmp_randstate_t rstate)
{
	mpz_t p,g,j,x,odr,exposant;
	mpz_inits(p,g,j,x,odr,exposant,NULL);
	int tps =0;//u_1,u_2;

	while(tps< 3600*24)	// tant que le temps d'exécution est moins d'1 jour
	{

		// use first fonction to generate p et g

		// p est next prime 2^i et le next de ce prime
		tps=rho(p,g,odr,j,rstate,& exposant,1);
	}
	



	
}

void Euristique_Pohlig(gmp_randstate_t rstate)

{

	
}


/****************************************************************************
*								Main										*
***************************************************************************/
/*int main ()
{

	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,time(NULL));
	test1(rstate);
	
	test_Pohlig_Hellman(rstate);

	return 0;
}
*/
