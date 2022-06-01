#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>


gmp_randstate_t rstate;

//Affichage
void aff(mpz_t* List, int taille,int type)
{
	int i;
	if(type==0)//affichage d'une liste d'entier classique
	{
		printf("la taille %i \n",taille);
		
		printf("( ");	
		for(i=0;i<taille-1;i++)
		{
			printf("(i =%d ) , ",i);
			gmp_printf(" %Zd \n ",List[i]);
		}
		if (taille!=0)
		{
			printf("(i =%d ) , ",taille-1);
			gmp_printf(" %Zd  )\n \n  ",List[taille-1]);
		}	
		else
			printf(")\n ");
	}
	else //affichage d'une liste de type:[facteur1,exposant1,facteur2,exposant2,...]
	{
		printf("( ");	
		for(i=0;i<taille-2;i=i+2)
		{
			printf("(i =%d ), ",i/2);
			gmp_printf(" %Zd , %Zd \n ",List[i],List[i+1]);
		}
		printf("(i =%d ), ",i/2);
		gmp_printf(" %Zd , %Zd ",List[i],List[i+1]);
		printf(")\n ");
	}
}
//******************************************************************************************************
int AlgoRS(mpz_t **Result,int taille,mpz_t m)
{
	/*Etant donné un tableau de grands entiers,sa taille et un grand entier m, on génère aléatoirement des  grands entiers que l'on stocke dans le tableau.
 La taille du tabeau étant variable, on renvoie sa taille à la fin.
	*/
	mpz_t* ListEntier;
	ListEntier=(mpz_t*)malloc(taille*sizeof(mpz_t));
	int newtaille=taille;//taille courante du tableau
	int k=0,i;
	mpz_init(ListEntier[0]);
	mpz_set(ListEntier[0],m);
	
	while(mpz_cmp_ui(ListEntier[k],1)!=0)//tant que l'on ne génère pas 1 on continue
	{
				
		k=k+1;
		if(k>=newtaille)//vérification pour savoir si on va dépasser la taille du tableau
		{
			ListEntier=(mpz_t*)realloc(ListEntier,(newtaille+newtaille/2)*sizeof(mpz_t));//on réalloue la taille du tableau
			newtaille= newtaille + newtaille/2;//on actualise la taille
				
		}		
		mpz_init(ListEntier[k]);
		mpz_urandomm(ListEntier[k],rstate,ListEntier[k-1]);//on génère n(k)
		mpz_add_ui(ListEntier[k],ListEntier[k],1);//on ajoute 1 à n(k) car la génèration se fait entre 0 et n(k-1)-1 et on veut entre 1 et n(k-1)
		
	}
	*Result=(mpz_t*)realloc(*Result,(k)*sizeof(mpz_t));//on copie les n(1),...,n(k) valeurs dans Result[]
	for(i=1;i<=k;i++)
	{
		mpz_init((*Result)[i-1]);
		mpz_set((*Result)[i-1],ListEntier[i]);
	}
	for(i=0;i<=k;i++)
	{
		mpz_clear(ListEntier[i]);
	}
	free(ListEntier);

	/*******************************************/
	//printf("\n *******True RS length %i*******\n", k);
	/********************************************/

	return k;
}
//********************************************************************************************************

int AlgoRFN(mpz_t **decomp,int taille,mpz_t m,mpz_t p)
{
	
	int i,l;//variable de boucle
	int k;//taille de la liste des premiers
	int t=0;
	mpz_t x;// random <=m
	mpz_init(x);

	mpz_t y;//produit de nombres premiers
	mpz_init(y);
	
	
	mpz_t *ListP;

	mpz_t *ListEntier;//sortie de AlgoRS: liste d'entiers
	ListEntier=(mpz_t*)malloc(taille*sizeof(mpz_t));
	while(1)
	{											
		t=AlgoRS(&ListEntier,taille,m);//en sortie on a (n(1) ,.., n(k)) dans ListEntier , et actualisation de taille
		/****************/
		//printf("aff Liste d'entiers (rfn)\n");
		/****************/
		//aff(ListEntier,t,0);		
		ListP=(mpz_t *)malloc(t*sizeof(mpz_t));// on alloue la mémoire du tableau des premiers
		k=0;
		for(i=0;i<t;i++)//on parcours la liste des entiers
		{
			if((mpz_probab_prime_p(ListEntier[i],5))!=0)//si n(i) est premier (ou probablement) on le place dans ListP
			{
				mpz_init(ListP[k]);
				mpz_set(ListP[k],ListEntier[i]);
				k++;
			}
		}
		ListP=(mpz_t*)realloc(ListP,(k)*sizeof(mpz_t));  //on réduit le tableau pour avoir que les p(0),...,p(k-1) valeurs
		/****************/
		//printf("aff List de premier (rfn)\n");
		/****************/
		//aff(ListP,k,0);
		mpz_set_ui(y,1);
		
		for(i=0;i<k;i++)//calcul de y
		{
			mpz_mul(y,y,ListP[i]);
			
		}
		//gmp_printf("m= %Zd \n",m);		
		//gmp_printf("y= %Zd \n",y);
		if((mpz_cmp(m,y))>=0)// y<=m ?
		{
			mpz_urandomm(x,rstate,m);  //on génère 1<= x <=m
			//gmp_printf("x= %Zd \n",x);
			if((mpz_cmp(y,x)) >= 0) // x<= y ?
			{
				int cpt=1;
				
				(*decomp)=(mpz_t*)realloc((*decomp),2*k*sizeof(mpz_t));
				//mpz_set(decomp[0],ListP[0]);				
				l=0;
				for(i=0;i<k-1;i++)//fabrication de la liste [facteur,exposant,...]
				{
					if(mpz_cmp(ListP[i],ListP[i+1])==0)
					{
						cpt++;
					}
					else
					{
						mpz_init((*decomp)[l]);
						mpz_init((*decomp)[l+1]);
						mpz_set((*decomp)[l],ListP[i]);
						mpz_set_ui((*decomp)[l+1],cpt);
						cpt=1;
						l=l+2;
					}
				}
				mpz_init((*decomp)[l]);
				mpz_init((*decomp)[l+1]);
				mpz_set((*decomp)[l],ListP[i]);
				mpz_set_ui((*decomp)[l+1],cpt);
				l=l+2;
				if(l<2*k)
					(*decomp)=(mpz_t*)realloc((*decomp),l*sizeof(mpz_t));
				mpz_set(p,y);
			 	return l;//taille de decomp
			}
		}
	}

}
//********************************************************************************************************
void gen(mpz_t p,mpz_t *ListP,int taille,mpz_t* g)
{
	int i;
	mpz_t var,beta,q,p_1;
	mpz_inits(var,beta,q,p_1,NULL);
	mpz_t *gamma_i;
	gamma_i=(mpz_t*)malloc((taille/2)*sizeof(mpz_t));
	//calcul des gamma_i=alpha^((p-1)/(q_i^e_i))
	for(i=0;i<taille/2;i++)
	{
		do{
			mpz_urandomm(var,rstate,p);
			if(mpz_cmp_ui(var,0)==0)
				mpz_add_ui(var,var,1);
			mpz_sub_ui(p_1,p,1);// calcul p-1
			mpz_divexact(q,p_1,ListP[i*2]);
			//mpz_cdiv_q(q,p_1,ListP[i*2]);//q = (p-1)/qi
			mpz_powm(beta,var,q,p);//beta = var^q
		}while(mpz_cmp_ui(beta,1)==0);

		mpz_init(gamma_i[i]);
		mpz_set(q,ListP[i*2]);//qi
		mpz_powm(q,q,ListP[i*2+1],p_1);//qi^ei
		mpz_divexact(q,p_1,q);//(p-1)/qi^ei
		//mpz_cdiv_q(q,p_1,q);
		mpz_powm(gamma_i[i],var,q,p);


		gmp_printf("gamma_i[%d]=%Zd\n",i,gamma_i[i]);
	}
	mpz_set_ui(*g,1);
	//gmp_printf("g =%Zd \n",*g);
	for(i=0;i<taille/2;i++)
	{
		mpz_mul(*g,*g,gamma_i[i]);
		mpz_mod(*g,*g,p);
	}
	//gmp_printf("g =%Zd \n",*g);
}
//**********************************************************************************************************
void gen_groupe (unsigned int k,mpz_t p, mpz_t* gamma)
{
	mpz_t m,q;
	mpz_inits(m,q,NULL);
	mpz_ui_pow_ui(m,10,k);
	
	int taille=15;
	mpz_t *ListP;
	ListP=(mpz_t*)malloc(sizeof(mpz_t));	
		
	do{
		//RFN + aff
		taille = AlgoRFN(&ListP,taille,m,p);
		//gmp_printf("p-1=%Zd \n",p);
		//printf("aff List de premier (gen_groupe)\n");
		//aff(ListP,taille,1);
	
		/*
		for(i=0;i<taille-1;i=i+2)
		{
			mpz_powm(q,ListP[i],ListP[i+1],p);
			mpz_mul(p,p,q);
		}*/
		mpz_add_ui(p,p,1);
		
	}while(mpz_probab_prime_p(p,5)==0);
	printf("aff List de premier (gen_groupe)\n");
	aff(ListP,taille,1);
	gen(p,ListP,taille,gamma);
	
}

void gen_groupe2 (unsigned int k,mpz_t p, mpz_t* gamma,mpz_t *ListP,int *sizeP)
{
	mpz_t m,q;
	mpz_inits(m,q,NULL);
	mpz_ui_pow_ui(m,10,k);
	
	int taille=*sizeP;
	/*mpz_t *ListP;
	ListP=(mpz_t*)malloc(sizeof(mpz_t));	
	*/	
	do{
		//RFN + aff
		taille = AlgoRFN(&ListP,taille,m,p);
		gmp_printf("p-1=%Zd \n",p);
		printf("aff List de premier (gen_groupe)\n");
		aff(ListP,taille,1);
	
		/*
		for(i=0;i<taille-1;i=i+2)
		{
			mpz_powm(q,ListP[i],ListP[i+1],p);
			mpz_mul(p,p,q);
		}*/
		mpz_add_ui(p,p,1);
		
	}while(mpz_probab_prime_p(p,5)==0);
	*sizeP=taille;
	gen(p,ListP,taille,gamma);
	
}

//*********************Baby step Giant step**************************************************************
int babystep_giantstep(mpz_t p,mpz_t alpha,mpz_t beta,mpz_t result)
{
	mpz_t m,rest,p_1;
	mpz_inits(m,rest,p_1,NULL);
	//mpz_sub_ui(p_1,p,1);
	mpz_sqrtrem(m,rest,p);
	if(mpz_cmp_ui(rest,0)!=0)
		mpz_add_ui(m,m,1);

	mpz_t *Elm;
	unsigned long int len;
	len=mpz_get_ui(m);
	Elm=(mpz_t*)malloc(len*sizeof(mpz_t));

	//printf("ok\n");

	int i,j;
	for(i=0;i<len;i++)//création de la "table" avec les alpha^j
	{
		mpz_init(Elm[i]);
		mpz_powm_ui(Elm[i],alpha,i,p);
		//gmp_printf("Elm[%d]=%Zd \n",i,Elm[i]);
	}
	mpz_t x,gamma;
	mpz_inits(x,gamma,NULL);
	/*mpz_powm(x,alpha,m,p);
	mpz_t inv;
	mpz_init(inv);	
	if(mpz_invert(inv,x,p)==0)	
		printf("L'inverse n'existe pas\n");
	*/
	mpz_t op;
	mpz_init(op);
	//gmp_printf("alpha = %Zd \n",alpha );
	if(mpz_invert(op,alpha,p)==0)
		gmp_printf("erreur \n");
	//mpz_neg(op,m);//op = -m
	mpz_powm(x,op,m,p);
	mpz_set(gamma,beta);
	
	//printf("ok\n");

	for(i=0;i<len;i++)
	{
		//gmp_printf("gamma=%Zd \n",gamma);
		for(j=0;j<len;j++)//on cherche si gamma est dans la table
		{
			if(mpz_cmp(Elm[j],gamma)==0)
			{
				mpz_mul_ui(result,m,i);
				mpz_add_ui(result,result,j);
				return 0;
			}
		}
		mpz_mul(gamma,gamma,x);//on augmente gamma
		mpz_mod(gamma,gamma,p);
	}
	return 1;
}

//**********************************************************************************************************
void verif(mpz_t alpha,mpz_t exp,mpz_t beta,mpz_t p)//vérifi si alpha^exp = beta mod p
{
	mpz_t ver;
	mpz_init(ver);
	mpz_powm(ver,alpha,exp,p);
	gmp_printf("verif = %Zd \n",ver);
	if(mpz_cmp(ver,beta)==0)
		printf("On a trouvé le log de beta !\n");
	else
		printf("Erreur le programme n'est pas bon =(\n");

}
