#include <stdio.h>          /* TOUT COMMENTER */
#include <stdlib.h>         
#include <time.h>

#include "GJB_Extraction.h"

/* Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S*/

/* Pour l"instant n<=32, donc le type sera int. Si on veut n>32, il faudra mofier en long long*/


int* Resultat; // Liste qui va stocker chi_a_S
int** L; // Liste qui contiendra tout les sous espacess de dimension d
int n; // dimension des vecteurs
int nombre_de_vect;


typedef struct Int_Matrice Int_Matrice;
struct Int_Matrice
{
    int x;
    int** y;
};


int chi_a_S(int a, int* S, int* Resultat)   /* Version Facile de chi_a_S sans trier S*/
{
    int k=0;
    for(int i=0; i<nombre_de_vect; i++)
    {
        if ((S[i] > a) & (S[i]^a > S[i]))
        {
            for(int j=0; j<nombre_de_vect; j++)
            {
                if (S[i]^a == S[j])
                {
                    
                    Resultat[k] = S[i];
                    k++;
                    if(k%1024 == 0){
                        Resultat=(int*) realloc(Resultat,(k+1024)*sizeof(int));
                    }
                }
            }
        }
    }
    return(k);
}

Int_Matrice GJBExtraction(int* S, int d) /* A finir.*/
{
    Int Matrice Sortie;
    int Taille_prime;
    int compteur;
    int Taille;
    int ** L_prime;
    for (int i=0; i<nombre_de_vect; i++)
    {
        Taille= chi_a_S(S[i], S, Resultat);
        if((d-1)!=0){
            if (Taille>= ((1<<(d-1))-1))
            {
                Taille_prime= GJBExtraction(Resultat, d-1).x; 
                L_prime= GJBExtraction(Resultat, d-1).y;
                compteur= compteur+Taille_prime;
                for(int j=0; j<Taille_prime; j++)
                {
                    if(j%1024 == 0){
                        L=(int **)realloc(L,(j+1024)*sizeof(long long));
                        for(int a=0; a<1024; a++){
                            L[j+a]=malloc(d*sizeof(long long));
                        }
                    L[j+i*Taille_prime][d]=S[i];
                    for( int p=0; p<d-1; p++){
                        L[j+i*Taille_prime][p]=L_prime[j][p];
                    }

                    }
                    
                }
                }
        }

    }
    Sortie.x=compteur;
    Sortie.y=L;
    return(Sortie);
}





    
