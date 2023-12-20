#include <stdio.h>          /* TOUT COMMENTER */
#include <stdlib.h>         
#include <time.h>

#include "Generation.h"
/* Ce fichier contient les fonctions necessaires à la génération de S et à l'affichage */


void Generation_Vecteur(unsigned int *pt, int l)  /* Attention rand() renvoie un int donc n<= 32. A modifier, si on veut travailler avec n>32.*/
{
    for(int i=0; i<l; i++)
    {
        pt[i]=rand()&((1<<n)-1);
    }
}

void Affiche_Vecteur(unsigned int y) /*Affiche un entier 32 bits en binaire*/
{
    int x;
    for(int i=0,i<n;i++){
        x=y&1;
        printf("%d",x);
        y=y>>1;
    }
    printf("\n");
}

void Affiche_Ensemble(unsigned int *pt, int l)
{
    for(int i=0;i<l;i++){
        Affiche_Vecteur(pt[i]);
    }

}
