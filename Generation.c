#include <stdio.h>
#include <stdlib.h>         
#include <time.h>

#include "GJB_Extraction.h"

/* Ce fichier contient les fonctions necessaires à la génération de S et à l'affichage */


/* 
Attention rand() renvoie un int donc n<= 32. A modifier, si on veut travailler avec n>32.
pt : Pointeur vers notre vecteur.
length : taille de ce dernier.
*/
void Generation_Vecteur(unsigned int *pt, unsigned int length){
    srand(time(NULL));
    for(int i=0; i<length; i++){
        pt[i] = rand(); /* faut voir comment RAND_MAX change pour voir se qu'il faut faire pour des unsigned long long, OU FAIRE NOTRE PROPRE RANDOM en X bits */
    }
}

/*
Affiche un entier 32 bits en binaire
*/
void Affiche_Vecteur(unsigned int n) {
    int bit;
    for(int i=0; i<32; i++){
        bit = n&1;
        printf("%d", bit);
        n = n>>1;
    }
    printf("\n");
}
/*
Affiche un ensemble d'entiers 32 bits en binaire
*/
void Affiche_Ensemble(unsigned int *pt, unsigned int length){
    printf("{");
    for(int i=0; i<length; i++){
        printf("{");
        Affiche_Vecteur(pt[i]);
        printf("}");
    }
    printf("}");
}
