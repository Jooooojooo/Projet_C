#include <stdio.h>
#include <stdlib.h>         
#include <time.h>

#include "Generation.h"

/* Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S*/

/*
Version Facile de chi_a_S sans trier S
*/
int chi_a_S(int a, int* S, int* pt, int nombre_de_vect){
    int k = 0;
    for(int i = 0; i<nombre_de_vect; i++){
        if (((S[i]) > a) & (((S[i])^a)> S[i])){
            for(int j=0; j<nombre_de_vect; j++){
                if (((S[i])^a) == S[j]){
                    pt[k] = S[i];
                    k++;
                    if(k%1024 == 0){
                        pt = (int*) realloc(pt,(k+1024)*sizeof(int));
                    } 
                }
            }
        }
    }
    return(k);
}

/*
newL est la pour stocke la nouvelle adresse de L en sortie de fonction car realloc peut modifier l'addresse de.
*/
int GJBExtraction(int* S, int d, int nombre_de_vect,int*pt, int** L, int** newL){
    int Taille_L_new=0;
    int Taille_L_old=0;
    int Taille_S_a;
    if((d-1)!=0){
        for (int i=0; i<nombre_de_vect; i++){
        Taille_S_a= chi_a_S(S[i],S,pt,nombre_de_vect);
            if (Taille_S_a> ((1<<(d-1))-1)){
                Taille_L_old=GJBExtraction(pt, d-1,Taille_S_a,pt,L,newL);
                for(int j=Taille_L_new; j<Taille_L_new+Taille_L_old; j++){
                    if(j%1024 == 0){
                        L=(int **)realloc(L,(j+1024)*sizeof(int*));
                        for(int a=0; a<1024; a++){
                            L[j+a]=malloc(d*sizeof(int));
                        }
                        L[j][d-1]=S[i];
                    }
                Taille_L_new=Taille_L_new+Taille_L_old;    
                }
            }
        }
    }
    for(int l=0;l<nombre_de_vect;l++){
        if((l)%1024 ==0) L=(int **)realloc(L,(l+1024)*sizeof(int*));
            for(int q=0;q<1024;q++){
                L[l]=malloc(d*sizeof(int));
            }
        L[l][0]=S[l];
        Taille_L_new=Taille_L_new + nombre_de_vect;
    }
    newL= L;
    return(Taille_L_new);
}
   



int main(){
    int* S=malloc(6*sizeof(int));
    int** newL=malloc(1024*sizeof(int*));
    for(int i=0;i<1024;i++){
        newL[i]=malloc(2*sizeof(int));}
    int** L=malloc(1024*sizeof(int*));
    for(int i=0;i<1024;i++){
        L[i]=malloc(2*sizeof(int));
    }
    S[0]=0;
    S[1]=1;
    S[2]=2;
    S[3]=3;
    S[4]=6;
    S[5]=7;
    int *pt=malloc(1024*sizeof(int));
    int x=GJBExtraction(S,2,6,pt,L,newL);
    for(int j=0;j<x;j++){
        Affiche_Ensemble(newL[j],3);
    
    }

}