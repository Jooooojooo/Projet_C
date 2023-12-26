#include <stdio.h>
#include <stdlib.h>         
#include <time.h>


/* Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S*/

/*
Version Facile de chi_a_S sans trier S
*/

int chi_a_S(int a, int* S, int** pt , int Taille){  /* Prend en entrée a, S, un vecteur à modifier, et le nombre d'éléments de S */  
    int k = 0;  /* Initialisation de k entier, qui sera la taille de l'ensemble \chi_a(S) */
    for(int i = 0; i<Taille; i++){      /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
        if (((S[i]) > a) & (((S[i])^a)> S[i])){ /* Si la 1ère condition d'appartenance à S est vérifiée */
            for(int j=0; j<Taille; j++){        /* Initialisation de la 2ème for loop qui itère sur le nombre d'éléments de S */
                if (((S[i])^a) == S[j]){        /* Si x XOR a est égal à un élément arbitraire de S */
                    pt[0][k] = S[i];            /* Alors on définit x comme le k-ième élément de pt[0] (== "x est dans chi_a(S)") */
                    k++;                        /* Augmente la valeur de k de 1 pour continuer l'itération */
                    if(k%1024 == 0){            /* Si k est un multiple de 1024 */
                        pt[0] = (int*) realloc(pt[0],(k+1024)*sizeof(int)); /* ??? */
                    } 
                }
            }
        }
    }
    return(k);  /* Après avoir modifié le vecteur d'entrée, la fonction retourne k, la longueur de ce vecteur*/
}

int GJBExtraction(int* S, int d, int Taille,int**pt, int*** L){ /* Prend en entrée S, d, le nombre d'éléments de S, un tableau 
(à utiliser comme entrée pour Chi_a_S) et un autre tableau vide qui sera modifié par la fonction (la liste L)*/
    int Taille_L_new=0;     /* Initialisation de la taille de L_new, la liste L post-modification */
    int Taille_L_old=0;     /* Initialisation de la taille de L_old, la liste L pre-modification */
    int Taille_S_a;         /* Initialisation de la taille de \Chi_a(S) */
    if((d-1)!=0){           /* Si d-1 ≠ 0 */
        for (int i=0; i<Taille; i++){   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
        Taille_S_a= chi_a_S(S[i],S,pt,Taille);  /* Taille_S_a est définie comme la taille de \Chi_a(S) où a est le i-ème
        élément de S. La fonction modifie également pt, qui devient la liste des éléments de \Chi_a(S). */
            if (Taille_S_a> ((1<<(d-1))-1)){    /* Si la taille de \Chi_a(S) est au moins 2^{d-1} */
                Taille_L_old=GJBExtraction(pt[0], d-1,Taille_S_a,pt,L); /* Taille_L_old est définie comme la taille de 
                la liste L', la liste des GJB obtenues à partir de \Chi_a(S) de taille d-1 */
                /* ERREUR POTENTIELLE: ICI Taille_L_new EST TOUJOURS 0 DONC LA FOR LOOP NE SE LANCE JAMAIS!!! */
                for(int j=Taille_L_new; j<Taille_L_new+Taille_L_old; j++){  /* Initialisation de la 1ère for loop 
                qui itère sur la position des nouveaux éléments ajoutés à L à cette étape (EN THÉORIE, CF ERREUR PRÉCÉDENTE) */
                    if(j%1024 == 0){    /* Si l'itérante j est un multiple de 1024 */
                        L[0]=(int **)realloc(L[0],(j+1024)*sizeof(int*)); /* ERREUR ICI */
                        for(int a=0; a<1024; a++){      /* On alloue de la mémoire pour 1024 cases supplémentaires pour L*/
                            L[0][j+a]=malloc(d*sizeof(int));
                        }}
                        L[0][j][d-1]=S[i];  /* L'on définit la j-ème coordonnée de L[0] comme la "a" d'origine, soit
                        le i-ème élément de S, pour une base de taille d-1 */
                    
                Taille_L_new=Taille_L_new+Taille_L_old;     /* PROBABLEMENT LIÉE À L'ERREUR D'ORIGINE - À VÉRIFIER */
                }
            }
        }
    }
    /* Le cas où d-1 ≠ 0 n'est plus vérifié, c'est à dire le cas où d = 1*/
    for(int l=0;l<Taille;l++){      /* Initialisation de la for loop qui itère sur le nombre d'éléments de S */
        if((l)%1024 ==0){           /* Méthode de réalloument d'espace pour L comme dans le cas précédent, en blocs de 1024 bits*/
         L[0]=(int **)realloc(L,(l+1024)*sizeof(int*));
            for(int q=0;q<1024;q++){
                L[0][l]=malloc(d*sizeof(int));
            }
        }
        L[0][l][0]=S[l];}           /* Définir le l-ème élément de L[0] comme le l-ème élément de S, vu qu'en dimension 1 
        tout élément NON-NUL de S forme la base d'un sous-espace de dimension 1*/
        /* ATTENTION AU CAS NON-NUL, ICI L'ON PEUT RAJOUTER 0 PAR INADVERTANCE - À CORRIGER */
    Taille_L_new=Taille_L_new + Taille;     /* Rajoute la taille de S à la taille de L quand d-1 ≠ 0, afin de bien 
    avoir la taille finale de L; ATTENTION AU CAS OÙ 0 EST DANS S, DANS CE CAS-LÀ IL FAUT FAIRE 
    Taille_L_new = Taille_L_new + Taille - 1 POUR AVOIR LA VRAIE TAILLE */
    return(Taille_L_new);}      /* Après avoir modifié L, la fonction retourne la taille finale de L*/


/* Version où l'affichage est inversé par rapport à la lecture des bits*/  

void Affiche_Vecteur(int v, int n) {    /* Fonction d'affichage standard de int/long int en base 2 */
    int bit;
    for(int i=0; i<n; i++){
        bit = v&1;
        printf("%d", bit);
        v = v>>1;
    }
}

void Affiche_Ensemble(int *pt, int length, int n){  /* Fonction d'affichage standard d'un ensemble de nombres en base 2 */
    for(int i=0; i<length; i++){
        printf("{");
        Affiche_Vecteur(pt[i],n);
        printf("}");
  
    }

    if(length ==0){
        printf("{}\n");
    }
    else {
        printf("\n");
}
}
    
/* Version où l'affichage est dans le bon ordre par rapport à la lecture des bits*/
/* ATTENTION: possède _ord dans le nom mais n'ordonne pas (encore) l'ensemble donné en entrée */
void Affiche_Vecteur_ord(int v, int n) {    /* Fonction d'affichage standard de int/long int en base 2 à l'envers*/
    int bit;
    for(int i=0; i<n; i++){
        bit = (v>>(n-1-i))&1;
        printf("%d", bit);
    }
}


void Affiche_Ensemble_ord(int *pt, int length, int n){  /* Fonction d'affichage standard d'un ensemble de nombres en base 2 à l'envers*/
    for(int i=0; i<length; i++){
        printf("{");
        Affiche_Vecteur_ord(pt[i],n);
        printf("}");
  
    }

    if(length ==0){
        printf("{}\n");
    }
    else {
        printf("\n");
}
}




int main(){
    /*Test pour n=3, d=2, S={0,1,2,3,5,7}, |S|=6, le résultat attendu est 4. 
    Chi_a_S est vide pour a dans{0,2,3,4,5}  et possède deux éléments pour a=1.*/
    int* S=malloc(6*sizeof(int));
    S[0]=0;
    S[1]=1;
    S[2]=2;
    S[3]=3;
    S[4]=6;
    S[5]=7;

    int n=3;
    int d=2;
    int Taille = 6;

    int*** L=malloc(sizeof(int**));
    L[0]=malloc(1024*sizeof(int*));    
    for(int j=0;j<1024;j++){
        L[0][j]=malloc(d*sizeof(int));
    }
    
    int** pt=malloc(sizeof(int*));
    pt[0]=malloc(1024*sizeof(int));
    
    
    int T;
    printf("S :");
    Affiche_Ensemble_ord(S,Taille,n);
    for(int a=0;a<6;a++){
    T=chi_a_S(a,S,pt,Taille);
    printf("chi_%d_S : ", a);
    Affiche_Ensemble_ord(pt[0],T,n);}


    int x=GJBExtraction(S,2,6,pt,L);
    for(int j=0;j<x;j++){
        printf("Le %d-ème sous espace de dimension %d est : ",j+1,d);
        Affiche_Ensemble_ord(L[0][j],d,n);
    }

    

    
    
}
