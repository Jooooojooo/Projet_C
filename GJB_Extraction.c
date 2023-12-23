#include <stdio.h>
#include <stdlib.h>         
#include <time.h>


/* Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S*/

/*
Version Facile de chi_a_S sans trier S
*/

int chi_a_S(int a, int* S, int** pt , int Taille){  
    int k = 0;
    for(int i = 0; i<Taille; i++){
        if (((S[i]) > a) & (((S[i])^a)> S[i])){
            for(int j=0; j<Taille; j++){
                if (((S[i])^a) == S[j]){
                    pt[0][k] = S[i];
                    k++;
                    if(k%1024 == 0){
                        pt[0] = (int*) realloc(pt[0],(k+1024)*sizeof(int));
                    } 
                }
            }
        }
    }
    return(k);
}

int GJBExtraction(int* S, int d, int Taille,int**pt, int*** L){
    int Taille_L_new=0;
    int Taille_L_old=0;
    int Taille_S_a;
    if((d-1)!=0){
        for (int i=0; i<Taille; i++){
        Taille_S_a= chi_a_S(S[i],S,pt,Taille);
            if (Taille_S_a> ((1<<(d-1))-1)){
                Taille_L_old=GJBExtraction(pt[0], d-1,Taille_S_a,pt,L);
                for(int j=Taille_L_new; j<Taille_L_new+Taille_L_old; j++){
                    if(j%1024 == 0){
                        L[0]=(int **)realloc(L[0],(j+1024)*sizeof(int*));
                        for(int a=0; a<1024; a++){
                            L[0][j+a]=malloc(d*sizeof(int));
                        }
                        L[0][j][d-1]=S[i];
                    }
                Taille_L_new=Taille_L_new+Taille_L_old;    
                }
            }
        }
    }
    for(int l=0;l<Taille;l++){
        if((l)%1024 ==0){
         L[0]=(int **)realloc(L,(l+1024)*sizeof(int*));
            for(int q=0;q<1024;q++){
                L[0][l]=malloc(d*sizeof(int));
            }
        }
        L[0][l][0]=S[l];}
    
    Taille_L_new=Taille_L_new + Taille;
    return(Taille_L_new);}


/* Version où l'affichage est inversé par rapport à la lecture des bits*/  

void Affiche_Vecteur(int v, int n) {
    int bit;
    for(int i=0; i<n; i++){
        bit = v&1;
        printf("%d", bit);
        v = v>>1;
    }
}

void Affiche_Ensemble(int *pt, int length, int n){ 
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
void Affiche_Vecteur_ord(int v, int n) {
    int bit;
    for(int i=0; i<n; i++){
        bit = (v>>(n-1-i))&1;
        printf("%d", bit);
    }
}


void Affiche_Ensemble_ord(int *pt, int length, int n){ 
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
    for(int j=0;j<1024;d++){
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

    int y=x%1024;

    free(pt[0]);
    free(pt);


    if(y =!0){
        for(int k=0;k<(x+(1024-(x%1024)));k++){
        free(L[0][k]);}}
    if(y==0){
        for(int k=0; k<x;k++ )
        { free(L[0][k]);
        }
    }
    free(L[0]);
    free(L);
    
}
