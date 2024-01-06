/* Inclusion des bibliothèques standard de C */
#include <stdio.h>
#include <stdlib.h>         
#include <time.h>

/* Variable globale nécessaire à GJBExtraction - à éviter si possible */
int compteur_global=0;

/* Déclaration préalable des fonctions d'affichage de résultats */
void Affiche_Vecteur(int v, int n);
void Affiche_Ensemble(int *pt, int length, int n);
void Affiche_Ensemble_ord(int *pt, int length, int n);
void Affiche_Vecteur_ord(int v, int n);

/* 
Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S
VERSION NAIVE: PAS DE TRI DES ELEMENTS DE S, FONCTION PHI BLOQUÉE À L'IDENTITÉ, 
REPRÉSENTATION FIDÈLE DE L'ALGO DONNÉ EN PSEUDO-CODE.
*/


/* Fonctions d'affichage de vecteurs */

/* Version où l'affichage est inversé par rapport à la lecture des bits*/  
void Affiche_Vecteur(int v, int n)
{
    int bit;
    for(int i=0; i<n; i++)
    {
        bit = v&1;
        printf("%d", bit);
        v = v>>1;
    }
}

void Affiche_Ensemble(int *pt, int length, int n)
{ 
    for(int i=0; i<length; i++)
    {
        printf("{");
        Affiche_Vecteur(pt[i],n);
        printf("}");
  
    }

    if(length ==0)
    {
        printf("{}\n");
    }
    else
    {
        printf("\n");
    }
}

/* Version où l'affichage est dans le bon ordre par rapport à la lecture des bits*/
void Affiche_Vecteur_ord(int v, int n)
{
    int bit;
    for(int i=0; i<n; i++)
    {
        bit = (v>>(n-1-i))&1;
        printf("%d", bit);
    }
}

void Affiche_Ensemble_ord(int *pt, int length, int n)
{ 
    for(int i=0; i<length; i++)
    {
        printf("{");
        Affiche_Vecteur_ord(pt[i],n);
        printf("}");
  
    }

    if(length ==0)
    {
        printf("{}\n");
    }
    else
    {
        printf("\n");
    }
}



/*
Version Facile de chi_a(S), sans altérer S avant le passage dans la fonction
*/

int chi_a_S(int a, int* S, int** pt , int Taille){
    /*
    ENTRÉE: 
    - un élément a de F_2^n
    - un sous-ensemble S de F_2^n, préalablement défini
    - un vecteur arbitraire pt, mais avec pt[0] préalablement malloqué à 1024*sizeof(int*)
    - un entier Taille, le nombre d'éléments de S, préalablement connu ou calculé

    SORTIE:
    - Un entier, le nombre d'élements de l'ensemble chi_a(S)

    FAIT ÉGALEMENT:
    - modifie les éléments de pt[0] 1 par 1 en les éléments de chi_a(S)
    - réalloque de la mémoire à pt[0] par blocs de 1024 si la mémoire précédente est complètement utilisée
    */

    int k = 0;      /* Initialisation de k entier, qui sera la valeur de sortie */
    for(int i = 0; i<Taille; i++)   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
    {
        if (((S[i]) > a) & (((S[i])^a)> S[i]))  /* Si la 1ère condition d'appartenance à S est vérifiée */
        {
            for(int j=0; j<Taille; j++) /* Initialisation de la 2ème for loop qui itère sur le nombre d'éléments de S */
            {
                if (((S[i])^a) == S[j]) /* Si x XOR a est égal à un élément arbitraire de S (== appartient à S) */
                {
                    pt[0][k] = S[i];    /* Alors on définit x comme le k-ième élément de pt[0] (== "x est dans chi_a(S)") */
                    k++;                /* Augmente la valeur de k de 1 pour continuer l'itération */
                    if(k%1024 == 0)     /* Si k est un multiple de 1024, réallouer de la mémoire supplémentaire à pt[0] */
                    {
                        pt[0] = (int*) realloc(pt[0],(k+1024)*sizeof(int));
                    } 
                }
            }
        }
    }
    return(k);  /* Après avoir fini de modifier pt[0], la fonction retourne k, la longueur de ce vecteur */
}

int GJBExtraction(int* S, int d, int Taille, int**pt, int*** L){
    /*
    ENTRÉE: 
    - un sous-ensemble S de F_2^n, préalablement défini
    - un entier d, avec 0 < d <= n
    - un entier Taille, le nombre d'éléments de S, préalablement connu ou calculé
    - un vecteur arbitraire pt, mais avec pt[0] préalablement malloqué à 1024*sizeof(int*)
    - une matrice arbitraire L, mais avec un malloc préalablement établi à:
        L[0]=malloc(1024*sizeof(int*));    
        for(int j=0;j<1024;j++){
            L[0][j]=malloc(d*sizeof(int));
        }

    SORTIE:
    - Un entier, le nombre de bases de Gauss-Jordan distinctes de S de dimension d

    FAIT ÉGALEMENT:
    - modifie les éléments de L[0] afin que chaque "ligne" L[0][i] corresponde à une base de GJ
    - réalloque de la mémoire à L[0] par blocs de 1024 si la mémoire précédente est complètement utilisée
    (== rajoute de la mémoire s'il y a besoin de plus de 'lignes').
    */
    int Taille_L_new=0; /* Initialisation de la taille de L_new, la liste L post-modification */
    int Taille_L_old=0; /* Initialisation de la taille de L_old, la liste L pre-modification */
    int Taille_S_a;     /* Initialisation de la taille de \Chi_a(S) */
    if((d-1)!=0){       /* Si d ≠ 1 */
        for (int i=0; i<Taille; i++){   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
        Taille_S_a= chi_a_S(S[i],S,pt,Taille);
        /*
        Taille_S_a est définie comme la taille de \Chi_a(S) où a est le i-ème élément de S. 
        La fonction modifie également pt, qui devient la liste des éléments de \Chi_a(S).
        */
            if (Taille_S_a>=((1<<(d-1))-1)){    /* Si \Chi_a(S) a au moins 2^{d-1} éléments */
                Taille_L_old=GJBExtraction(pt[0], d-1,Taille_S_a,pt,L)+Taille_L_new;
                /* ÉTAPE DE RÉCURSIVITÉ:
                Taille_L_old est définie comme la taille de la liste L' 
                (la liste des GJB de dimension d-1 obtenues à partir de \Chi_a(S)),
                ajoutée à Taille_L_new, la taille de L' à la précédente étape de récursion.
                */
                for(int j=Taille_L_new; j<Taille_L_old; j++){   /* Pour tous les nouveaux indices j ajoutés */
                    if(((j+1) % 1024 == 0)&((j+1)<Taille_L_old)){   /* Réallocation de mémoire */
                        L[0]=(int **)realloc(L[0],(j+1+1024)*sizeof(int*)); 
                        for(int a=0; a<1024; a++){
                            L[0][j+1+a]=malloc(d*sizeof(int));
                        }}

                        L[0][j][d-1]=S[i];}     /* La position (j,d-1) de L est changée en le i-ème élément de S,
                        qui par construction est un élément de la base de GJ correspondante */
                        Taille_L_new=Taille_L_old;  /* Taille_L_new est actualisée pour refléter l'ajout d'élements */
                               }}
                               return(Taille_L_new);}   /* retourne le nombre de "lignes" total à ce stade */

    if(d ==1){      /* Dans le cas où d=1 (condition d'arrêt de récurrence, vu qu'à chaque appel, d est réduit de 1) */
        for(int l=compteur_global;l<(Taille+compteur_global);l++){  /* Pour tous les nouveaux indices l ajoutés */
        if(((l+1)%1024 ==0) & ((l+1)<Taille)){      /* Réallocation de mémoire */
            L[0]=(int **)realloc(L[0],(l+1+1024)*sizeof(int*));
            for(int q=0;q<1024;q++){
                L[0][l+1+q]=malloc(d*sizeof(int));
            }
        }
        L[0][l][0]=S[l-compteur_global];}   /* Le premier élément de la l-ème ligne de L[0] est changé en le l-ème
        élément de S. En clair, quand d=1, tout élément forme une base de GJ du sous-espace de dimension 1 qu'il génère,
        et donc tout élément de S doit être ajouté */
        Taille_L_new=Taille;    /* Taille_L_new est actualisée pour refléter l'ajout d'élements */
        compteur_global += Taille;  /* compteur_global est actualisé pour refléter l'ajout d'élements */
        return(Taille_L_new);}  /* retourne le nombre de "lignes" ajoutées */
    }
/*
EXEMPLE:

S = {001, 010, 100, 110, 111} = {1,2,4,6,7}
CLAIREMENT DEUX SOUS-ESPACES DE DIMENSION 2, {0, 010, 100, 110} ET {0, 001, 110, 111}, DONC DEUX BASES
ON DEVRAIT AVOIR:
chi_001(S) = {110}, chi_010(S) = {100}, autres chi_a vides
DONC pour a=001 le code doit faire:
    L' <- GJBExtraction(chi_001(S), d=1) = GJBExtraction(110, 1) = {110}
    et ajouter {001}u{110} A L
puis pour a=010 faire:
    L' <- GJBExtraction(chi_010(S), d=1) = GJBExtraction(100, 1) = {100}
    et ajouter {010}u{100} A L
ET DONC OBTENIR L = {{001}{110}; {010}{100}}

*/




int main(int argc, char const *argv[]){
    /* POUR L'INSTANT: pas de contraintes sur argc, et on s'attend à ce que argv[1] = d */
    /* AMÉLIORATION POTENTIELLE: pouvoir mettre S en argument (mais inutile pour plus gros exemples) */

    /* Création de S */
    int* S=malloc(9*sizeof(int));
    S[0]=1;
    S[1]=2;
    S[2]=4;
    S[3]=6;
    S[4]=7;
    S[5]=8;
    S[6]=9;
    S[7]=11;
    S[8]=13;

    /* DÉFINITION DES VALEURS DE n, d, Taille (= nombre d'éléments de S) */
    int n=4;
    int d=atoi(argv[1]);
    int Taille = 9;

    /* Allocation initiale de mémoire pour L[0] */
    int*** L=malloc(sizeof(int**));
    L[0]=malloc(1024*sizeof(int*));    
    for(int j=0;j<1024;j++){
        L[0][j]=malloc(d*sizeof(int));
    }
    
    /* Allocation initiale de mémoire pour pt[0] */
    int** pt=malloc(sizeof(int*));
    pt[0]=malloc(1024*sizeof(int));
    
    /* Affiche l'ensemble S "proprement"*/
    printf("S :");
    Affiche_Ensemble_ord(S,Taille,n);

    /*
    ANCIEN CODE DE TEST DE LA FONCTION CHI - À SUPPRIMER UNE FOIS CERTAIN DE LA VALIDITÉ

    int T;
    int J;
    for(int a=0;a<Taille2;a++)
    {
        T=chi_a_S(S[a],S,pt,Taille2);
        printf("chi_%d_S : ", S[a]);

        Affiche_Ensemble_ord(pt[0],T,n);
    }
    T=chi_a_S(1,S,pt,Taille2);
    printf("chi_%d_S : ", 1);
    Affiche_Ensemble_ord(pt[0],T,n);
    printf("S : ");
    Affiche_Ensemble_ord(pt[1],Taille2,n);

    J=chi_a_S(1,pt[0],pt,T);
    printf("chi_1_%d_S : ",1);
    Affiche_Ensemble_ord(pt[0],J,n);
    printf("S : ");
    Affiche_Ensemble_ord(pt[1],T,n);
    */
    
    /* x est défini comme le nombre total de bases de GJ de S de dimension d 
    Cette étape modifie également pt de nombreuses fois (sans toutefois garder un intérêt à la fin), 
    et surtout définit L comme la matrice de taille x fois d, où les lignes sont les bases distinctes de GJ */
    int x=GJBExtraction(S,d,Taille,pt,L);

    /* Affichage de x et de toutes les bases de GJ calculées à l'étape précédente */
    printf("Le nombre de GJB vaut %d\n",x);
    for(int j=0;j<x;j++){
        printf("Le %d-ème sous espace de dimension %d est : ",j+1,d);
        Affiche_Ensemble_ord(L[0][j],d,n);
    }
    
}
