/* Inclusion des bibliothèques standards de C */
#include <stdio.h>
#include <stdlib.h>         
#include <time.h>

/* Variable globale nécessaire à GJBExtraction - à éviter si possible */
int compteur_global = 0;

/* Déclaration préalable des fonctions d'affichage de résultats */
void Affiche_Vecteur(int v, int n);
void Affiche_Ensemble(uint64_t *pt, int length, int n);
void Affiche_Ensemble_ord(uint64_t *pt, uint64_t length, int n);
void Affiche_Vecteur_ord(uint64_t v, int n);

/* 
Ce fichier contiendra la fonction qui calcule les différents sous espace de dimension d d'un ensemble S
VERSION AMÉLIORÉE: TRI DES ELEMENTS DE S, FONCTION PHI IMPLÉMENTÉE
*/


/*
Fonctions d'affichage de vecteurs 
Version où l'affichage est inversé par rapport à la lecture des bits
*/  
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

void Affiche_Ensemble(uint64_t *pt, int length, int n)
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
void Affiche_Vecteur_ord(uint64_t v, int n){
    uint64_t i; 
    i = 1ULL<<(n-1);
    while(i>0){
        if(v&i)
            printf("1"); 
        else 
            printf("0"); 
        i >>= 1;
    }
}

void Affiche_Ensemble_ord(uint64_t *pt, uint64_t length, int n)
{ 
    for(uint64_t i=0; i<length; i++)
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

uint64_t rnd64(uint64_t n)
{
    const uint64_t z = 0x9FB21C651E98DF25;

    n ^= ((n << 49) | (n >> 15)) ^ ((n << 24) | (n >> 40));
    n *= z;
    n ^= n >> 35;
    n *= z;
    n ^= n >> 28;

    return n;
}


/* Fonction retournant le Most Significant Bit d'un élément a de F_2^n */
int MSB(uint64_t a, int n){
    int msb = 0;
    for(int k = 0; k<n; k++)
        if (((a>>k)&1) == 1)  msb = k;
    return msb;
}

void phi_d(uint64_t* S, int d, uint64_t* pt, int n, uint64_t Taille){
    for (uint64_t j = 0; j < Taille; j++)
    {
        int msb = MSB(S[j], n);
        int msb_liste[n]; 
        for(uint64_t k = 0; k<Taille; k++) msb_liste[MSB(S[k], n)]++;
        int i = 1;
        for (int k = msb+1; k < n; k++)
        {
            if (msb_liste[k]>(1<<i)) i++;        
        }
        if (i>d-1) pt[j] = 1;
        else pt[j] = 0;
    }
}


void permuter(uint64_t *a, uint64_t *b)
{
    uint64_t temp = *a;
    *a = *b;
    *b = temp;
}

void triRapid(uint64_t tab[], uint64_t first, uint64_t last) {
    uint64_t pivot, i, j;
    if(first < last) {
        pivot = first+1;
        i = first;
        j = last;
        while (i < j) {
            while(tab[i] <= tab[pivot] && i < last)
                i++;
            while(tab[j] > tab[pivot])
                j--;
            if(i < j) {
                permuter(&tab[i], &tab[j]);
            }
        }
        permuter(&tab[pivot], &tab[j]);
        triRapid(tab, first, j - 1);
        triRapid(tab, j + 1, last);
    }
}



/*
Version Facile de chi_a(S), sans altérer S avant le passage dans la fonction

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
uint64_t chi_a_S(uint64_t a, uint64_t* S, uint64_t** pt , uint64_t Taille){


    uint64_t k = 0;      /* Initialisation de k entier, qui sera la valeur de sortie */
    for(uint64_t i = 0; i<Taille; i++)   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
    {
        if (((S[i]) > a) & (((S[i])^a)> S[i]))  /* Si la 1ère condition d'appartenance à S est vérifiée */
        {
            for(uint64_t j=0; j<Taille; j++) /* Initialisation de la 2ème for loop qui itère sur le nombre d'éléments de S */
            {
                if (((S[i])^a) == S[j]) /* Si x XOR a est égal à un élément arbitraire de S (== appartient à S) */
                {
                    pt[0][k] = S[i];    /* Alors on définit x comme le k-ième élément de pt[0] (== "x est dans chi_a(S)") */
                    k++;                /* Augmente la valeur de k de 1 pour continuer l'itération */
                    if(k%1024 == 0)     /* Si k est un multiple de 1024, réallouer de la mémoire supplémentaire à pt[0] */
                    {
                        pt[0] = (uint64_t*) realloc(pt[0],(k+1024)*sizeof(uint64_t));
                    } 
                }
            }
        }
    }
    return(k);  /* Après avoir fini de modifier pt[0], la fonction retourne k, la longueur de ce vecteur */
}

/*
Version élaborée de chi_a(S), où l'on suppose que S est préalablement triée

ENTRÉE: 
- un élément a de F_2^n
- un sous-ensemble S de F_2^n, préalablement défini et trie.
- un vecteur arbitraire pt, mais avec pt[0] préalablement malloqué à 1024*sizeof(int*)
- un entier Taille, le nombre d'éléments de S, préalablement connu ou calculé

SORTIE:
- Un entier, le nombre d'élements de l'ensemble chi_a(S)

FAIT ÉGALEMENT:
- modifie les éléments de pt[0] 1 par 1 en les éléments de chi_a(S)
- réalloque de la mémoire à pt[0] par blocs de 1024 si la mémoire précédente est complètement utilisée
*/
uint64_t chi_a_S_Tri(uint64_t a, uint64_t* S, uint64_t** pt , uint64_t Taille){
    
    uint64_t m = 0;      /* Indice a partir du quel les elements de la liste sont plus grand que a */
    uint64_t k = 0;      /* Initialisation de k entier, qui sera la valeur de sortie */
    while ((S[m]) <= a)
    {
        m++;
    }
    for(uint64_t i = m; i<Taille; i++)   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S*/
    {
        if (((S[i])^a) > S[i])  /* Si la 1ère condition d'appartenance à S est vérifiée */
        {
            for(uint64_t j = i+1 ; j<Taille; j++) /* Initialisation de la 2ème for loop qui itère sur le nombre d'éléments de S */
            {
                if (((S[i])^a) == S[j]) /* Si x XOR a est égal à un élément arbitraire de S (== appartient à S) */
                {
                    pt[0][k] = S[i];    /* Alors on définit x comme le k-ième élément de pt[0] (== "x est dans chi_a(S)") */
                    k++;                /* Augmente la valeur de k de 1 pour continuer l'itération */
                    if(k%1024 == 0)     /* Si k est un multiple de 1024, réallouer de la mémoire supplémentaire à pt[0] */
                    {
                        pt[0] = (uint64_t*) realloc(pt[0],(k+1024)*sizeof(uint64_t));
                    } 
                }
            }
        }
    }
    return(k);  /* Après avoir fini de modifier pt[0], la fonction retourne k, la longueur de ce vecteur */
}

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
uint64_t GJBExtraction(uint64_t* S, int d, uint64_t Taille, uint64_t**pt, uint64_t*** L, int n, uint64_t* phi){

    uint64_t Taille_L_new=0; /* Initialisation de la taille de L_new, la liste L post-modification */
    uint64_t Taille_L_old=0; /* Initialisation de la taille de L_old, la liste L pre-modification */
    uint64_t Taille_S_a;     /* Initialisation de la taille de \Chi_a(S) */
    phi_d(S,d,phi,n,Taille);
    if((d-1)!=0){       /* Si d ≠ 1 */
        for (uint64_t i=0; i<Taille; i++){   /* Initialisation de la 1ère for loop qui itère sur le nombre d'éléments de S */
        if (phi[i]==1){
        Taille_S_a = chi_a_S(S[i],S,pt,Taille);
        /*
        Taille_S_a est définie comme la taille de \Chi_a(S) où a est le i-ème élément de S. 
        La fonction modifie également pt, qui devient la liste des éléments de \Chi_a(S).
        */
            if (Taille_S_a >= ((1<<(d-1))-1)){    /* Si \Chi_a(S) a au moins 2^{d-1} éléments */
                Taille_L_old=GJBExtraction(pt[0], d-1,Taille_S_a,pt,L,n,phi)+Taille_L_new;
                /* ÉTAPE DE RÉCURSIVITÉ:
                Taille_L_old est définie comme la taille de la liste L' 
                (la liste des GJB de dimension d-1 obtenues à partir de \Chi_a(S)),
                ajoutée à Taille_L_new, la taille de L' à la précédente étape de récursion.
                */
                for(uint64_t j=Taille_L_new; j<Taille_L_old; j++){   /* Pour tous les nouveaux indices j ajoutés */
                    if(((j+1) % 1024 == 0)&((j+1)<Taille_L_old)){   /* Réallocation de mémoire */
                        L[0]=(uint64_t **)realloc(L[0],(j+1+1024)*sizeof(uint64_t*)); 
                        for(int a=0; a<1024; a++){
                            L[0][j+1+a]=malloc(d*sizeof(uint64_t));
                        }}

                        L[0][j][d-1]=S[i];}     /* La position (j,d-1) de L est changée en le i-ème élément de S,
                        qui par construction est un élément de la base de GJ correspondante */
                        Taille_L_new=Taille_L_old;  /* Taille_L_new est actualisée pour refléter l'ajout d'élements */
                               }}}
                               return(Taille_L_new);}   /* retourne le nombre de "lignes" total à ce stade */

    if(d ==1){      /* Dans le cas où d=1 (condition d'arrêt de récurrence, vu qu'à chaque appel, d est réduit de 1) */
        for(int l=compteur_global;l<(Taille+compteur_global);l++){  /* Pour tous les nouveaux indices l ajoutés */
        if(((l+1)%1024 ==0) & ((l+1)<Taille)){      /* Réallocation de mémoire */
            L[0]=(uint64_t **)realloc(L[0],(l+1+1024)*sizeof(uint64_t*));
            for(int q=0;q<1024;q++){
                L[0][l+1+q]=malloc(d*sizeof(uint64_t));
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

    /* DÉFINITION DES VALEURS DE n ET d */
    int n=64;
    int d=atoi(argv[1]);

    /* Mise en place de l'horloge pour le temps de calcul */
    clock_t sec = clock();

    /* Création de S */
    uint64_t Taille = strtoull(argv[2], NULL, 10);
    uint64_t* S=malloc(Taille*sizeof(uint64_t));
    
    uint64_t state = 1;

    for (uint64_t i=0;i<Taille;i++)
    {
        const uint64_t n = rnd64(state++);
        S[i] = n;
    }

    /* Affichage du temps de calcul de S */
    double temps = ((double)(clock() - sec)/CLOCKS_PER_SEC);
    printf("Le temps de calcul de S est: %f secondes\n", temps);

    /* Allocation initiale de mémoire pour L[0] */
    uint64_t*** L=malloc(sizeof(uint64_t**));
    L[0]=malloc(1024*sizeof(uint64_t*));    
    for(int j=0;j<1024;j++){
        L[0][j]=malloc(d*sizeof(uint64_t));
    }
    
    /* Allocation initiale de mémoire pour pt[0] */
    uint64_t** pt=malloc(sizeof(uint64_t*));
    pt[0]=malloc(1024*sizeof(uint64_t));

    /* Allocation initiale de mémoire pour phi */
    uint64_t* phi=malloc(Taille*sizeof(uint64_t));
    
    /* Affiche l'ensemble S "proprement"
    printf("S :");
    Affiche_Ensemble_ord(S,Taille,n);
    
    printf("S en décimale :");
    for(uint64_t i=0; i<Taille; i++)
    {
        printf("{%llu}", S[i]);
    }
    printf("\n");
    */
    sec = clock();

    triRapid(S, 0, (Taille-1));

    temps = ((double)(clock() - sec)/CLOCKS_PER_SEC);
    printf("Le temps de calcul du tri de S est: %f secondes\n", temps);
    /*
    printf("S trié:");
    Affiche_Ensemble_ord(S,Taille,n);
    
    printf("S trié en décimale :");
    for(uint64_t i=0; i<Taille; i++)
    {
        printf("{%llu}", S[i]);
    }
    printf("\n");
    */
    


    
    /* x est défini comme le nombre total de bases de GJ de S de dimension d 
    Cette étape modifie également pt de nombreuses fois (sans toutefois garder un intérêt à la fin), 
    et surtout définit L comme la matrice de taille x fois d, où les lignes sont les bases distinctes de GJ */
    sec = clock();

    uint64_t x=GJBExtraction(S,d,Taille,pt,L,n,phi);

    temps = ((double)(clock() - sec)/CLOCKS_PER_SEC);
    printf("Le temps de calcul de GJBExtraction est: %f secondes\n", temps);

    /* Affichage de x et de toutes les bases de GJ calculées à l'étape précédente */
    printf("Le nombre de GJB de dimension %d dans S vaut %llu\n",d,x);
    for(uint64_t j=0;j<x;j++){
        printf("Le %llu-ème sous espace de dimension %d est : ",j+1,d);
        Affiche_Ensemble_ord(L[0][j],d,n);
    }
}
