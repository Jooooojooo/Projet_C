/*_________________________________________________________________________________________________
Ce fichier contient la fonction qui calcule les différents sous espace de dimension d d'un ensemble S
VERSION CONTENANT: Fonction GJBExtraction appliquée sur la boite S russe avec parallélisation.

LA PARALLÉLISATION a été effectuée avec la bibliothèque OpenMP.

SELON VOTRE MACHINE ET VERSION DE C: l'affichage des GJB dans le main peut être soit réalisée avec
"%lu", soit avec "%llu". Vérifiez et essayez les deux versions en cas de message d'erreur lié.
_______________________________________________________________________________________________________________________
POUR COMPILER SUR MAC (via homebrew):

TERMINAL GLOBAL: brew install libomp

PUIS:
gcc -O3 -Xpreprocessor -fopenmp -lomp Version_SBox.c -I/opt/homebrew/Cellar/libomp/17.0.6/include -o Version_SBox -L/opt/homebrew/Cellar/libomp/17.0.6/lib

PUIS:
./Version_SBox (d)

<!> attention: changement potentiellement nécessaire des -I/ et -L/ (et du numéro de version) selon lieu du téléchargement de libomp <!>

_______________________________________________________________________________________________________________________
POUR COMPILER SUR WINDOWS (via ubuntu):

gcc -O3 -fopenmp -o Version_SBox Version_SBox.c

PUIS:
./Version_SBox (d)
___________________________________________________________________________________________________*/

/* Inclusion des bibliothèques standards de C */
#include <stdio.h>
#include <stdlib.h>         
#include <time.h>
#include <stdint.h>
#include <omp.h>


/* 
Fonction d'affichage d'élément de F_2^n - dans le bon ordre par rapport à la lecture des bits (poids faible à droite).
_____________________________________________________________________________
UTILISE: /
ENTRÉE: un entier "n" <= 64 et un élément "v" de F_2^n.
SORTIE: aucune (fonction void)
FAIT: imprime v en représentation binaire. 
EXEMPLE: n=3, v=7, la fonction imprime 111.
_____________________________________________________________________________*/
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

/*
Fonction d'affichage d'un ensemble de vecteurs.
_____________________________________________________________________________
UTILISE: Affiche_Vecteur_ord()
ENTRÉE: un entier "n" <= 64, un ensemble de vecteurs de F_2^n, la taille de cet ensemble.
SORTIE: aucune (fonction void)
FAIT: imprime l'ensemble sous la forme {v_0}...{v_k}\n avec v_i imprimé comme dans Afficher_Vecteur_ord.
Exemple: si n=3, S = {6,7} avec 2 éléments, alors la fonction retourne {110}{111}\n
_____________________________________________________________________________*/
void Affiche_Ensemble_ord(uint64_t *pt, uint64_t length, int n){ 
    for(uint64_t i=0; i<length; i++){
        printf("{");
        Affiche_Vecteur_ord(pt[i],n);
        printf("}");
    }
    if(length ==0){
        printf("{}\n");
    }
    else{
        printf("\n");
    }
}

/*
Fonction simple retournant le Most Significant Bit d'un élément a de F_2^n.
_____________________________________________________________________________
UTILISE: /
ENTRÉE: un entier "n" <= 64 et un élément "a" de F_2^n
SORTIE: un entier "msb" <= 64, le Most Significant Bit de a
_____________________________________________________________________________*/
int MSB(uint64_t a, int n){
    int msb = 0;
    for(int k = 0; k<n; k++)
        if (((a>>k)&1) == 1)  msb = k;
    return msb;
}

/*
Fonction phi_d(S) sous forme booléenne, améliorée pour S trié.
_____________________________________________________________________________
UTILISE: MSB()
ENTRÉE: 
    - un entier "n" <= 64
    - un entier "d" <= n
    - un ensemble "S" de vecteurs de F_2^n
    - un entier "Taille", le nombre d'éléments dans S
    - une matrice '3D' d'entiers arbitraires "pt" de taille (#coeurs x d x Taille)
    - un entier "thread_id", le coeur en cours d'utilisation
SORTIE: aucune (fonction void)
FAIT:
    - calcule la liste intermédiaire msb_liste, de taille n, dont le k-ième élément est le nombre d'éléments de S 
    ayant k pour Most Significant Bit
    - utilise cette liste pour vérifier, pour chaque élément S[j] de S, la condition du Lemme 2.3.1.
    - si la condition est vérifiée, change pt[thread_id][d-1][j] en 1, sinon, le change en 0.
NOTES: 
    - ce choix de présentation booléenne plutôt que de produire l'ensemble phi_d(S) permettera paradoxalement de faciliter
    la vérification de l'appartenance à cet ensemble dans GJBExtraction: il suffira en effet d'appeler phi_d() sur une 
    matrice préalablement définie, puis de faire for(i){ if(pt[thread_id][d-1][i] == 1), do ...}
    - Bien qu'un appel unique de phi_d() ne change qu'une seule ligne de la matrice pt, nous le gardons sous forme 
    matricielle afin de permettre à la fonction de se "souvenir" des valeurs précédentes lors de la récursivité, comme
    expliqué pour chi_a(S) à la fin de la Section 3.2.
_____________________________________________________________________________*/
void phi_d_Tri(uint64_t* S, int d, uint64_t*** pt, int n, uint64_t Taille, int thread_id){
    uint64_t msb_liste[n]; 
    short prev_msb = -1;
    short prev_result = -1;
    for(int i = 0; i < n; i++){
        msb_liste[i] = 0;
    }
    for(uint64_t k = 0; k < Taille; k++) msb_liste[MSB(S[k], n)]++;
    for(uint64_t j = 0; j < Taille; j++){
        short msb = MSB(S[j], n);
        short result;  
        if(msb == prev_msb){
            result = prev_result;
        }
        else{
            int i = 1;
            for(int k = msb+1; k < n; k++){
                if (msb_liste[k]>(1<<i&1ULL)) i++;        
            }
            result = (i > d - 1) ? 1 : 0;
            prev_msb = msb;
            prev_result = result; 
        }
        pt[thread_id][d-1][j] = result;
    }    
}

/*
Fonction de permutation de deux éléments dans un vecteur, nécessaire pour le tri rapide avec pivot
_____________________________________________________________________________
UTILISE: /
ENTRÉE: deux entiers, *a et *b
SORTIE: aucune (fonction void)
FAIT: inverse *a et *b
_____________________________________________________________________________*/
void permuter(uint64_t *a, uint64_t *b){
    uint64_t temp = *a;
    *a = *b;
    *b = temp;
}

/*
Fonction de tri rapide d'un ensemble, de type QuickSort avec 1 pivot
_____________________________________________________________________________
UTILISE: permuter()
ENTRÉE: un ensemble d'entiers "tab", deux entiers "first" et "last", les bornes de la partie de S à trier.
SORTIE: aucune (fonction void).
FAIT: trie rapidement tab à l'aide d'un pivot via la méthode standard de QuickSort.
NOTE: en pratique l'on désire toujours trier entièrement l'ensemble, donc first=0 et last=(taille de l'ensemble)-1
_____________________________________________________________________________*/
void triRapid(uint64_t tab[], uint64_t first, uint64_t last){
    uint64_t pivot, i, j;
    if(first < last){
        pivot = last;
        i = first;
        j = last;
        while (i < j){
            while(tab[i] <= tab[pivot] && i < last)
                i++;
            while(tab[j] > tab[pivot])
                j--;
            if(i < j){
                permuter(&tab[i], &tab[j]);
            }
        }
        permuter(&tab[pivot], &tab[j]);
        triRapid(tab, first, j - 1);
        triRapid(tab, j + 1, last);
    }
}

/* Fonction Chi_a(S) définie en Section 2.2, améliorée pour S trié.
_____________________________________________________________________________
UTILISE: /
ENTRÉE: 
    - un entier "a", élément de F_2^n
    - un sous-ensemble "S" de vecteurs de F_2^n
    - un entier "Taille", le nombre d'éléments de S
    - un entier "d"
    - une 'matrice 3D' d'entiers arbitraires "pt" de taille (#coeurs x d x Taille)
    - un entier "thread_id", le coeur en cours d'utilisation

SORTIE: Un entier "k", le nombre d'éléments contenus dans l'ensemble chi_a(S)
FAIT:
    - Vérifie pour chaque élément de S la condition d'appartenance à Chi_a(S) expliquée en Section 2.2
    - À chaque nouvel élément vérifiant ces conditions, définit pt[thread_id][d-1][k] comme cet élément puis incrémente k
NOTES:
    - Bien qu'un appel unique de chi_a_S() ne change qu'une seule ligne de la matrice pt, nous le gardons sous forme 
    matricielle afin de permettre à la fonction de se "souvenir" des valeurs précédentes lors de la récursivité, comme
    expliqué à la fin de la Section 3.2.
    - Si l'on désire calculer et utliser la borne "optimale" de la fin de la Section 2.2 afin de comparer les temps 
    de calculs, il faut rajouter un argument "int n" à la fonction. Le calcul et l'utilisation de cette borne sont 
    présents en commentaires dans la fonction - il suffit de la décommenter et de commenter la boucle au dessus pour 
    passer à l'autre version.
_____________________________________________________________________________*/
uint64_t chi_a_S_Tri(uint64_t a, uint64_t* S, uint64_t*** pt , uint64_t Taille, int d, int thread_id){
    uint64_t k = 0;
    uint64_t m = 0;
    while ((S[m]) <= a)     /* Méthode de recherche linéaire du 1er indice m tel que S[m] > a*/
    {
        m++;
    }
    /*  MÉTHODE DE RECHERCHE DICHOTOMIQUE & BORNE OPTIMALE - BIEN PLUS LENTE EN PRATIQUE
    if (MSB(a,n) == n-1){
        return 0;
    }
    int msb_power = (1ULL<<(MSB(a, n) + 1));
    uint64_t low = 0;
    uint64_t m = Taille;

    while(low<m){           
        uint64_t mid = low +(m-low)/2.0;
        if(S[mid]<msb_power){
            low = mid + 1;
       }
        else m = mid;
    }
    */
    for(uint64_t i = m; i<Taille; i++){
        if (((S[i])^a) > S[i]){
            for(uint64_t j = i+1 ; j<Taille; j++){
                if (((S[i])^a) == S[j]){
                    pt[thread_id][d-1][k] = S[i];
                    k++;
                }
            }
        }
    }
    return(k);
}

/* Fonction de Somme de valeurs de compteurs, utilisée pour uniformiser les compteurs selon les threads
_____________________________________________________________________________
UTILISE: /
ENTRÉE: une matrice d'entiers L et un entier d
SORTIE: la somme des éléments de la d-ème colonne de L
_____________________________________________________________________________*/
uint64_t Somme(uint64_t **L,int d){
    uint64_t S=0;
    for(int i=0;i<omp_get_max_threads();i++ ){
        S = S + L[i][d];
    }
    return(S);
}

/*
Fonction GJBExtraction définie en Section 3.2, fonctionnant pour un ensemble S général (pas nécessairement trié)
_____________________________________________________________________________
UTILISE: phi_d(), chi_a_S()
ENTRÉE: 
    - un sous-ensemble "S" de vecteurs de F_2^n
    - un entier "n" <= 64
    - un entier "d", avec 0 < d <= n
    - un entier "Taille", le nombre d'éléments de S
    - une 'matrice 3D' d'entiers arbitraires "pt" de taille (#coeurs x d x Taille)
    - une 'matrice 3D' "L", sans éléments prédéfinis, mais avec L[thread_id] de taille (1024 x d), pour chaque thread_id (ie,
    nombre de coeurs total). 
    - une matrice d'entiers "compteurs", tous initialisés à 0, de taille (#coeurs x d)
    - une 'matrice 3D' d'entiers arbitraires "phi" de taille (#coeurs x d x Taille)
    - un entier "multithreads" égal à 1 au premier appel
    - un entier "d_const" égal à d
SORTIE:
    - Un entier "Taille_L_new", le nombre de bases de Gauss-Jordan distinctes de S de dimension d
FAIT:
    - Au premier appel, applique la parallélisation et change le multithreads en 0 (ce qui veut dire que la parallélisation
    n'est appliquée qu'au niveau le plus élevé). L'ensemble S se voit donc découpé et réparti dans chaque thread_id
    - pour chaque thread_id, calcule phi_d(S) et regarde la valeur de d.
    - si d=1, ajoute chaque élément de S dans la colonne 0 de L[thread_id], et incrémente compteurs[thread_id][0] du nombre 
    d'éléments ajoutés. Dans le cas où ce nombre atteint un multiple de 1024, réalloque également 1024 "lignes" supplémentaires 
    de mémoire à L[thread_id]. Retourne ensuite le nombre d'éléments ajoutés.
    - si d≠1, pour chaque élément a de S, calcule chi_a(S). Si sa taille est au moins 2^{d-1}-1, appelle GJBExtraction
    à nouveau, pour cet ensemble et pour d-1, et rajoute cet élément dans la colonne d-1 de L[thread_id] un certain nombre de fois
    (dépendant du nombre d'éléments ajoutés pour tous ses chi respectifs) et incrémente compteurs[d-1] de ce nombre.
    Retourne ensuite le nombre de "nouveaux" éléments ajoutés
NOTES:
    - Cette allocation "dynamique" de mémoire pour L est nécessaire afin de ne pas gâcher trop de mémoire, vu que l'on a
    aucune information sur le nombre de GJB contenues dans un ensemble avant leurs calculs.
    - d_const ne sert uniquement qu'à la réallocation de mémoire de L[thread_id] (pour que la fonction sache combien de 
    "colonnes" L[thread_id] possède)
    - La matrice compteurs[] sert à savoir où mettre le prochain élément dans L[thread_id], et doit donc pour cela être
    "en dehors" de la fonction afin de ne pas être affecté par la récursivité.
_____________________________________________________________________________*/
uint64_t GJBExtraction_Tri(uint64_t* S, int d, uint64_t Taille, uint64_t***pt, uint64_t*** L, int n, uint64_t** compteurs, uint64_t*** phi, int multithreads, int d_const){
    uint64_t Taille_L_new=0;
    uint64_t Taille_L_old=0;
    uint64_t Taille_S_a;
    int thread_id;
    if((d-1)!=0){
        if(multithreads==1){     
            phi_d_Tri(S,d,phi,n,Taille,0);    
            for(int i=1;i<omp_get_max_threads();i++){
                for(uint64_t k=0; k<Taille;k++){
                    phi[i][d-1][k] = phi[0][d-1][k];} 
            }
            #pragma omp parallel for private(Taille_S_a, thread_id, Taille_L_new, Taille_L_old) schedule(dynamic)
            for (uint64_t i=0; i<Taille; i++){   
                thread_id=omp_get_thread_num();
                if (phi[thread_id][d-1][i]==1){    
                    Taille_S_a = chi_a_S_Tri(S[i],S,pt,Taille,d,thread_id);
                    if (Taille_S_a >= ((1ULL<<(d-1))-1)){
                        Taille_L_old=GJBExtraction_Tri(pt[thread_id][d-1], d-1,Taille_S_a,pt,L,n,compteurs,phi,0,d_const);
                        for(uint64_t j=(compteurs[thread_id][d-1]); j<(Taille_L_old+compteurs[thread_id][d-1]); j++){ 
                            L[thread_id][j][d-1]=S[i];
                        } 
                        compteurs[thread_id][d-1]=Taille_L_old+compteurs[thread_id][d-1];
                    }
                }
            }          
            uint64_t S=Somme(compteurs,d-1);
            return(S);
        }
        else{
            thread_id=omp_get_thread_num();
            phi_d_Tri(S,d,phi,n,Taille,thread_id);
            for (uint64_t i=0; i<Taille; i++){          
                if (phi[thread_id][d-1][i]==1){    
                    Taille_S_a = chi_a_S_Tri(S[i],S,pt,Taille,d,thread_id);
                    if (Taille_S_a >= ((1ULL<<(d-1))-1)){    
                        Taille_L_old=GJBExtraction_Tri(pt[thread_id][d-1], d-1,Taille_S_a,pt,L,n,compteurs,phi,0,d_const)+Taille_L_new;
                        for(uint64_t j=(compteurs[thread_id][d-1]+Taille_L_new); j<(Taille_L_old+compteurs[thread_id][d-1]); j++){  
                            L[thread_id][j][d-1]=S[i];
                        }              
                        Taille_L_new = Taille_L_old;    
                    }
                }
            } 
            compteurs[thread_id][d-1]=Taille_L_old+compteurs[thread_id][d-1]; 
            return(Taille_L_new);
        }
    }
    if(d ==1){   
        thread_id=omp_get_thread_num();  
        for(uint64_t l=compteurs[thread_id][0];l<(Taille+compteurs[thread_id][0]);l++){
            if(((l+1)%1024 == 0) & ((l+1)<=(Taille+compteurs[thread_id][0]))){ 
                L[thread_id]=(uint64_t **)realloc(L[thread_id],(l+1+1024)*sizeof(uint64_t*));
                for(int q=0;q<1024;q++){
                    L[thread_id][l+1+q]=malloc(d_const*sizeof(uint64_t));
                }
            }
            L[thread_id][l][0]=S[l-compteurs[thread_id][0]];
        }      
        compteurs[thread_id][0] = Taille+compteurs[thread_id][0]; 
        return(Taille);
    } 
}

/* Fonction de calcul produit scalaire d'éléments de F_2^8 */
int scalar_prod(int x, int y){
    int temp_1=x;
    int temp_2=y;
    int prod = 0;
    for(int i=0;i<8;i++){
        prod =prod ^ ((temp_1&1)*(temp_2&1));
        temp_1= temp_1>>1;
        temp_2= temp_2>>1;
    }
    return(prod);
}

/* Fonction générant dans une matrice L de taille 16 x 16 la boite S russe */
void S_box_Russe(int **L){
    int valeurs[] = {
        0xFC, 0xEE, 0xDD, 0x11, 0xCF, 0x6E, 0x31, 0x16, 0xFB, 0xC4, 0xFA, 0xDA, 0x23, 0xC5, 0x04, 0x4D,
        0xE9, 0x77, 0xF0, 0xDB, 0x93, 0x2E, 0x99, 0xBA, 0x17, 0x36, 0xF1, 0xBB, 0x14, 0xCD, 0x5F, 0xC1,
        0xF9, 0x18, 0x65, 0x5A, 0xE2, 0x5C, 0xEF, 0x21, 0x81, 0x1C, 0x3C, 0x42, 0x8B, 0x01, 0x8E, 0x4F,
        0x05, 0x84, 0x02, 0xAE, 0xE3, 0x6A, 0x8F, 0xA0, 0x06, 0x0B, 0xED, 0x98, 0x7F, 0xD4, 0xD3, 0x1F,
        0xEB, 0x34, 0x2C, 0x51, 0xEA, 0xC8, 0x48, 0xAB, 0xF2, 0x2A, 0x68, 0xA2, 0xFD, 0x3A, 0xCE, 0xCC,
        0xB5, 0x70, 0x0E, 0x56, 0x08, 0x0C, 0x76, 0x12, 0xBF, 0x72, 0x13, 0x47, 0x9C, 0xB7, 0x5D, 0x87,
        0x15, 0xA1, 0x96, 0x29, 0x10, 0x7B, 0x9A, 0xC7, 0xF3, 0x91, 0x78, 0x6F, 0x9D, 0x9E, 0xB2, 0xB1,
        0x32, 0x75, 0x19, 0x3D, 0xFF, 0x35, 0x8A, 0x7E, 0x6D, 0x54, 0xC6, 0x80, 0xC3, 0xBD, 0x0D, 0x57,
        0xDF, 0xF5, 0x24, 0xA9, 0x3E, 0xA8, 0x43, 0xC9, 0xD7, 0x79, 0xD6, 0xF6, 0x7C, 0x22, 0xB9, 0x03,
        0xE0, 0x0F, 0xEC, 0xDE, 0x7A, 0x94, 0xB0, 0xBC, 0xDC, 0xE8, 0x28, 0x50, 0x4E, 0x33, 0x0A, 0x4A,
        0xA7, 0x97, 0x60, 0x73, 0x1E, 0x00, 0x62, 0x44, 0x1A, 0xB8, 0x38, 0x82, 0x64, 0x9F, 0x26, 0x41,
        0xAD, 0x45, 0x46, 0x92, 0x27, 0x5E, 0x55, 0x2F, 0x8C, 0xA3, 0xA5, 0x7D, 0x69, 0xD5, 0x95, 0x3B,
        0x07, 0x58, 0xB3, 0x40, 0x86, 0xAC, 0x1D, 0xF7, 0x30, 0x37, 0x6B, 0xE4, 0x88, 0xD9, 0xE7, 0x89,
        0xE1, 0x1B, 0x83, 0x49, 0x4C, 0x3F, 0xF8, 0xFE, 0x8D, 0x53, 0xAA, 0x90, 0xCA, 0xD8, 0x85, 0x61,
        0x20, 0x71, 0x67, 0xA4, 0x2D, 0x2B, 0x09, 0x5B, 0xCB, 0x9B, 0x25, 0xD0, 0xBE, 0xE5, 0x6C, 0x52,
        0x59, 0xA6, 0x74, 0xD2, 0xE6, 0xF4, 0xB4, 0xC0, 0xD1, 0x66, 0xAF, 0xC2, 0x39, 0x4B, 0x63, 0xB6
    };
    int k = 0;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            L[i][j] = valeurs[k++];
        }
    }
}

/* Fonction calculant la valeur de W_F(a,b) où F est implémentée dans L */
int W_F(int a,int b,int**L){
    int Resultat=0;
    int x;
    int temp_1;
    int temp_2;
    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            x=i+16*j;
            temp_1= scalar_prod(a,x);
            temp_2= scalar_prod(b,L[j][i]);
            temp_1= temp_1 ^ temp_2 ;
            temp_2=temp_1 & 1;
            if(temp_2==0){
                Resultat = Resultat +1;}
            else{
                Resultat = Resultat -1; 
            }
        }
    }
    return(Resultat);
}


/* Fonction de création de Z_F, mis dans le vecteur S, et qui retourne sa taille. Vu que Z_F est dans F_2^8xF_2^8,
pour le transformer en sous ensemble de F_2^16, on fait le choix de coordonné suivant : (a||b)*/
int Creation_S(int** S_0,uint64_t* S){
    int compteurs=0;
    int temp;
    for(int a=0;a<256;a++){
        for(int b=0;b<256;b++){
            temp=W_F(a,b,S_0);
            if(temp==0){
                S[compteurs]=(a<<8)^b;
                compteurs+=1;
            }
        }
    }
    return(compteurs);
}

int main(int argc, char const *argv[]){
    /* UN SEUL ARGUMENT ICI: d = argv[1] */
    /* Déclaration table de la S box russe et remplissage */
    int** S_0=malloc(16*sizeof(int*));
    for(int i=0;i<16;i++){
        S_0[i]=malloc(16*sizeof(int));
    }
    S_box_Russe(S_0);

    /* Déclaration et construction de S */
    uint64_t* S=malloc((1<<16)*sizeof(uint64_t));

    /* Déclaration des variables n, Taille, d comme entrées de argv */
    int n=16;
    uint64_t Taille = Creation_S(S_0,S);
    int d=atoi(argv[1]);

    /* Affichage des zéros de Walsh
    printf("L'ensemble des zéros de Walsh de la LAT de la S Box est: \n");
    for(int i=0; i < Taille; i++){
        printf("(%llu,%llu) ",(S[i] >> 8), ((S[i])%256));
    }
    */
    printf("La Taille de S est: %llu\n", Taille);
    /* Déclaration de compteurs, L, pt, phi*/
    uint64_t** compteurs=malloc(omp_get_max_threads()*sizeof(uint64_t*));
    for(int j=0;j<omp_get_max_threads();j++){
        compteurs[j]=malloc(d*sizeof(uint64_t*));
        for(int i=0; i<d; i++){
            compteurs[j][i] = 0;
        }
    }

 
    uint64_t*** L=malloc(omp_get_max_threads()*sizeof(uint64_t**));
    for(int i=0;i<omp_get_max_threads();i++){
        L[i]=malloc(1024*sizeof(uint64_t*));    
        for(int j=0;j<1024;j++){
            L[i][j]=malloc(d*sizeof(uint64_t));
        }
    }

  
    uint64_t*** pt=malloc(omp_get_max_threads()*sizeof(uint64_t**));
    for(int j=0;j<omp_get_max_threads();j++){
        pt[j]=malloc(d*sizeof(uint64_t*));
        for(int i=0; i<d; i++){
            pt[j][i]=malloc(Taille*sizeof(uint64_t));
        }
    }


    uint64_t*** phi=malloc(omp_get_max_threads()*sizeof(uint64_t**));
    for(int j=0;j<omp_get_max_threads();j++){
        phi[j]=malloc(d*sizeof(uint64_t*));
        for(int i=0; i<d; i++){
            phi[j][i]=malloc(Taille*sizeof(uint64_t));
            for(uint64_t k=0;k<Taille;k++){
                phi[j][i][k]=0;
                }
            }
    }


    /* Tri de S - temps négligeable */
    triRapid(S, 0, Taille-1);
    

    /* Déclaration du temps */
    double temps = omp_get_wtime();

    uint64_t x=GJBExtraction_Tri(S,d,Taille,pt,L,n,compteurs,phi,1,d);

    /* Obtention et affichage du temps de calcul de la fonction GJBExtraction (en secondes)*/
    temps -= omp_get_wtime();
    temps *= -1; 
    printf("Le temps de calcul de GJBExtraction_Tri est: %f secondes\n", temps);

    /* Affichage de x et de toutes les bases de GJ calculées à l'étape précédente */
    
    for(int i=0;i<omp_get_max_threads();i++){
        int somme_totale = 0;
        for(int k = 0; k<i; k++){
            somme_totale += compteurs[k][d-1];
        }
        for(uint64_t j=0;j<compteurs[i][d-1];j++){
            printf("Le %llu-ème sous espace de dimension %d associé au coeur %d est : ",j+1+somme_totale,d,i);
            Affiche_Ensemble_ord(L[i][j],d,n);
        }
    }

    printf("Le nombre de GJB de dimension %d dans S vaut %llu\n",d,x);


    free(S);

    for(int i=0;i<16;i++){
        free(S_0[i]);
    }
    free(S_0);
    
    
    for(int i=0;i<omp_get_max_threads();i++){
        for(uint64_t j=0;j< ((compteurs[i][0] / 1024) + 1)*1024;j++){
            free(L[i][j]);
        }
        free(L[i]);
    }
    free(L);
    
    for(int i=0;i<omp_get_max_threads();i++){
        free(compteurs[i]);
    }
    free(compteurs);

    for(int i=0;i<omp_get_max_threads();i++){
        for(int j=0;j<d;j++){
            free(pt[i][j]);
        }
        free(pt[i]);
    }
    free(pt);

    for(int i=0;i<omp_get_max_threads();i++){
        for(int j=0;j<d;j++){
            free(phi[i][j]);
        }
        free(phi[i]);
    }
    free(phi);
}