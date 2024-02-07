/*_________________________________________________________________________________________________
Ce fichier contient la fonction qui calcule les différents sous espace de dimension d d'un ensemble S
VERSION CONTENANT: Fonction GJBExtraction optimisée pour le cas où S est trié avec parallélisation.

LA PARALLÉLISATION a été effectuée avec la bibliothèque OpenMP.

DANS CE FICHIER, S est supposé contenir 0 à chaque fois par définition pour pouvoir avoir des sous-espaces
vectoriels. Afin de faciliter l'écriture des fonctions, nous avons fait le choix de ne jamais écrire 0
dans les éléments de S (il sera toujours "implicitement" supposé être là). Si vous décidez d'utiliser un 
S personnalisé, définissez le bien comme l'ensemble des éléments non nuls qu'il contient.

SELON VOTRE MACHINE ET VERSION DE C: l'affichage des GJB dans le main peut être soit réalisée avec
"%lu", soit avec "%llu". Vérifiez et essayez les deux versions en cas de message d'erreur lié.
_______________________________________________________________________________________________________________________
POUR COMPILER SUR MAC (via homebrew):

TERMINAL GLOBAL: brew install libomp

PUIS:
gcc -O3 -Xpreprocessor -fopenmp -lomp Version_Parallelisee.c -I/opt/homebrew/Cellar/libomp/17.0.6/include -o Version_threads -L/opt/homebrew/Cellar/libomp/17.0.6/lib

PUIS:
./Version_Parallelisee (n) (Taille) (d)

<!> attention: changement potentiellement nécessaire des -I/ et -L/ (et du numéro de version) selon lieu du téléchargement de libomp <!>

_______________________________________________________________________________________________________________________
POUR COMPILER SUR WINDOWS (via ubuntu):

gcc -O3 -fopenmp -o Version_Parallelisee Version_Parallelisee.c

PUIS:
./Version_Parallelisee (n) (Taille) (d)
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
Générateur "rrmxmx" pseudo-aléatoire de nombres binaires de 64 bits pour mesure de performances de temps de calculs
_____________________________________________________________________________
UTILISE: /
ENTRÉE: un entier "n" de 64 bits maximum, la "graine"
SORTIE: un entier de 64 bits pseudo-aléatoire
NOTE: Cette fonction est déterministe et fonctionne comme rand(): une même graine génère toujours le même nombre.
_____________________________________________________________________________*/
uint64_t rnd64(uint64_t n){
    const uint64_t z = 0x9FB21C651E98DF25;
    n ^= ((n << 49) | (n >> 15)) ^ ((n << 24) | (n >> 40));
    n *= z;
    n ^= n >> 35;
    n *= z;
    n ^= n >> 28;
    return n;
}


/*
Fonction simple retournant le Most Significant Bit d'un élément a de F_2^n.
_____________________________________________________________________________
UTILISE: /
ENTRÉE: un entier "n" <= 64 et un élément "a" de F_2^n
SORTIE: un entier "msb" <= 64, le Most Significant Bit de a
_____________________________________________________________________________*/
int MSB(uint64_t a, int n){
    short msb = 0;
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
uint64_t Somme(uint64_t **L, int d){
    uint64_t S = 0;
    for(int i=0; i < omp_get_max_threads(); i++){
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
        for(uint64_t l = compteurs[thread_id][0]; l < (Taille + compteurs[thread_id][0]); l++){
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

int main(int argc, char const *argv[]){
    /*_____________________________________________________________________________
    Pas de contraintes sur argc, et on s'attend à ce que: 
    - argv[1] = n; 
    - argv[2] = Taille de S; 
    - argv[3] = d
    _____________________________________________________________________________*/
    int n=atoi(argv[1]);
    uint64_t Taille = atoi(argv[2]);
    int d=atoi(argv[3]);
    
    /* Création de S */
    uint64_t* S=malloc((Taille)*sizeof(uint64_t));

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
            for(uint64_t k = 0; k < Taille; k++){
                phi[j][i][k]=0;
            }
        }
    }
    
    /*_____________________________________________________________________________
    Ci-dessous, définition des éléments de S. Par défaut ici, S est défini comme l'ensemble {1, ..., Taille}
    Disponible également en dessous en commentaires, une définition de S utilisant la fonction aléatoire de 64 bits. 
    ATTENTION, utiliser cette fonction seulement avec n=64 pour éviter d'obtenir deux éléments égaux dans S
    Si vous décidez de définir un S personnalisé, définissez le ici et vérifiez bien que son nombre d'élément
    est égal à l'argument Taille = argv[2] du main.
    _____________________________________________________________________________*/
    for (uint64_t i=0;i<Taille;i++){
        S[i] = i+1;
    }
    /*
    uint64_t state = 1;
    for (uint64_t i=0;i<Taille;i++)
    {
        const uint64_t n = rnd64(state++);
        S[i] = n;
    }
    */


    /* Tri de S - temps négligeable */
    triRapid(S, 0, Taille-1);

    /* Déclaration du temps - OMP est nécessaire ici car time.h ne gère pas bien la parallélisation */
    double temps = omp_get_wtime();

    uint64_t x=GJBExtraction_Tri(S,d,Taille,pt,L,n,compteurs,phi,1,d);

    /* Obtention et affichage du temps de calcul de la fonction GJBExtraction (en secondes)*/
    temps -= omp_get_wtime();
    temps *= -1; 

    /* Affichage des bases trouvées - ATTENTION: peut tuer votre terminal si vous avez utilisé la
    version S[i] = {1, ..., Taille} avec n et d trop grands */
    for(int i=0;i<omp_get_max_threads();i++){
        int somme_totale = 0;
        for(uint64_t k = 0; k<i; k++){
            somme_totale += compteurs[k][d-1];
        }
        for(uint64_t j=0;j<compteurs[i][d-1];j++){
            printf("Le %llu-ème sous espace de dimension %d associé au coeur %d est : ",j+1+somme_totale,d,i);
            Affiche_Ensemble_ord(L[i][j],d,n);}
    }

    printf("Le nombre de GJB de dimension %d dans S vaut %llu\n",d,x);

    /*
    ATTENTION - SI UTILISATION MULTIPLE DE LA FONCTION GJBEXTRACTION DANS LE MAIN - 
    BIEN REMETTRE TOUTES LES VALEURS DE COMPTEURS À ZÉRO ENTRE CHAQUE APPEL !!
    */

    /* Libération de mémoire des vecteurs et matrices malloquées au début du main */
    free(S);
    
    for(int i=0;i<omp_get_max_threads();i++){
        for(uint64_t j = 0; j < ((compteurs[i][0] / 1024) + 1)*1024; j++){
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