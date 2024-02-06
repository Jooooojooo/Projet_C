/*_________________________________________________________________________________________________
Ce fichier contient la fonction qui calcule les différents sous espace de dimension d d'un ensemble S
VERSION CONTENANT: Fonction GJBExtraction standard, Fonction GJBExtraction_Tri optimisée 
pour le cas où S est trié, pas de parallélisation.

DANS CE FICHIER, S est supposé contenir 0 à chaque fois par définition pour pouvoir avoir des sous-espaces
vectoriels. Afin de faciliter l'écriture des fonctions, nous avons fait le choix de ne jamais écrire 0
dans les éléments de S (il sera toujours "implicitement" supposé être là). Si vous décidez d'utiliser un 
S personnalisé, définissez le bien comme l'ensemble des éléments non nuls qu'il contient.

SELON VOTRE MACHINE ET VERSION DE C: l'affichage des GJB dans le main peut être soit réalisée avec
"%lu", soit avec "%llu". Vérifiez et essayez les deux versions en cas de message d'erreur lié. 
___________________________________________________________________________________________________*/


/* Inclusion des bibliothèques standards de C */
#include <stdio.h>
#include <stdlib.h>         
#include <time.h>
#include <stdint.h>

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
    i = 1ULL << (n - 1);
    while(i > 0){
        if(v & i)
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
ENTRÉE: un entier "n" <= 64, un ensemble de vecteurs de F_2^n et la taille de cet ensemble.
SORTIE: aucune (fonction void)
FAIT: imprime l'ensemble sous la forme {v_0}...{v_k}\n avec v_i imprimé comme dans Afficher_Vecteur_ord.
Exemple: si n = 3, S = {6,7} avec 2 éléments, alors la fonction retourne {110}{111}\n
_____________________________________________________________________________*/
void Affiche_Ensemble_ord(uint64_t *pt, uint64_t size, int n){ 
    for(uint64_t i = 0; i < size; i++){
        printf("{");
        Affiche_Vecteur_ord(pt[i], n);
        printf("}");
    }
    if(size == 0){
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
    for(int k = 0; k < n; k++)
        if (((a >> k) & 1ULL) == 1)  msb = k;
    return msb;
}

/*
Fonction phi_d(S) sous forme booléenne, fonctionnant pour un ensemble S général (pas nécessairement trié)
_____________________________________________________________________________
UTILISE: MSB()
ENTRÉE: 
    - un entier "n" <= 64
    - un entier "d" <= n
    - un ensemble "S" de vecteurs de F_2^n
    - un entier "Taille", le nombre d'éléments dans S
    - une matrice d'entiers arbitraires "pt" de taille (d x Taille)
SORTIE: aucune (fonction void)
FAIT:
    - calcule la liste intermédiaire msb_liste, de taille n, dont le k-ième élément est le nombre d'éléments de S 
    ayant k pour Most Significant Bit
    - utilise cette liste pour vérifier, pour chaque élément S[j] de S, la condition du Lemme 2.3.1.
    - si la condition est vérifiée, change pt[d-1][j] en 1, sinon, le change en 0.
NOTES: 
    - ce choix de présentation booléenne plutôt que de produire l'ensemble phi_d(S) permettera paradoxalement de faciliter
    la vérification de l'appartenance à cet ensemble dans GJBExtraction: il suffira en effet d'appeler phi_d() sur une 
    matrice préalablement définie, puis de faire for(i){ if(pt[d-1][i] == 1), do ...}
    - Bien qu'un appel unique de phi_d() ne change qu'une seule ligne de la matrice pt, nous le gardons sous forme 
    matricielle afin de permettre à la fonction de se "souvenir" des valeurs précédentes lors de la récursivité, comme
    expliqué pour chi_a(S) à la fin de la Section 3.2.
_____________________________________________________________________________*/
void phi_d(uint64_t* S, int d, uint64_t** pt, int n, uint64_t Taille){
    uint64_t msb_liste[n]; 
    for(int i = 0; i < n; i++){
        msb_liste[i] = 0;
    }
    for(uint64_t k = 0; k < Taille; k++) msb_liste[MSB(S[k], n)]++;
    for(uint64_t j = 0; j < Taille; j++){
        short msb = MSB(S[j], n);
        int i = 1;
        for (int k = msb + 1; k < n; k++){
            if (msb_liste[k] > (1 << i & 1ULL)) i++;        
        }
        if (i > d-1) pt[d-1][j] = 1;
        else pt[d-1][j] = 0;
    }    
}

/*
Fonction phi_d(S) sous forme booléenne, améliorée pour S trié (même structure et remarques que phi_d)
_____________________________________________________________________________*/
void phi_d_Tri(uint64_t* S, int d, uint64_t** pt, int n, uint64_t Taille){
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
            for (int k = msb + 1; k < n; k++){
                if (msb_liste[k] > ( 1<< i & 1ULL)) i++;        
            }
            result = (i > d-1) ? 1 : 0;
            prev_msb = msb;
            prev_result = result; 
        }
        pt[d-1][j] = result;
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

/*
Fonction Chi_a(S) définie en Section 2.2, fonctionnant pour un ensemble S général (pas nécessairement trié)
_____________________________________________________________________________
UTILISE: /
ENTRÉE: 
    - un entier "a", élément de F_2^n
    - un sous-ensemble "S" de vecteurs de F_2^n
    - un entier "Taille", le nombre d'éléments de S
    - un entier "d"
    - une matrice d'entiers arbitraires "pt" de taille (d x Taille)

SORTIE: Un entier "k", le nombre d'éléments contenus dans l'ensemble chi_a(S)
FAIT:
    - Vérifie pour chaque élément de S la condition d'appartenance à Chi_a(S) expliquée en Section 2.2
    - À chaque nouvel élément vérifiant ces conditions, définit pt[d-1][k] comme cet élément puis incrémente k
NOTES:
    - Bien qu'un appel unique de chi_a_S() ne change qu'une seule ligne de la matrice pt, nous le gardons sous forme 
    matricielle afin de permettre à la fonction de se "souvenir" des valeurs précédentes lors de la récursivité, comme
    expliqué à la fin de la Section 3.2.
_____________________________________________________________________________*/
uint64_t chi_a_S(uint64_t a, uint64_t* S, uint64_t** pt, uint64_t Taille, int d){
    uint64_t k = 0;
    for(uint64_t i = 0; i < Taille; i++){
        if (((S[i]) > a) & (((S[i]) ^ a) > S[i])){     /* Condition a < x < x xor a */
            for(uint64_t j = 0; j < Taille; j++){
                if (((S[i]) ^ a) == S[j]){            /* Condition (S[i] xor a) appartient à S */
                    pt[d-1][k] = S[i];
                    k++;
                }
            }
        }
    }
    return(k);
}

/*
Fonction Chi_a(S) définie en Section 2.2, améliorée pour S trié (même structure et remarques que chi_a_S())
_____________________________________________________________________________
NOTES: L'argument supplémentaire "n" n'est pas utilisé tel quel, mais est nécessaire pour faire appel à la fonction
MSB(), qui est elle même nécessaire si l'on désire calculer et utliser la borne "optimale" de la fin de la Section 2.2
afin de comparer les temps de calculs. Le calcul et l'utilisation de cette borne sont présents en commentaires dans
la fonction - il suffit de la décommenter et de commenter la boucle au dessus pour passer à l'autre version.
_____________________________________________________________________________*/
uint64_t chi_a_S_Tri(uint64_t a, uint64_t* S, uint64_t** pt , uint64_t Taille, int d, int n){
    uint64_t k = 0;
    uint64_t m = 0;
    while ((S[m]) <= a){     /* Méthode de recherche linéaire du 1er indice m tel que S[m] > a*/
        m++;
    }
    /*  MÉTHODE DE RECHERCHE DICHOTOMIQUE & BORNE OPTIMALE - BIEN PLUS LENTE EN PRATIQUE 
    if (MSB(a,n) == n - 1){
        return 0;
    }
    int msb_power = (1ULL << (MSB(a, n) + 1));
    uint64_t low = 0;
    uint64_t m = Taille;

    while(low < m){           
        uint64_t mid = low + (m - low)/2.0;
        if(S[mid] < a){
            low = mid + 1;
       }
        else m = mid;
    }
    */
    for(uint64_t i = m; i < Taille; i++){
        if (((S[i]) ^ a) > S[i]){                     /* Condition x < x xor a */
            for(uint64_t j = i + 1 ; j < Taille; j++){
                if (((S[i]) ^ a) == S[j]){            /* Condition (S[i] xor a) appartient à S */
                    pt[d-1][k] = S[i];
                    k++;
                }
            }
        }
    }
    return(k);
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
    - une matrice d'entiers arbitraires "pt" de taille (d x Taille)
    - une 'matrice 3D' "L", sans éléments prédéfinis, mais avec L[0] de taille (1024 x d). Même si l'on utilise pas 
    la 1ère des 3 dimensions de L (on n'utilise que L[0][..][..]), il est défini ainsi pour permettre de garder son 
    adresse en mémoire dans la récursivité lors d'une potentielle réallocation de mémoire. 
    - un vecteur de d entiers "compteurs", tous initialisés à 0
    - une matrice d'entiers arbitraires "phi" de taille (d x Taille)
    - un entier "d_const" égal à d
SORTIE:
    - Un entier "Taille_L_new", le nombre de bases de Gauss-Jordan distinctes de S de dimension d
FAIT:
    - calcule phi_d(S) et regarde la valeur de d.
    - si d=1, ajoute chaque élément de S dans la colonne 0 de L[0], et incrémente compteurs[0] du nombre d'éléments ajoutés. 
    Dans le cas où ce nombre atteint un multiple de 1024, réalloque également 1024 "lignes" supplémentaires de mémoire à L[0].
    Retourne ensuite le nombre d'éléments ajoutés.
    - si d≠1, pour chaque élément a de S, calcule chi_a(S). Si sa taille est au moins 2^{d-1}-1, appelle GJBExtraction
    à nouveau, pour cet ensemble et pour d-1, et rajoute cet élément dans la colonne d-1 de L[0] un certain nombre de fois
    (dépendant du nombre d'éléments ajoutés pour tous ses chi respectifs) et incrémente compteurs[d-1] de ce nombre.
    Retourne ensuite le nombre de "nouveaux" éléments ajoutés
NOTES:
    - Cette allocation "dynamique" de mémoire pour L est nécessaire afin de ne pas gâcher trop de mémoire, vu que l'on a
    aucune information sur le nombre de GJB contenues dans un ensemble avant leurs calculs.
    - d_const ne sert uniquement qu'à la réallocation de mémoire de L (pour que la fonction sache combien de 
    "colonnes" L possède)
    - Le vecteur compteurs[] sert à savoir où mettre le prochain élément dans L[0], et doit donc pour cela être
    "en dehors" de la fonction afin de ne pas être affecté par la récursivité.
_____________________________________________________________________________*/
uint64_t GJBExtraction(uint64_t* S, int d, uint64_t Taille, uint64_t**pt, uint64_t*** L, int n, uint64_t* compteurs, uint64_t** phi, int d_const){
    uint64_t Taille_L_new = 0;
    uint64_t Taille_L_old = 0;
    uint64_t Taille_S_a;
    phi_d(S, d, phi, n, Taille);
    if((d-1) != 0){
        for (uint64_t i = 0; i < Taille; i++){
            if (phi[d-1][i] == 1){
                Taille_S_a = chi_a_S(S[i], S, pt, Taille, d);
                if (Taille_S_a >= ((1 << (d-1)) - 1)){
                    Taille_L_old = GJBExtraction(pt[d-1], d-1, Taille_S_a, pt, L, n, compteurs, phi, d_const) + Taille_L_new;
                    for(uint64_t j = (Taille_L_new + compteurs[d-1]); j < (Taille_L_old + compteurs[d-1]); j++){
                        L[0][j][d-1]=S[i];
                    }
                    Taille_L_new=Taille_L_old;
                }
            }
        }
        compteurs[d-1] += Taille_L_old;
        return(Taille_L_new);
    }
    if(d == 1){
        for(int l = compteurs[0]; l < (Taille + compteurs[0]); l++){
            if(((l + 1) % 1024 == 0) & ((l + 1) <= (Taille + compteurs[0]))){      /* Réallocation dynamique de mémoire pour L */
                L[0]=(uint64_t **)realloc(L[0],(l+1+1024) * sizeof(uint64_t*));
                for(int q=0; q < 1024; q++){
                    L[0][l+1+q] = malloc(d_const * sizeof(uint64_t));
                }
            }
            L[0][l][0] = S[l - compteurs[0]];
        }
        Taille_L_new = Taille;
        compteurs[0] += Taille;
        return(Taille_L_new);
    }
}

/*
Fonction GJBExtraction définie en Section 3.2, améliorée pour S trié (même structure et remarques que GJBExtraction())
_____________________________________________________________________________*/
uint64_t GJBExtraction_Tri(uint64_t* S, int d, uint64_t Taille, uint64_t**pt, uint64_t*** L, int n, uint64_t* compteurs, uint64_t** phi, int d_const){
    uint64_t Taille_L_new = 0;
    uint64_t Taille_L_old = 0;
    uint64_t Taille_S_a;
    phi_d_Tri(S, d, phi, n, Taille);
    if((d-1) != 0){
        for (uint64_t i = 0; i < Taille; i++){
            if (phi[d-1][i] == 1){
                Taille_S_a = chi_a_S_Tri(S[i], S, pt, Taille, d, n);
                if (Taille_S_a >= ((1 << (d-1)) - 1)){
                    Taille_L_old = GJBExtraction_Tri(pt[d-1], d-1, Taille_S_a, pt, L, n, compteurs, phi, d_const) + Taille_L_new;
                    for(uint64_t j = (Taille_L_new + compteurs[d-1]); j < (Taille_L_old + compteurs[d-1]); j++){
                        L[0][j][d-1]=S[i];
                    }
                    Taille_L_new = Taille_L_old;
                }
            }
        }
        compteurs[d-1] += Taille_L_old;
        return(Taille_L_new);
    }
    if(d == 1){
        for(uint64_t l = compteurs[0]; l < (Taille + compteurs[0]); l++){
            if(((l + 1) % 1024 == 0) & ((l + 1) <= (Taille + compteurs[0]))){      /* Réallocation de mémoire */
                L[0]=(uint64_t **)realloc(L[0], (l+1+1024) * sizeof(uint64_t*));
                for(int q=0; q < 1024; q++){
                    L[0][l+1+q]=malloc(d_const * sizeof(uint64_t));
                }
            }
            L[0][l][0] = S[l-compteurs[0]];
        }
        Taille_L_new = Taille;
        compteurs[0] += Taille;
        return(Taille_L_new);
    }
}

/*_____________________________________________________________________________
EXEMPLE:

S = {001, 010, 100, 110, 111} = {1,2,4,6,7}
CLAIREMENT DEUX SOUS-ESPACES DE DIMENSION 2, {0, 010, 100, 110} ET {0, 001, 110, 111}, DONC DEUX BASES
ON A:
chi_001(S) = {110}, chi_010(S) = {100}, autres chi_a vides
DONC, pour a=001, le code fait:
    L' <- GJBExtraction(chi_001(S), d=1) = GJBExtraction(110, 1) = {110}
    et ajouter {001}u{110} A L
puis, pour a=010, fait:
    L' <- GJBExtraction(chi_010(S), d=1) = GJBExtraction(100, 1) = {100}
    et ajouter {010}u{100} A L
ET DONC OBTENIR L = {{001}{110}; {010}{100}}
_____________________________________________________________________________*/

int main(int argc, char const *argv[]){
    /*_____________________________________________________________________________
    Pas de contraintes sur argc, et on s'attend à ce que: 
    - argv[1] = n; 
    - argv[2] = Taille de S; 
    - argv[3] = d
    _____________________________________________________________________________*/

    /* DÉFINITION DES VALEURS DE n, d, Taille (= nombre d'éléments de S) */
    int n = atoi(argv[1]);
    uint64_t Taille = atoi(argv[2]);
    int d = atoi(argv[3]);

    /* Création de S */
    uint64_t* S = malloc((Taille) * sizeof(uint64_t));

    /* Allocation initiale de mémoire pour L */
    uint64_t*** L = malloc(sizeof(uint64_t**));
    L[0] = malloc(1024 * sizeof(uint64_t*));    
    for(int j = 0; j < 1024; j++){
        L[0][j] = malloc(d * sizeof(uint64_t));
    }
    
    /* Allocation initiale de mémoire pour pt */
    uint64_t** pt = malloc(d * sizeof(uint64_t*));
    for(int i = 0; i < d; i++){
        pt[i] = malloc(Taille * sizeof(uint64_t));
    }
    /* Allocation initiale de mémoire pour phi */
    uint64_t** phi = malloc(d * sizeof(uint64_t*));
    for(int i = 0; i < d; i++){
        phi[i] = malloc(Taille * sizeof(uint64_t));
    }
    /* Création et initialisation du vecteur compteurs */
    uint64_t* compteurs = malloc(d * sizeof(uint64_t));
    for(int i = 0; i < d; i++){
        compteurs[i] = 0;
    }
    
    /*_____________________________________________________________________________
    Ci-dessous, définition des éléments de S. Par défaut ici, S est défini comme l'ensemble {1, ..., Taille}
    Disponible également en dessous en commentaires, une définition de S utilisant la fonction aléatoire de 64 bits. 
    ATTENTION, utiliser cette fonction seulement avec n=64 pour éviter d'obtenir deux éléments égaux dans S
    Si vous décidez de définir un S personnalisé, définissez le ici et vérifiez bien que son nombre d'élément
    est égal à l'argument Taille = argv[2] du main.
    _____________________________________________________________________________*/
    
    for (uint64_t i = 0; i < Taille; i++){
        S[i] = (i+1);
    }
    
    /*
    uint64_t state = 1;
    for (uint64_t i = 0; i < Taille; i++){
        const uint64_t n = rnd64(state++);
        S[i] = n;
    }
    */

    /*_____________________________________________________________________________
    Ci-dessous l'utilisation de GJBExtraction() et l'affichage des bases trouvées. Par défaut la version avec tri
    est utilisée. Si vous désirez utiliser la version sans tri, remplacez GJBExtraction_Tri par GJBExtraction, et 
    commentez les 4 premières lignes suivantes. N'oubliez pas de re-spécifier le type de "sec" et "temps" dans ce cas là
    _____________________________________________________________________________*/
    clock_t sec = clock();
    triRapid(S, 0, Taille-1);
    double temps = ((double)(clock() - sec)/CLOCKS_PER_SEC);
    printf("Le temps de calcul du tri de S est: %f secondes\n", temps);
    

    sec = clock();
    uint64_t x = GJBExtraction_Tri(S, d, Taille, pt, L, n, compteurs, phi, d);
    temps = ((double)(clock() - sec)/CLOCKS_PER_SEC);

    /* Affichage des bases trouvées - ATTENTION: peut tuer votre terminal si vous avez utilisé la
    version S[i] = {1, ..., Taille} avec n et d trop grands */
    for(uint64_t j = 0; j < x; j++){
        printf("Le %llu-ème sous espace de dimension %d est : ",j+1, d);
        Affiche_Ensemble_ord(L[0][j], d, n);
    }
    
    printf("Le temps de calcul de GJBExtraction_Tri est: %f secondes\n", temps);
    printf("Le nombre de GJB de dimension %d dans S vaut %llu\n", d, x);

    /* Remise à 0 des compteurs - inutile si un seul appel, mais si vous décidez d'appeler la fonction GJBExtraction_Tri 
    ou GJBExtraction pluseurs fois d'affilé dans le main, NE JAMAIS OUBLIER DE FAIRE ENTRE CHAQUE APPEL */
    for(int i = 0; i < d; i++){
        compteurs[i] = 0;
    }
    
    /* Libération de mémoire des vecteurs et matrices malloquées au début du main */
    free(S);
    
    free(compteurs);
    
    for(int j = 0; j < ((x / 1024) + 1) * 1024; j++){
        free(L[0][j]);
    }
    free(L[0]);
    free(L);

    for(int i = 0; i < d; i++){
        free(pt[i]);
    }
    free(pt);

    for(int i = 0; i < d; i++){
        free(phi[i]);
    }
    free(phi);
}