/* Inclusion des bibliothèques standard de C */
#include <stdio.h>
#include <stdlib.h>         
#include <time.h>

void permuter(int *a, int *b) {
    *a = (*a)^(*b);
    *b = (*a)^(*b);
    *a =(*a)^(*b);}

void triRapid(int tab[], int first, int last) {
    int pivot, i, j;
    if(first < last) {
        pivot = first;
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




int main(int argc, char const *argv[]){
    int* L=malloc(6*sizeof(int));
    L[0]= 5;
    L[1]= 1;
    L[2]= 2;
    L[3]=7;
    L[4]=2;
    L[5]= 1;
    for(int i=0; i<6;i++){
        printf("%d",L[i]);
    }
    printf("\n");
    triRapid(L,0,5);
    for(int i=0; i<6;i++){
        printf("%d",L[i]);
    }

}