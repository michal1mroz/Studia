// Need to add the gameplay mechanic and all the logic

#include<stdio.h>
#include<stdlib.h>

#define DECK_SIZE 52

// Creating shuffeld deck bassing on the seed inputed
int rand_num(int a, int b){
    int num = rand()%(b-a+1)+1;
    return num;
}
void shuffle(int tab[],int seed){
    srand((unsigned) seed);
    int i;
    for(i=0;i<DECK_SIZE;i++){
        tab[i]=i;
    }
    for(i=0;i<DECK_SIZE-1;i++){
        int k = rand_num(i,DECK_SIZE-1);
        int tmp = tab[i];
        tab[i]=tab[k];
        tab[k]=tmp;
    }
}
// Queue deffinitions and functions
// With struct there is no need to remember the sizes of both queues
struct cbuff {
    int len;
    int out;
    int deck[DECK_SIZE];
};
int cbuff_norm(int r){
    if(r>=DECK_SIZE){r-=DECK_SIZE;}
    if(r<0){r+=DECK_SIZE;}
    return r;
}
void cbuff_push(struct cbuff *Player, int val){
    int r = cbuff_norm(Player->len+Player->out);
    Player ->deck[r]=val;
    Player->len++;
}
int cbuff_pop(struct cbuff *Player, int r){
    int val = Player->deck[Player->out];
    Player->len--;
    Player->out = cbuff_norm(Player->out+1);
    return val;
}
void cbuff_print(const struct cbuff *Player){
    int i;
    int r = Player->out;
    for(i=0;i<(Player->len);i++){
        printf("%d ",(Player->deck[i]));
        r = cbuff_norm(r+1);
    }
}

// Game stuff


int main(void){

    // Initialization
    int seed = 0;
    scanf("%d",&seed);
    int deck[DECK_SIZE];
    struct cbuff Player_A = {.len = 0, .out = 0};
    struct cbuff Player_B = {.len =0, .out = 0};
    shuffle(deck,seed);
    int i;
    // Dealing cards to both players
    for(i=0;i<DECK_SIZE/2;i++){
        cbuff_push(&Player_A ,deck[i]);
    }
    for(i;i<DECK_SIZE;i++){
        cbuff_push(&Player_B , deck[i]);
    }

    /*
    cbuff_print(&Player_A);
    printf("\n");
    cbuff_print(&Player_B);  
    */

    return 0;
}
