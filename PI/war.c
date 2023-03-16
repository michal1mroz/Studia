#include<stdio.h>
#include<stdlib.h>

#define DECK_SIZE 52

// Creating shuffeld deck bassing on the seed inputed
int rand_num(int a, int b){
    int num = rand()%(b-a+1)+a;
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
		// Some error here, but it should't happen >:(
		if(Player->len == DECK_SIZE){
			abort();
		}
    int r = cbuff_norm(Player->len + Player->out);
    Player->deck[r]=val;
    Player->len++;
}
int cbuff_pop(struct cbuff *Player){
    int val = Player->deck[Player->out];
    Player->len--;
    Player->out = cbuff_norm(Player->out+1);
    return val;
}
// Printing all the values stored in the player's deck
void cbuff_print(const struct cbuff *Player){
    int l = Player->len;
    int r = Player->out;
    for(int i=0;i<l;i++){
        printf("%d ",(Player->deck[r]));
        r = cbuff_norm(r+1);
    }
}
// Getting the value of the card under (r) depth into the queue for the card comparison
int cbuff_get_card(const struct cbuff *Player,int r){
	r = cbuff_norm(Player->out + r);
	int val = Player->deck[r];
	return val;
}
// Passing cards from Givers queue to Receivers queue
void cbuff_pass_cards(struct cbuff *Giver, struct cbuff *Receiver,int depth){
	int i;
	for(i=0;i<depth;i++){
		int val = cbuff_pop(Giver);
		cbuff_push(Receiver,val);
	}
}

// Game loop | Oompa-loompa 
void game_loop(int MAX_GAMES,int GAME_TYPE, struct cbuff *Player_A,struct cbuff *Player_B){
	int game_counter = 0; // Counting number of card comparisons
	int cards_dealt = 0; // keeping the depth for which comparisons go into the queue

	while(1){
		// Out of moves
		if(++game_counter>MAX_GAMES){
			printf("0 %d %d %d",Player_A->len,Player_B->len);
			return;
		}
		// Player A wins
		if(Player_B->len==0){
			printf("2 %d",game_counter-1);
			return;
		}
		// Player B wins
		if(Player_A->len==0){
			printf("3\n");
			cbuff_print(Player_B);
			return;
		}
		cards_dealt++;
		// Out of cards
		if(Player_A->len < cards_dealt || Player_B->len < cards_dealt){
			printf("1 %d %d",Player_A->len,Player_B->len);
			return;
		}
		// Getting the first card from the players
		int card_A = cbuff_get_card(Player_A,cards_dealt-1);
		int card_B = cbuff_get_card(Player_B,cards_dealt-1);
		card_A = card_A >> 2; // Same as division by 4. 
		card_B = card_B >> 2; 
		// A wins fight
		if(card_A>card_B){
			cbuff_pass_cards(Player_A,Player_A,cards_dealt);
			cbuff_pass_cards(Player_B,Player_A,cards_dealt);
			cards_dealt = 0;
		}
		else if(card_A==card_B){
			if(GAME_TYPE==0){
				cards_dealt+=1;
			}
			else if(GAME_TYPE==1){
				cbuff_pass_cards(Player_A,Player_A,cards_dealt);
				cbuff_pass_cards(Player_B,Player_B,cards_dealt);
				cards_dealt = 0;
			}
		}// B wins fight
		else{
			cbuff_pass_cards(Player_B,Player_B,cards_dealt);
			cbuff_pass_cards(Player_A,Player_B,cards_dealt);
			cards_dealt = 0;
		}
	}
}

int main(void){

    // Initialization
    int seed;
		int GAME_TYPE,MAX_GAMES;
    scanf("%d",&seed);
		scanf("%d",&GAME_TYPE);
		scanf("%d",&MAX_GAMES);
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
	 	// Running the main game function
    game_loop(MAX_GAMES,GAME_TYPE,&Player_A,&Player_B);
    
    return 0;
}
