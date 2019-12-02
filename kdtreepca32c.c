/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2018/19
* 
* Progetto dell'algoritmo di Product Quantization for Nearest Neighbor Search
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 kdtreepca32.nasm && gcc -O0 -m32 -msse kdtreepca32.o kdtreepca32c.c -o kdtreepca32c && ./kdtreepca32c
* 
* oppure
* 
* ./runkdtreepca32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

#define	MATRIX		float*

typedef struct
{
 int indP;
 float P;
 int Livello;
 struct KDTREET *figlioSx, *figlioDx;
 }KDTREET;

#define	KDTREE		KDTREET* // modificare con il tipo di dato utilizzato


typedef struct {
    char* filename; //nome del file, estensione .ds per il data set, estensione .qs per l'eventuale query set
    MATRIX ds; //data set 
    MATRIX qs; //query set
    int n; //numero di punti del data set
    int k; //numero di dimensioni del data/query set
    int nq; //numero di punti del query set
    int h; //numero di componenti principali da calcolare 0 se PCA non richiesta
    int kdtree_enabled; //1 per abilitare la costruzione del K-d-Tree, 0 altrimenti
    KDTREE kdtree; //riferimento al K-d-Tree, NULL se costruzione non richiesta
    float r; //raggio di query, -1 se range query non richieste
    int silent; //1 per disabilitare le stampe, 0 altrimenti
    int display; //1 per stampare i risultati, 0 altrimenti

//    MATRIX U; //matrice U restituita dall'algoritmo PCA
//    MATRIX V; //matrice V restituita dall'algoritmo PCA

    //strutture nostre
    float* H; // matrice n*2*k delle regioni dei punti del dataset
    float* U; // matrice degli score n*h
    float* V; // matrice dei load d*h
    float* u; // vettori colonna di D n*1
    float* v; // vettore colonna load k*1
    float* puntoP; // vettore di costruzione punto più vicino
    float* dsTras;   // dataset trasposto temporaneo k*n

    int* vetTmp; // vettore usabile di n posizioni

// strutture per centrare su media
    float* vetMediaDS; // vettore con ogni cella la media di una colonna del ds n*1

// strutture per ricerca range
    float* Point;   // vettore k*1 che in caso di pca viene riempito fino a h*1
    float* Qoint;   // vettore k*1 che in caso di pca viene riempito fino a h*1
    float* qsRidotto; // dataset query ridotto a n*h

    //STRUTTURE OUTPUT MODIFICABILI
    int* QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/


void* get_block(int size, int elements) { 
    return _mm_malloc(elements*size,16); 
}


void free_block(void* p) { 
    _mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(double),rows*cols);
}


void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}


/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;
    
    fp = fopen(filename, "rb");
    
    if (fp == NULL){
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }
    
    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    
    MATRIX data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(float), rows*cols, fp);
    fclose(fp);
    
    *n = rows;
    *k = cols;
    
    return data;
}

/*
* 
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
* 
*/
void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
        fwrite(&n, 4, 1, fp);
        fwrite(&k, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, 4, k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += 4*k;
        }
    }
    fclose(fp);
}



// PROCEDURE ASSEMBLY
extern void prova(params* input);



/*
*	Metodi di supporto
*/

void stampaMatrice(float* m, int r, int c){
    printf("\n");
    int i,j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            if(j==c-1) printf("%.2f", m[i*c+j]);
            else printf("%.2f\t", m[i*c+j]);
        }
        printf("\n");
    }
    
}

void divisioneMatriceScalare( float* m , float s, int nRighe, int nColonne, float* risultato){
    int i,j;
        for(i=0; i<nRighe; i++){
            for(j=0; j<nColonne; j++){
                risultato[i*nColonne+j]= m[ i*nColonne +j ]/s;
            }
        }
}

void trasponi(float* m, int nRighe, int nColonne, float* risultato){
    int i, j;
    for(i=0; i<nRighe; i++){
        for(j=0; j<nColonne; j++){
            risultato[ j*nRighe+i ] = m[ i*nColonne+j ];
        }
    }
}

void prodMatrVett(float* m, float* v, int nRighe, int nColonne, int lengthVettore, float* risultato){
    if( nColonne==lengthVettore){
        int i,j;
        float sum;
        for(i=0; i<nRighe; i++){
            sum=0.0;
            for(j=0; j<nColonne; j++){
                sum += m[i*nColonne+j] * v[j];
            }
            risultato[i] = sum;
        }
    }
}

float prodScalare(float* v1, int dim1, float* v2, int dim2){
    int i;
    float r=0.0;

    if(dim1=dim2){
        for(i=0; i<dim1; i++)
            r+= v1[i]*v2[i];
    
    return r;
    }    
}

float calcolaNorma( float* v, int dim){
    float r=0;
    int i;
    for(i=0; i<dim;i++){
        r+= powf(v[i],2.0);
    }
    return sqrtf(r);
}

void sottrazioneMatrici(float* m1, float* m2, int nRighe1, int nColonne1, int nRighe2, int nColonne2, float* risultato){

    if (nRighe1==nRighe2 && nColonne1==nColonne2){
        int i,j;

        for(i=0; i<nRighe1; i++){
            for(j=0; j<nColonne1; j++){
                risultato[ i*nColonne1+j ]= m1[i*nColonne1 +j] - m2[i*nColonne1 +j] ;
            }
        }
    }
}

void prodMatrici(float* m1, int nRighe1, int nColonne1, float* m2, int nRighe2, int nColonne2, float* risultato){

    if(nColonne1==nRighe2){
        int i,j,l,z;
        float sum;
        for(i=0; i<nRighe1; i++){
            for(j=0; j<nColonne2; j++){
                sum=0.0;
                for(l=0; l<nColonne1; l++){
                    for(z=0; z<nRighe2; z++){
                        sum += m1[i*nColonne1+l] * m2[z*nColonne2+j];
                    }
                }
                risultato[i*nColonne2+j] = sum;
            }
        }
    
    }
}

void prodVettori(float* v1, int dim1, float* v2, int dim2, float* risultato){
    int i,j;
    for (i = 0; i < dim1; i++){
        for(j=0; j<dim2; j++){
            risultato[ i*dim2 + j ] = v1[i] * v2[j];
        }
    }
    
}

void divisioneVettoreScalare(float* v, float s, int dim, float* risultato){
    int i;
    for(i=0; i<dim; i++){
        risultato[i] = v[i]/s;
    }
}

float distanzaEuclidea(float* P,float* Q, int dimen){
    int i;
    float dist=0;

    for(i=0; i<dimen; i++){
        dist+= powf(P[i]-Q[i],2);    
    }
    return sqrtf(dist);

}


/*
*   ===================================================================================================================
*	PCA
* 	===================================================================================================================
*/

void centraMediaDS(params* input){
    int dimensione, punto;
    float media;
    for(dimensione=0; dimensione<input->k; dimensione++){
        media=0.0;
        for(punto=0; punto<input->n; punto++){
            media += input->ds[punto*input->k+dimensione];
        }
        media /= input->n;
        for(punto=0; punto<input->n; punto++){
            input->ds[punto*input->k+dimensione] -= media;
        }
        input->vetMediaDS[dimensione] = media;
    }
}


void nipals(params *input){
    float soglia= 1*expf(-8);
    centraMediaDS(input);
    int i;
    for(i=0; i<input->n; i++){
        input->u[i]=input->ds[input->k*i]; //ci prendiamo gli elementi della colonna 0
    }
    for(i=0; i<input->h; i++){
        start_for: 
        trasponi(input->ds, input->n, input->k, input->dsTras);
        float t= prodScalare(input->u, input->n, input->u, input->n);
        prodMatrVett(input->dsTras, input->u, input->k, input->n, input->n, input->v);
        divisioneVettoreScalare(input->v, t, input->k, input->v);

        divisioneVettoreScalare(input->v, calcolaNorma(input->v, input->k), input->k, input->v);

        prodMatrVett(input->ds, input->v, input->n, input->k, input->k, input->u);
        divisioneVettoreScalare(input->u, prodScalare(input->v, input->k,input->v, input->k), input->n, input->u);

        float t1= prodScalare(input->u, input->n, input->u, input->n);

        if( fabsf( t1 - t) >= soglia*t1)
            goto start_for;

        else{
            int j;
            for(j=0; j<input->n; j++){
                input->U[j*input->h+i]=input->u[j];
            }
            for(j=0; j<input->k; j++){
                input->V[j*input->h+i]=input->v[j];
            }

            prodVettori(input->u, input->n, input->v, input->k, input->dsTras);
            sottrazioneMatrici(input->ds, input->dsTras, input->n, input->k, input->n, input->k, input->ds);
            //stampaMatrice(input->ds, input->n, input->k);
            //scanf("%d",&j);
            //printf("%d\t",i);
        }
    }

}

void pca(params* input) {
    // Codificare qui l'algoritmo PCA
    //prova(input);
    // Calcola le matrici U e V
    nipals(input);
}

/*
/*
* 	===================================================================================================================
* 	===================================================================================================================
*	K-d-Tree
* 	===================================================================================================================
* 	===================================================================================================================
*/

void merge(float* arr, int l, int m, int r, params* input) 
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    /* create temp arrays */
    float L[n1], R[n2]; 
    int Lind[n1], Rind[n2];
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){
        L[i] = arr[l + i];
        Lind[i] = input->vetTmp[l+i]; 
    }
    for (j = 0; j < n2; j++) {
        R[j] = arr[m + 1+ j];
        Rind[j] = input->vetTmp[m+1+j]; 
    }
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) { 
        if (L[i] <= R[j]) { 
            arr[k] = L[i];
            input->vetTmp[k] = Lind[i];
            i++; 
        } 
        else{ 
            arr[k] = R[j]; 
            input->vetTmp[k] = Lind[j];
            j++; 
        } 
        k++; 
    } 
  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) { 
        arr[k] = L[i]; 
        input->vetTmp[k] = Lind[i];
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) { 
        arr[k] = R[j];
        input->vetTmp[k] = Lind[j];
        j++; 
        k++; 
    } 
} 

void mergeSort(float* arr, int l, int r, params* input){
    if (l < r){
    // Finding mid element
    int m = l+(r-l)/2;
    // Recursively sorting both the halves
    mergeSort(arr, l, m, input);
    mergeSort(arr, m+1, r, input);

    // Merge the array
    merge(arr, l, m, r, input);
    }
}


int cercaMediano(float* dataset, int dimensioneTaglio, int nCol, params* input){
    int i;
    float* arr = (float*) malloc(sizeof(float)*input->n);
    for(i=0; i<input->n; i++){
        arr[i] = dataset[i*nCol+dimensioneTaglio];
        input->vetTmp[i] = i;
    }
    mergeSort(arr, 0, input->n-1, input);
    
    if(input->n%2!=0){
        return input->vetTmp[(input->n+1)/2];
    }else if(input->n%2==0){
        return input->vetTmp[((input->n/2))];
    }
}

float** creaDataset(params *input, float* dataset, int col, int c, int indP){
    int i;
    int j1=0,j2=0;
    int z;
    float **ris = (float **) malloc( 2 * sizeof(float*));
    float* dsMinore= (float* ) malloc((input->n/2 * col) *sizeof(float));
    float* dsMaggiore= (float* ) malloc((input->n/2 * col) *sizeof(float));    

    ris[0]= dsMinore;
    ris[1]= dsMaggiore;

    for(i=0; i<input->n; i++){
        if(i!= indP){

            if (dataset[i * col +c] <  dataset[indP * col + c]){

                for(z=0; z< col; z++){
                    dsMinore[j1*col + z] =dataset[i*col +z];
                }
                j1++;
            }
            else if (dataset[i * col +c] >=  dataset[indP * col + c]){
                for(z=0; z< col; z++){
                    dsMaggiore[j2*col + z] =dataset[i*col +z];
                }
                j2++;
            }
        }
        else if( input-> n%2==0){
               for(z=0; z< col; z++){
                    dsMinore[j1*col + z] =dataset[i*col +z];
                }
                j1++;
        }
    }
    return ris;
}

KDTREE buildTree(MATRIX dataset,int livello, int col, params *input){
    if( dataset == NULL) return NULL;
    int c= livello%col;
    int indicePunto= cercaMediano(dataset,c,col,input);

    float** dueDataset= creaDataset( input, dataset, col, c , indicePunto);

    //MATRIX D1=creaDataset(dataset, c, indicePunto,0);
    //MATRIX D2=creaDataset(dataset, c, indicePunto,1);
    KDTREE n= malloc(sizeof(KDTREET));
    n->indP=indicePunto;
    n->P= dataset[indicePunto*col+c];   
    n->figlioSx= buildTree(dueDataset[0], livello+1, col, input);
    n->figlioDx= buildTree(dueDataset[1], livello+1, col, input);
    return n;
}

void kdtree(params* input) {
    
    // -------------------------------------------------
    // Codificare qui l'algoritmo di costruzione
    // -------------------------------------------------

    if(input->h > 0){
        input->kdtree= buildTree(input->U, 0, input->h, input);
    }else{    
        input->kdtree= buildTree(input->ds, 0, input->k, input);
    }
}


/*
* 	===================================================================================================================
* 	===================================================================================================================
* 	===================================================================================================================
*	Range Query Search
* 	===================================================================================================================
* 	===================================================================================================================
* 	===================================================================================================================
*/

float distance( int indQ, int indP, params* input) {
    
    int nColonne = input->h>0? input->h : input->k;

    int j;
    for(j=0; j< nColonne; j++){
        input->Qoint[j]=input->qs[indQ * input->k +j];
        if(input->qs[indQ*nColonne +j] <= input->H[ indP * nColonne *2 + j*2])
            input->Point[j]= input->H[ indP * input->k *2 + j*2];
        else if ( input->qs[indQ*input->k +j] >= input->H[ indP * input->k *2 + j*2 +1])
            input->Point[j]= input->H[ indP * input->k *2 + j*2 +1];
        else 

            input->Point[j]= input->qs[indQ * input->k +j];        
    }
    return distanzaEuclidea(input->Point, input->Qoint, input->k);
}


void ricercaRange(float* dataSet, float* querySet, int nColonne, KDTREE n, int indQ, params* input){

    if( distance(indQ, n->indP, input) > input->r) return 0;

    int i;

    for(i=0; i<nColonne; i++){
        input-> Point[i]= dataSet[n->indP*nColonne+ i];
        input-> Qoint[i]= querySet[indQ*nColonne+ i];
    }
    
    if( distanzaEuclidea(input-> Point, input-> Qoint, nColonne) <= input->r ){
        input->QA[input->nQA * 2]=indQ;
        input->QA[input->nQA * 2 +1]=n->indP;
        input->nQA+=1;
    }
    if( n->figlioSx != NULL){
        ricercaRange(dataSet, querySet, nColonne, n->figlioSx, indQ, input);
    }
    if( n->figlioDx != NULL){
        ricercaRange(dataSet, querySet, nColonne, n->figlioDx, indQ, input);
    }

}

void centraMediaQS(params* input){
    int dimensione, punto;
    for(dimensione=0; dimensione<input->k; dimensione++){
        for(punto=0; punto<input->nq;punto++){
            input->qs[punto*input->k+dimensione] -= input->vetMediaDS[dimensione];
        }
    }
}

void range_query(params* input) {
    
    // Codificare qui l'algoritmo di ricerca
    if(input->h>0){
        centraMediaQS(input);
        prodMatrici(input->qs, input->nq, input->k, input->V, input->k, input->h, input->qsRidotto);

    // Calcola il risultato come una matrice di nQA coppie di interi
    // (id_query, id_vicino)
    // o in altro formato
        int i;

        for(i=0; i<input->nq; i++){
            ricercaRange(input->U, input->qsRidotto, input->h, input->kdtree, i, input);
        }
    }else{
        int i;

        for(i=0; i<input->nq; i++){
            ricercaRange(input->ds, input->qs, input->k, input->kdtree, i, input);
        }
    }
}



// 0 Hmin , 1 Hmax, nP numerorigaPunto==(indiceDataset), 
float distanza(float* Q, int dimen, int nP, params* input){
    int j;
    

    for(j=0; j<dimen; j++){
        if(Q[j]<= input-> H[nP*(2*dimen)+j+0]){ // prendo min della coordinata della regione
            input->puntoP[j]= input-> H[nP*(2*dimen)+j+0];
        }

        else if (Q[j]>= input-> H[nP*(2*dimen)+j+1]){
            input->puntoP[j]= input-> H[nP*(2*dimen)+j+1];
        }
        else
        {
            input->puntoP[j]=Q[j];
        }         
    }
    return distanzaEuclidea(input->puntoP, Q, dimen);
}

//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************
//******************************************************************************************************************

int main(int argc, char** argv) {
    
    char fname[256];
    int i, j, k;
    clock_t t;
    float time;
    
    //
    // Imposta i valori di default dei parametri
    //
    
    params* input = malloc(sizeof(params));
    
    input->filename = NULL;
    input->h = 0;
    input->kdtree = NULL;
    input->r = -1;
    input->silent = 0;
    input->display = 1;
    input->QA = NULL;
    input->nQA = 0;
    
    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //
    
    if (argc <= 1 && !input->silent) {
        printf("Usage: %s <data_name> [-pca <h>] [-kdtree [-rq <r>]]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\t-d: display query results\n");
        printf("\t-s: silent\n");
        printf("\t-pca <h>: h-component PCA enabled\n");
        printf("\t-kdtree: kdtree building enabled\n");
        printf("\t-rq <r>: range query search with radius r enabled\n");
        printf("\n");
        exit(0);
    }
    
    //
    // Legge i valori dei parametri da riga comandi
    //
    
    int par = 1;
    while (par < argc) {
        if (par == 1) {
            input->filename = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-pca") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing h value!\n");
                exit(1);
            }
            input->h = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-kdtree") == 0) {
            input->kdtree_enabled = 1;
            par++;
            if (par < argc && strcmp(argv[par],"-rq") == 0) {
                par++;
                if (par >= argc) {
                    printf("Missing radius value!\n");
                    exit(1);
                }
                input->r = atof(argv[par]);
                if(input->r < 0){
                    printf("Range query radius must be non-negative!\n");
                    exit(1);
                }
                par++;
            }
        } else{
            printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
            par++;
        }
    }

    //
    // Legge i dati e verifica la correttezza dei parametri
    //
    
    if(input->filename == NULL || strlen(input->filename) == 0){
        printf("Missing input file name!\n");
        exit(1);
    }
    
    sprintf(fname, "%s.ds", input->filename);
    input->ds = load_data(fname, &input->n, &input->k);


// AGGIUNTO IO ****************************************************
// AGGIUNTO IO ****************************************************
// AGGIUNTO IO ****************************************************
// AGGIUNTO IO ****************************************************
// AGGIUNTO IO ****************************************************
    input->n = 80;
    input->k = 3;
    input->nq = 20;

    if(input->h < 0){
        printf("Invalid value of PCA parameter h!\n");
        exit(1);
    }
    if(input->h > input->k){
        printf("Value of PCA parameter h exceeds data set dimensions!\n");
        exit(1);
    }
    
    if(input->r >= 0){
        sprintf(fname, "%s.qs", input->filename);
        int k;
        input->qs = load_data(fname, &input->nq, &k);
        if(input->k != k){
            printf("Data set dimensions and query set dimensions are not compatible!\n");
            exit(1);
        }
    }
    
    //
    // Visualizza il valore dei parametri
    //
    
    if(!input->silent){
        printf("Input file name: '%s'\n", input->filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [k]: %d\n", input->k);
        if(input->h > 0){
            printf("PCA search enabled\n");
            printf("Number of principal components [h]: %i\n",input->h);
        }else{
            printf("PCA search disabled\n");
        }
        if(input->kdtree_enabled){
            printf("Kdtree building enabled\n");
            if(input->r >= 0){
                printf("Range query search enabled\n");
                printf("Range query search radius [r]: %f\n",input->r);
            }else{
                printf("Range query search disabled\n");
            }
        }else{
            printf("Kdtree building disabled\n");
        }
    }

    //allocazione strutture dati
    input-> H= (float*) malloc(input->n*input->k*2 * sizeof(float));
    input-> U= (float*) malloc(input->n*input->h * sizeof(float));
    input-> V= (float*) malloc(input->k*input->h * sizeof(float));
    input-> u= (float*) malloc(input->n * sizeof(float));
    input-> v= (float*) malloc(input->k* sizeof(float));
    input-> puntoP= (float*) malloc(input->k * sizeof(float));
    input-> Point= (float*) malloc(sizeof(float)*input->k);
    input-> Qoint= (float*) malloc(sizeof(float)*input->k);
    input-> vetMediaDS = (float*) malloc(sizeof(float) * input-> k);
    input-> dsTras = (float*) malloc( sizeof(float) * input->n * input-> k );
    input-> qsRidotto = (float*) malloc( sizeof(float) * input->nq * input-> h );
    
    input-> vetTmp = (int*) malloc(sizeof(int)*input->n);

    //
    // Calcolo PCA
    //
    
    if(input->h > 0){
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
        sprintf(fname, "%s.U", input->filename);
        save_data(fname, input->U, input->n, input->h);
        sprintf(fname, "%s.V", input->filename);
        save_data(fname, input->V, input->k, input->h);
    }else
        time = -1;
       
    if (!input->silent)
        printf("\nPCA time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    
//    stampaMatrice(input->U, input->n, input->h);
//        printf("\n----------------------------------------------------------\n");
//    stampaMatrice(input->V, input->k, input->h);

//    int x;
//    printf(">");
//    scanf("%d",&x);

    //
    // Costruzione K-d-Tree
    //
    
//AGGIUNTO IO******************************************************************************************
//stampaMatrice(input->U, input->n, input->h);

//int indP = cercaMediano(input->U, 0, input->h, input);
//printf("%d\n",indP);
float** punt = creaDataset(input, input->U, input->h, 0, 75);
float* d1 = punt[0];    // n/2 * h
float* d2 = punt[1];  
  
stampaMatrice(d1,input->n/2,input->h);
stampaMatrice(d2,input->n/2,input->h);  
  
  
    if(input->kdtree){
        t = clock();
        kdtree(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
        time = -1;
    if (!input->silent)
        printf("\nIndexing time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Range query search
    //
    
    if(input->r >= 0){
        t = clock();
        range_query(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
        time = -1;
    if (!input->silent)
        printf("\nQuerying time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    //
    // Salva il risultato delle query
    // da modificare se si modifica il formato delle matrici di output
    //
    
    if(input->r >= 0){
        if(!input->silent && input->display) {
            //NB: il codice non assume che QA sia ordinata per query, in caso lo sia ottimizzare il codice
            printf("\nQuery Answer:\n");
            for(i = 0; i < input->nq; i++){
                printf("query %d: [ ", i);
                for(j = 0; j < input->nQA; j++)
                    if(input->QA[j*2] == i)
                        printf("%d ", input->QA[j*2+1]);
                printf("]\n");
            }
            printf("\n");
        }
        sprintf(fname, "%s.qa", input->filename);
        save_data(fname, input->QA, input->nQA, 2);
    }
    
    if (!input->silent)
        printf("\nDone.\n");

    return 0;
}