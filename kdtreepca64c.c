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
#include <libgen.h>
#include <xmmintrin.h>

#define	MATRIX		float*
//#define	KDTREE		float* // modificare con il tipo di dato utilizzato



typedef struct KDTREE {
    float *P;
    struct KDTREE *figlioSx;
    struct KDTREE *figlioDx;
    float *H;

 }KDTREE;


typedef struct {
   char* filename; //nome del file, estensione .ds per il data set, estensione .qs per l'eventuale query set
    MATRIX ds; //data set 
    MATRIX qs; //query set
    int n; //numero di punti del data set
    int k; //numero di dimensioni del data/query set
    int nq; //numero di punti del query set
    int h; //numero di componenti principali da calcolare 0 se PCA non richiesta
    int kdtree_enabled; //1 per abilitare la costruzione del K-d-Tree, 0 altrimenti
    KDTREE *kdtreeRoot; //riferimento al K-d-Tree, NULL se costruzione non richiesta
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
    float** pQA; //punto in posizione nQA associato al punto query in posizione nQA di QA
    int nQA; //numero di risposte alle query
    int *QAfin;
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
    return _mm_malloc(elements*size,32); 
}


void free_block(void* p) { 
    _mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(float),rows*cols);
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
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, 4, k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += 4*k;
        }
    }
    fclose(fp);
}


// PROCEDURE ASSEMBLY
extern void prova(float* y, float* x, float* ris);
extern void calcolaNorma_ass_64(float* arr, int length, float* norma);
//extern void prodMatrVett_ass_64(float* m, float* v, int nRighe, int nColonne, int lengthVettore, float* risultato);
extern void prodScalare_ass_64(float* v1,int dim1,float* v2,int dim2,float* rs);
extern void divisioneVettoreScalare_ass_64(float* v, float* s, int dim, float* risultato);
//extern void sottrazioneMatrici_ass_64(float* m1, float* m2, int nRighe1, int nColonne1, int nRighe2, int nColonne2, float* risultato);
extern void distanzaEuclidea_ass_64(float* P,float* Q, int dimen, float* dist);
//extern void prodMatr_ass_64(float* m1, int nRighe1, int nColonne1, float* m2, int nRighe2, int nColonne2, float* risultato);
//extern void trasponi_ass_64(float* m, int nRighe, int nColonne, float* risultato);



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

void stampaMatriceInt(int* m, int r, int c){
    printf("\n");
    int i,j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            if(j==c-1) printf("%d", m[i*c+j]);
            else printf("%d ", m[i*c+j]);
        }
        printf("\n");
    }
    
}


/*
*	Metodi di supporto
*/

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

void prodScalare(float* v1, int dim1, float* v2, int dim2, float *r){
    int i;
    *r=0.0;

    if(dim1==dim2){
        for(i=0; i<dim1; i++)
            *r+= v1[i]*v2[i];     
    }    
}

void calcolaNorma( float* v, int dim, float *r){
    *r=0;
    int i;
    for(i=0; i<dim;i++){
        *r+= powf(v[i],2.0);
    }
    *r=sqrtf(*r);
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
                    //for(z=0; z<nRighe2; z++){
                        sum += m1[i*nColonne1+l] * m2[l*nColonne2+j];
                    //}
                }
                risultato[i*nColonne2+j] = sum;
            }
        }
    
    }
}

void divisioneVettoreScalare(float* v, float s, int dim, float* risultato){
    int i;
    for(i=0; i<dim; i++){
        risultato[i] = v[i]/s;
    }
}

void distanzaEuclidea(float* P,float* Q, int dimen, float* dist){
    int i;
    *dist=0;
    for(i=0; i<dimen; i++){
        *dist+= powf(P[i]-Q[i],2);    
    }
    *dist=sqrtf(*dist);

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
        float t,t1;
        float r;
        do{
            trasponi(input->ds, input->n, input->k, input->dsTras);
            //trasponi_ass_64(input->ds, input->n, input->k, input->dsTras);

            prodScalare(input->u, input->n, input->u, input->n, &t);
            //prodScalare_ass_64(input->u, input->n, input->u, input->n, &t);

            prodMatrVett(input->dsTras, input->u, input->k, input->n, input->n, input->v);
            //prodMatrVett_ass_64(input->dsTras, input->u, input->k, input->n, input->n, input->v);

            divisioneVettoreScalare(input->v, t, input->k, input->v);
            //divisioneVettoreScalare_ass_64(input->v, &t, input->k, input->v);

            calcolaNorma(input->v, input->k, &r);
            //calcolaNorma_ass_64(input->v, input->k, &r);

            divisioneVettoreScalare(input->v, r, input->k, input->v);
            //divisioneVettoreScalare_ass_64(input->v, &r, input->k, input->v);

            prodMatrVett(input->ds, input->v, input->n, input->k, input->k, input->u);
            //prodMatrVett_ass_64(input->ds, input->v, input->n, input->k, input->k, input->u);

            prodScalare(input->v, input->k,input->v, input->k, &r);
            //prodScalare_ass_64(input->v, input->k,input->v, input->k, &r);
            
            divisioneVettoreScalare(input->u, r, input->n, input->u);
            //divisioneVettoreScalare_ass_64(input->u, &r, input->n, input->u);

            prodScalare(input->u, input->n, input->u, input->n, &t1);
            //prodScalare_ass_64(input->u, input->n, input->u, input->n, &t1);
            
        }while( fabsf( t1 - t) >= soglia*t1);
 
        int j;
        for(j=0; j<input->n; j++){
            input->U[j*input->h+i]=input->u[j];
        }
        for(j=0; j<input->k; j++){
            input->V[j*input->h+i]=input->v[j];
        }

        prodMatrici(input->u, input->n, 1, input->v, 1, input->k, input->dsTras);
        //prodMatr_ass_64(input->u, input->n, 1, input->v, 1, input->k, input->dsTras);

        sottrazioneMatrici(input->ds, input->dsTras, input->n, input->k, input->n, input->k, input->ds);
        //sottrazioneMatrici_ass_64(input->ds, input->dsTras, input->n, input->k, input->n, input->k, input->ds);
        
        }
    }

void pca(params* input) {
    // Codificare qui l'algoritmo PCA
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

void merge(float* a, int nElem, int p, int q, int r,params* input) {
    
    int i, j, k=0; 
    
    float b[nElem];
    int index[nElem];

    i = p;
    j = q+1;
    
    while (i<=q && j<=r) {
        if (a[i]<a[j]) {
            b[k] = a[i];
            index[k] = input->vetTmp[i];
            i++;
        } else {
            b[k] = a[j];
            index[k] = input->vetTmp[j];            
            j++;
        }
        k++;
    }
    while (i <= q) {
        b[k] = a[i];
        index[k] = input->vetTmp[i];
        i++;
        k++;
    }
    while (j <= r) {
        b[k] = a[j];
        index[k] = input->vetTmp[j];
        j++;
        k++;
    }
    for (k=p; k<=r; k++){
        a[k] = b[k-p];
        input->vetTmp[k] = index[k-p];
    }
    return; 
}

void mergeSort(float* a, int nElem, int p, int r, params* input) {
    int q;
    if (p < r) {
        q = (p+r)/2;
        mergeSort(a, nElem, p, q, input);
        mergeSort(a, nElem, q+1, r, input);
        merge(a, nElem, p, q, r, input);
    }
    return;
}


int cercaMediano(float* dataset, int nElem, int dimensioneTaglio,  int nCol, float* minEmax, params* input){
    int i;
    float* arr = (float*) malloc(sizeof(float)*nElem);
    for(i=0; i<nElem; i++){
        arr[i] = dataset[i*nCol+dimensioneTaglio];
        input->vetTmp[i] = i;
    }
    mergeSort(arr, nElem, 0, nElem-1, input);
    minEmax[0]= arr[0];
    minEmax[1]= arr[nElem-1];

    free(arr);
    
    if(nElem%2!=0){
        return input->vetTmp[ ((nElem+1)/2) -1 ];
    }else if(nElem%2==0){
        return input->vetTmp[(nElem/2)-1];
    }
}


float* creaDatasetMaggiore(params* input, float* dataset, int nElem, int nDim, int col, int c, int indP){
    int i;
    int j2=0;
    int z;

    float* dsMaggiore= (float* ) malloc((nDim * col) *sizeof(float));    

     for(i=0; i<nElem; i++){
        if(i!= indP){

            if (dataset[i * col +c] >=  dataset[indP * col + c]){
                for(z=0; z< col; z++){
                    dsMaggiore[j2 * col + z] =dataset[i * col + z];
                }
                j2++;
            }
        }
    }
    return dsMaggiore;    
}


float* creaDatasetMinore(params* input, float* dataset, int nElem, int nDim, int col, int c, int indP){
        int i;
    int j1=0;
    int z;
 
    float* dsMinore= (float* ) malloc((nDim * col) *sizeof(float));  

    for(i=0; i<nElem; i++){
        if(i!= indP){

            if (dataset[i * col +c] <  dataset[indP * col + c]){

                for(z=0; z< col; z++){
                    dsMinore[j1*col + z] =dataset[i*col +z];
                }
                j1++;
            }
        }
    }

    return dsMinore;
}

void calcolaDimDataset(float* dataset, int nElem, int col, int c, int mediano, int *dimMin, int *dimMagg){
    int i;
    for(i=0; i<nElem; i++){
        if(i != mediano){
            if (dataset[ i * col +c] < dataset[mediano*col + c]) *dimMin+=1;
            else *dimMagg+=1;
        }
    }
}
    
KDTREE* buildTree(float* dataset,int nElem, int livello, int col, params *input){

    if( nElem == 0) return NULL;
    int c= livello%col;
    int indicePunto;
    int i,j;
    KDTREE *curr = (KDTREE *) malloc (sizeof(KDTREE));
    curr->H = malloc (sizeof(float)* 2*col);


    float* minEmax= malloc (sizeof(float)*2); 

    for(i=0; i< col; i++){
        int med= cercaMediano(dataset, nElem, i, col,minEmax, input);        
        if( i== c){
            indicePunto=med;
        }
        curr->H[i*2] = minEmax[0];
        curr->H[i*2+1]= minEmax[1];
    }

    int dimMin=0, dimMagg=0;
    calcolaDimDataset(dataset, nElem, col, c, indicePunto, &dimMin, &dimMagg);

    float* datasetMagg = creaDatasetMaggiore(input,dataset,nElem, dimMagg, col,c,indicePunto);
    float* datasetMin = creaDatasetMinore(input,dataset,nElem, dimMin, col,c,indicePunto);

    curr->P = malloc(sizeof(float)); 
    curr->P = &dataset[indicePunto*col];


    curr->figlioSx = buildTree(datasetMin, dimMin, livello+1, col, input);
    curr->figlioDx = buildTree(datasetMagg, dimMagg, livello+1, col, input);

    return curr;
}

void kdtree(params* input) {
    // -------------------------------------------------
    // Codificare qui l'algoritmo di costruzione
    // -------------------------------------------------
    if(input->h > 0){
        input->kdtreeRoot= buildTree(input->U, input->n, 0, input->h, input);
    }else{    
        input->kdtreeRoot= buildTree(input->ds, input->n, 0, input->k, input);
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

void distance( float* querySet, int nColonne, int indQ, KDTREE *nodo, params* input, float* dist) {
    
    int j;

    for(j=0; j< nColonne; j++){
        input->Qoint[j]=querySet[indQ * nColonne +j];
        
        if(querySet[indQ*nColonne +j] <= nodo->H[j*2])
            input->Point[j]= nodo->H[j*2];
        
        else if ( querySet[indQ*nColonne +j] >= nodo->H[j*2 +1])
            input->Point[j]= nodo->H[j*2 +1];
        
        else 
            input->Point[j]= querySet[indQ * nColonne +j];        
    }
    distanzaEuclidea(input->Point, input->Qoint, nColonne, dist);
    //distanzaEuclidea_ass_64(input->Point, input->Qoint, nColonne, dist);
}


void ricercaRange(float* dataSet, float* querySet, int nColonne, KDTREE *n, int indQ, params* input){
    float dist;
    distance(querySet, nColonne, indQ, n, input, &dist);
    if( dist <= input->r){
        int i;

        for(i=0; i<nColonne; i++){
            //input-> Point[i]= dataSet[n->indP*nColonne+ i];
            input-> Qoint[i]= querySet[indQ*nColonne+ i];
        }


        distanzaEuclidea(n->P, input-> Qoint, nColonne, &dist);
        //distanzaEuclidea_ass_64(n->P, input-> Qoint, nColonne, &dist);

        if(  dist <= input->r ){
            
            input->QA[input->nQA]=indQ;
            
            input->pQA[input->nQA] = n->P;
            
            input->nQA+=1;
        }
        if( n->figlioSx != NULL){
            ricercaRange(dataSet, querySet, nColonne, n->figlioSx, indQ, input);
        }
        if( n->figlioDx != NULL){

            ricercaRange(dataSet, querySet, nColonne, n->figlioDx, indQ, input);
        }
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
        //prodMatr_ass_64(input->qs, input->nq, input->k, input->V, input->k, input->h, input->qsRidotto);

    // Calcola il risultato come una matrice di nQA coppie di interi
    // (id_query, id_vicino)
    // o in altro formato
        int i;

        for(i=0; i<input->nq; i++){
            ricercaRange(input->U, input->qsRidotto, input->h, input->kdtreeRoot, i, input);
        }
    }else{
        int i;

        for(i=0; i<input->nq; i++){
            ricercaRange(input->ds, input->qs, input->k, input->kdtreeRoot, i, input);
        }
    }
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


int main(int argc, char** argv) {
    
    int nCol1= 27,nCol2= 27;
    int nRig1= 1,nRig2= 1;

    float scalare = 2.2;

    float* p = (float*)malloc(sizeof(float)*nRig1*nCol1);
    float* p1 = (float*)malloc(sizeof(float)*nRig2*nCol2);

    for(int i=0; i<nRig1; i++){
        for(int j=0; j<nCol1;j++){
            p[i*nCol1+j] = 5.0 + j;
        }
    }

    for(int i=0; i<nRig2; i++){
        for(int j=0; j<nCol2;j++){
            p1[i*nCol1+j] = 5.0;
        }
    }

    stampaMatrice(p,nRig1,nCol1);
    //stampaMatrice(p1,nRig2,nCol2);

    float* ris_1 = (float*)malloc(sizeof(float)*nRig1*nCol1);
    float* ris_2 = (float*)malloc(sizeof(float)*nRig1*nCol1);

    float* punt1  = (float*)malloc(sizeof(float));
    float* punt2  = (float*)malloc(sizeof(float));

    divisioneVettoreScalare(p, scalare, nCol1, ris_1);
    divisioneVettoreScalare_ass_64(p, &scalare, nCol1, ris_2);

    printf("\n************************************\n");
    printf( "C: %.2f", 0.0 );
    stampaMatrice(ris_1, nRig1, nCol1);
    printf( "ASS: %.2f", 0.0 );
    stampaMatrice(ris_2, nRig1, nCol1);
    printf("\n**********************************\n");






    char fname[256];
    char* dsname;
    int i, j, k;
    clock_t t;
    float time;
    
    //
    // Imposta i valori di default dei parametri
    //
    
    params* input = malloc(sizeof(params));
    
    input->filename = NULL;
    input->h = 0;
    input->kdtreeRoot = NULL;
    input->r = -1;
    input->silent = 0;
    input->display = 0;
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
    dsname = basename(strdup(input->filename));
    input->ds = load_data(fname, &input->n, &input->k);

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
    input-> Point= (float*) malloc(sizeof(float)* (input->h>0? input->h : input->k) );
    input-> Qoint= (float*) malloc(sizeof(float)* (input->h>0? input->h : input->k) );
    input-> vetMediaDS = (float*) malloc(sizeof(float) * input-> k);
    input-> dsTras = (float*) malloc( sizeof(float) * input->n * input-> k );
    input-> qsRidotto = (float*) malloc( sizeof(float) * input->nq * input-> h );
    input-> QA = (int*) malloc(sizeof(int)* input->nq*input->n);
    input-> pQA = (float**) malloc(sizeof(float*)* input->nq*input->n);

    for(i=0; i<input->nq*input->n; i++){
        input->pQA[i] = (float*) malloc(sizeof(float));
    }

    input-> vetTmp = (int*) malloc(sizeof(int)*input->n);

    //
    // Calcolo PCA
    //
    
    if(input->h > 0){
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
        sprintf(fname, "%s.U", dsname);
        sprintf(fname, "%s.V", dsname);
    }else
        time = -1;
       
    if (!input->silent)
        printf("\nPCA time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    //
    // Costruzione K-d-Tree
    //
    
    if(input->kdtree_enabled){
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
    
    int z;
  
    float* dataS; float* queryS;
    int nCol;
    if(input->h>0){
         nCol = input->h;
         dataS = (float*) malloc(sizeof(float)* input->n * nCol );
         dataS = input->U;
         queryS = (float*) malloc(sizeof(float)* input->nq * nCol);
         queryS = input->qsRidotto;
    }else{
         nCol = input->k;
         dataS = (float*) malloc(sizeof(float)* input->n * nCol );
         dataS = input->ds;
         queryS = (float*) malloc(sizeof(float)* input->nq * nCol);
         queryS = input->qs;
    }

    int* queryAnswer = (int*)malloc(sizeof(int)*input->nQA);
    for(i=0; i<input->nQA; i++){
        queryAnswer[i]=-1;
    }
    int booleanUguali;
    for(i=0; i<input->nQA; i++){
        for(j=0; j<input->n; j++){
            booleanUguali=0;
            if( input->pQA[i][0] == dataS[j*nCol+0] ){
                booleanUguali=1;
                for(z=1; z<nCol && booleanUguali==1; z++){
                    if( input->pQA[i][z] != dataS[j*nCol+z] ){
                        booleanUguali=0;
                    }
                }
                if(booleanUguali==1){
                    queryAnswer[i] = j;
                }
            }
        }
    }

    int boolean = 0;
    if(input->r >= 0){
        if(!input->silent && input->display) {
            //NB: il codice non assume che QA sia ordinata per query, in caso lo sia ottimizzare il codice
            printf("\nQuery Answer:\n");
            for(i = 0; i < input->nq; i++){
                
            j=0;
        boolean = 0;
        while(j<input->nQA && boolean==0 ){
            if(i == input->QA[j] ){
                boolean=1;
            }
            j+=1;
        }
        if(boolean==0){
            printf("query %d: [ ]\n", i);
        }else{
            printf("query %d: [ ", i);
            j=0;
            while( j<input->nQA){
                if( i == input->QA[j] ){
                    printf("%d ", queryAnswer[j]);
                }
                j+=1;
            }
            printf("]\n");        
        }


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
