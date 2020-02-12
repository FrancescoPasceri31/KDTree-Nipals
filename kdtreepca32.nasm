; ---------------------------------------------------------
; PQNN con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 pqnn32.nasm 
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
esad: db 'valore: %x',10,0
intero: db 'valore: %d',10,0
duevalori: db '%d %d',10,0
trevalori: db '%d %d %x',10,0
stampaXMM: db '%x %x %x %x',10,0
quattro: dd 4
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4
tmp: resd 1
count: resd 1
xmmTMP: resd 4

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block
extern printf
extern getchar
extern scanf

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global prova

global calcolaNorma_ass
global prodMatrVett_ass
global prodScalare_ass
global divisioneVettoreScalare_ass
global sottrazioneMatrici_ass
global distanzaEuclidea_ass


input		equ		8

prova:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante



prodMatr_ass:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------
	
	
	
	
	
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante


distanzaEuclidea_ass:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; float* P +8 ,float* Q +12 , int dimen +16, float* dist +20
		xorps xmm0,xmm0
		mov ecx,[ebp+20]
		movss [ecx],xmm0
		
		mov edx,0
		mov eax,[ebp+16]
		div dword[quattro]
		mov edi, eax
		
		xorps xmm7,xmm7

		xor esi,esi
		ciclo_quoziente_DE:
			cmp esi,edi
			jge fine_quoziente_DE

			mov edx,0
			mov eax,16
			mul esi
			mov edx,[ebp+8]
			mov ebx, [ebp+12]
			movups xmm0, [edx+eax]
			movups xmm1, [ebx+eax]
			; ho caricato i 4 elementi dei due vettori

			subps xmm0,xmm1
			mulps xmm0,xmm0

			haddps xmm0,xmm0
			haddps xmm0,xmm0

			addss xmm7,xmm0

			add esi,1
			jmp ciclo_quoziente_DE

		fine_quoziente_DE:
		mov edx,0
		mov eax,4
		mul esi
		mov esi,eax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

		mov edi, [ebp+16]
		cmp esi,edi
		je fine_ciclo_resto_DE

		ciclo_resto_DE:
			cmp esi,edi
			jge fine_ciclo_resto_DE

			mov edx,0
			mov eax,4
			mul esi
			mov ebx,[ebp+8]
			mov ecx,[ebp+12]
			movss xmm0, [ebx+eax]
			movss xmm1, [ecx+eax]
			; ho caricato i due elementi

			subss xmm0,xmm1
			mulss xmm0,xmm0

			addss xmm7,xmm0

			add esi,1
			jmp ciclo_resto_DE

		fine_ciclo_resto_DE:
		
		mov eax,[ebp+20] ; ho preso il valore del puntatore del risultato
		sqrtss xmm7,xmm7
		movss [eax],xmm7

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante

sottrazioneMatrici_ass:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------

	; float* m1 +8, float* m2 +12, int nRighe1 +16, int nColonne1 +20, int nRighe2 +24, int nColonne2 +28
	; float* risultato +32
	
	mov eax, [ebp+16]
	cmp eax,[ebp+24]
	jne uscita_sottrazioneMatrice_ass_SM
	mov eax, [ebp+20]
	cmp eax,[ebp+28]
	jne uscita_sottrazioneMatrice_ass_SM

	mov edx,0
	mov eax, [ebp+16]
	div dword[quattro]
	mov edi, eax ; ho ottenuto il nRighe/4 in EDI

	xor esi,esi
	ciclo_esterno_1_SM:
		cmp esi,edi
		jge fine_ciclo_esterno_1_SM
		; altrimenti prendo 4 righe alla volta

		mov edx,0
		mov eax, [ebp+20]
		div dword[quattro]	
		mov ebx,eax	; in ebx ho nColonne/4

		push edi

		mov edx,0
		mov eax,4
		mul esi
		mov esi,eax
		; inserisco in ESI = ESI*4 per trovare il primo indice riga delle 4 da analizzare

		xor ecx,ecx
		ciclo_interno_1_SM:
			cmp ecx,ebx
			jge fine_ciclo_interno_1_SM

			; in xmm0-1-2-3 mantengo la prima matrice
			; in xmm4-5-6-7 mantengo la seconda matrice

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+8]
			movups xmm0, [edx+eax]	; ho caricato la prima riga
			mov edx, [ebp+12]
			movups xmm4, [edx+eax]

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm1, [edx+eax]
			mov edx, [ebp+12]
			movups xmm5, [edx+eax]
	
			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm2, [edx+eax]
			mov edx, [ebp+12]
			movups xmm6, [edx+eax]

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm3, [edx+eax]
			mov edx, [ebp+12]
			movups xmm7, [edx+eax]
			; ho caricato le 4 righe in xmm0,xmm1,xmm2,xmm3

			subps xmm0,xmm4
			subps xmm1,xmm5
			subps xmm2,xmm6
			subps xmm3,xmm7
			; ho sottratto gli elementi e devo inserirli in memoria

			sub esi,3
			mov edx,0
			mov eax,esi
			mov esi,eax			

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+32]
			movups [edx+eax], xmm0	; ho caricato la prima riga

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+32]
			movups [edx+eax],xmm1
	
			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+32]
			movups [edx+eax],xmm2

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+32]
			movups [edx+eax],xmm3
			; ho caricato le 4 righe da xmm0,xmm1,xmm2,xmm3
			; in memoria

			add ecx,1

			sub esi,3
			mov edx,0
			mov eax,esi
			;div dword[quattro]
			mov esi,eax
			; riporto esi al valore del ciclo esterno 1
			jmp ciclo_interno_1_SM

		fine_ciclo_interno_1_SM:
		
		xor edi,edi
		; ci manca il ciclo interno 2 in cui cicliamo sulle singole colonne dopo aver
		; controllato che non abbiamo raggiunto il nColonne giusto (ecx*4 == nColonne)
		mov edx,0
		mov eax,4
		mul ecx
		mov ecx,eax
		mov ebx, [ebp+20]

		cmp ecx,ebx
		je concluso_4_righe_SM

		ciclo_interno_2_SM:
			mov ebx, [ebp+20]
			cmp ecx,ebx
			jge concluso_4_righe_SM

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm0, [edx+eax]
			mov edx, [ebp+12]
			movss xmm4, [edx+eax]
 
			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm1, [edx+eax]
			mov edx, [ebp+12]
			movss xmm5, [edx+eax]

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm2, [edx+eax]
			mov edx, [ebp+12]
			movss xmm6, [edx+eax]

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm3, [edx+eax]
			mov edx, [ebp+12]
			movss xmm7, [edx+eax]
			; nei primi 32 bit di xmm0,xmm1,xmm2,xmm3 ho i valori delle colonne
			
			subss xmm0,xmm4
			subss xmm1,xmm4
			subss xmm2,xmm4
			subss xmm3,xmm4
			; ho diviso tutti gli elementi per lo scalare

			sub esi,3
			mov edx,0
			mov eax,esi
			mov esi,eax

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+32]
			movss [edx+eax],xmm0
 
			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+32]
			movss [edx+eax],xmm1

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+32]
			movss [edx+eax],xmm2

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+32]
			movss [edx+eax],xmm3

			sub esi,3
			mov edx,0
			mov eax,esi
			mov esi,eax

			add ecx,1
			jmp ciclo_interno_2_SM

		concluso_4_righe_SM:
		pop edi

		mov edx,0
		mov eax,esi
		div dword[quattro]
		mov esi,eax

		add esi,1
		jmp ciclo_esterno_1_SM

	fine_ciclo_esterno_1_SM:
	mov edx,0
	mov eax,esi
	mul dword[quattro]
	mov esi,eax
	
	mov edi,[ebp+16]
	cmp esi,edi ; controllo se ho fatto tutte le righe con 4 alla volta
	je uscita_sottrazioneMatrice_ass_SM

	ciclo_esterno_2_SM:	; svolgiamo una riga alla volta
		cmp esi,edi
		jge uscita_sottrazioneMatrice_ass_SM
		
		push edi

		; in esi ho l'indice riga
		mov edx,0
		mov eax, [ebp+20]
		div dword[quattro] ; eax = nColonne/4
		mov ebx, eax

		xorps xmm7,xmm7

		xor ecx,ecx
		ciclo_quoziente_2_SM:
			cmp ecx,ebx
			jge fine_ciclo_quoziente_2_SM

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+8]
			movups xmm0, [edx+eax]	; ho caricato la prima riga
			mov edx, [ebp+12]
			movups xmm4, [edx+eax]

			subps xmm0,xmm4

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+32]
			movups [edx+eax],xmm0	; ho caricato la prima riga

			add ecx,1
			jmp ciclo_quoziente_2_SM

		fine_ciclo_quoziente_2_SM:
		; ho terminato le colonne prese a 4 ora conto le singole

		mov edx,0
		mov eax,4
		mul ecx
		mov ecx,eax
		mov ebx, [ebp+20]

		cmp ecx,ebx
		jge fine_ciclo_resto_2_SM

		ciclo_resto_2_SM:
			mov ebx,[ebp+20]
			cmp ecx,ebx
			jge fine_ciclo_resto_2_SM

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna	
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm0, [edx+eax]
			mov edx, [ebp+12]
			movss xmm4, [edx+eax]
			; nei primi 32 di xmm0 ho inserito l'elemento della matrice

			subss xmm0,xmm4

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna	
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+32]
			movss [edx+eax],xmm0

			add ecx,1
			jmp ciclo_resto_2_SM

		fine_ciclo_resto_2_SM:
		
		pop edi

		add esi,1
		jmp ciclo_esterno_2_SM

	uscita_sottrazioneMatrice_ass_SM:
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante




divisioneVettoreScalare_ass:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------

	; v, scalare, nColonne, rs
	; v = +8
	; sca = +12
	; nCol = +16
	; rs = +20 -> è un puntatore

	mov edx,0
	mov eax,[ebp+16]
	div dword[quattro]
	mov edi,eax ; in edi ho il numero di volte che posso ciclare con xmm0 = nColonne/4

	xor esi,esi
	ciclo_quoziente_DVS:
		cmp esi,edi
		jge fine_ciclo_quoziente_DVS

		;xorps xmm0,xmm0
		mov edx,0
		mov eax,16
		mul esi
		mov edx,[ebp+8]
		movups xmm0, [edx+eax]
		; ho caricato i 4 elementi del vettore

		movss xmm1, [ebp+12]
		shufps xmm1,xmm1,0

		divps xmm0,xmm1

		mov edx,[ebp+20]
		movups [edx+eax],xmm0

		add esi, 1
		jmp ciclo_quoziente_DVS

	fine_ciclo_quoziente_DVS:
	mov edx,0
	mov eax,4
	mul esi
	mov esi,eax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

	mov edi, [ebp+16]
	cmp esi,edi
	je fine_ciclo_resto_DVS

	ciclo_resto_DVS:
		cmp esi,edi
		jge fine_ciclo_resto_DVS

		mov edx,0
		mov eax,4
		mul esi
		mov ebx,[ebp+8]
		movss xmm0, [ebx+eax]
		; ho caricato l'elemento del vettore

		movss xmm1, [ebp+12]
		; ho caricato lo scalare

		divss xmm0,xmm1

		mov edx,[ebp+20]
		movss [edx+eax],xmm0

		add esi,1
		jmp ciclo_resto_DVS

	fine_ciclo_resto_DVS:

	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante



prodScalare_ass:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------


		; (float* v1,int dim1,float* v2,int dim2,float* rs)
		; v1 -> ebp+8
		; dim1 -> ebp+12
		; v2 -> ebp+16
		; dim2 -> ebp+20
		; rs -> ebp+24 è un puntatore per cui il valore va in [[ebp+24]]

		mov eax, [ebp+12]
		cmp eax,[ebp+20]
		jne uscita_prodottoScalare_ass_PS

		mov edx,0
		mov eax,[ebp+12]
		div dword[quattro]
		mov edi,eax ; in edi ho il numero di volte che posso ciclare con xmm0 = nColonne/4

		xorps xmm7,xmm7

		xor esi,esi
		ciclo_quoziente_PS:
			cmp esi,edi
			jge fine_ciclo_quoziente_PS

			mov edx,0
			mov eax,16
			mul esi
			mov edx,[ebp+8]
			mov ebx, [ebp+16]
			movups xmm0, [edx+eax]
			movups xmm1, [ebx+eax]
			; ho caricato i 4 elementi dei due vettori

			mulps xmm0,xmm1

			haddps xmm0,xmm0
			haddps xmm0,xmm0

			addss xmm7,xmm0

			add esi, 1
			jmp ciclo_quoziente_PS

		fine_ciclo_quoziente_PS:
		mov edx,0
		mov eax,4
		mul esi
		mov esi,eax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

		mov edi, [ebp+12]
		cmp esi,edi
		je fine_ciclo_resto_PS

		ciclo_resto_PS:
			cmp esi,edi
			jge fine_ciclo_resto_PS

			mov edx,0
			mov eax,4
			mul esi
			mov ebx,[ebp+8]
			mov ecx,[ebp+16]
			movss xmm0, [ebx+eax]
			movss xmm1, [ecx+eax]
			; ho caricato i due elementi

			mulss xmm0,xmm1

			addss xmm7,xmm0

			add esi,1
			jmp ciclo_resto_PS

		fine_ciclo_resto_PS:
		
		mov eax,[ebp+24] ; ho preso il valore del puntatore del risultato
		movss [eax],xmm7

		uscita_prodottoScalare_ass_PS:
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante






prodMatrVett_ass:

	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------

	; elaborazione

	; float* m = ebp+8 , 
	; float* v = ebp+12, 
	; int nRighe = ebp+16, 
	; int nColonne = ebp+20,
	; int lengthVettore = ebp+24,
	; float* risultato = ebp+28

	mov eax, [ebp+20]
	mov ebx, [ebp+24]
	cmp eax,ebx
	jne uscita_prodMatrVett_ass_PMV

	mov edx,0
	mov eax, [ebp+16]
	div dword[quattro]
	mov edi, eax ; ho ottenuto il nRighe/4 in EDI

	xor esi,esi
	ciclo_esterno_1_PMV:
		cmp esi,edi
		jge fine_ciclo_esterno_1_PMV
		; altrimenti prendo 4 righe alla volta

		mov edx,0
		mov eax, [ebp+20]
		div dword[quattro]	
		mov ebx,eax	; in ebx ho nColonne/4

		xorps xmm7,xmm7 ; sum in xmm7

		push edi

		mov edx,0
		mov eax,4
		mul esi
		mov esi,eax
		; inserisco in ESI = ESI*4 per trovare il primo indice riga delle 4 da analizzare

		xor ecx,ecx
		ciclo_interno_1_PMV:
			cmp ecx,ebx
			jge fine_ciclo_interno_1_PMV

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+8]
			movups xmm0, [edx+eax]	; ho caricato la prima riga

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm1, [edx+eax]
	
			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm2, [edx+eax]

			add esi,1
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			mov edx, [ebp+8]
			movups xmm3, [edx+eax]
			; ho caricato le 4 righe in xmm0,xmm1,xmm2,xmm3

			mov edx,0
			mov eax, 16
			mul ecx
			; in eax ho ecx*16 = indice del vettore
			mov edx, [ebp+12]
			movups xmm4, [edx+eax]
			; ho caricato il vettore in xmm4

			mulps xmm0,xmm4
			mulps xmm1,xmm4
			mulps xmm2,xmm4
			mulps xmm3,xmm4
			
			haddps xmm0,xmm0
			haddps xmm1,xmm1
			haddps xmm2,xmm2
			haddps xmm3,xmm3

			haddps xmm0,xmm0
			haddps xmm1,xmm1
			haddps xmm2,xmm2
			haddps xmm3,xmm3
			; ho la somma dei 4 elementi della matrice moltiplcati per i rispettivi del vettore nei primi 32 bit

			shufps xmm0,xmm1,0
			shufps xmm2,xmm3,0	; [00-00-00-00]
			shufps xmm0,xmm2,136 ;[10-00-10-00]
			; ho ottenuto il vettore con le 4 somme [a,b,c,d]

			addps xmm7,xmm0
			
			add ecx,1

			sub esi,3
			mov edx,0
			mov eax,esi
			;div dword[quattro]
			mov esi,eax
			; riporto esi al valore del ciclo esterno 1

			jmp ciclo_interno_1_PMV

		fine_ciclo_interno_1_PMV:
		
		xor edi,edi
		; ci manca il ciclo interno 2 in cui cicliamo sulle singole colonne dopo aver
		; controllato che non abbiamo raggiunto il nColonne giusto (ecx*4 == nColonne)
		mov edx,0
		mov eax,4
		mul ecx
		mov ecx,eax
		mov ebx, [ebp+20]

		cmp ecx,ebx
		je concluso_4_righe_PMV

		ciclo_interno_2_PMV:
			mov ebx, [ebp+20]
			cmp ecx,ebx
			jge concluso_4_righe_PMV

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm0, [edx+eax]
 
			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm1, [edx+eax]

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm2, [edx+eax]

			xor eax,eax
			add esi,1
			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna		
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm3, [edx+eax]
			; nei primi 32 bit di xmm0,xmm1,xmm2,xmm3 ho i valori delle colonne
			
			mov edx,0
			mov eax,4
			mul ecx
			mov edx,[ebp+12]
			movss xmm4, [edx+eax]
			; ho in xmm4 l'elemento del vettore

			mulss xmm0,xmm4
			mulss xmm1,xmm4
			mulss xmm2,xmm4
			mulss xmm3,xmm4

			shufps xmm0,xmm1,0
			shufps xmm2,xmm3,0
			shufps xmm0,xmm2,136 ; come prima

			addps xmm7,xmm0

			sub esi,3
			mov edx,0
			mov eax,esi
			;_53:div dword[quattro]
			mov esi,eax

			add ecx,1
			jmp ciclo_interno_2_PMV

		concluso_4_righe_PMV:
		pop edi

		mov edx,0
		mov eax,esi
		div dword[quattro]
		mov esi,eax

		mov ebx,[ebp+28] ; prendo il puntatore del vettore risultato
		mov edx,0
		mov eax,16
		mul esi

		movups [ebx+eax],xmm7

		add esi,1
		jmp ciclo_esterno_1_PMV

	fine_ciclo_esterno_1_PMV:
	mov edx,0
	mov eax,esi
	mul dword[quattro]
	mov esi,eax
	mov edi,[ebp+16]
	cmp esi,edi ; controllo se ho fatto tutte le righe con 4 alla volta
	je uscita_prodMatrVett_ass_PMV

	ciclo_esterno_2_PMV:	; svolgiamo una riga alla volta
		cmp esi,edi
		jge uscita_prodMatrVett_ass_PMV
		
		push edi

		; in esi ho l'indice riga
		mov edx,0
		mov eax, [ebp+20]
		div dword[quattro] ; eax = nColone/4
		mov ebx, eax

		xorps xmm7,xmm7

		xor ecx,ecx
		ciclo_quoziente_2_PMV:
			cmp ecx,ebx
			jge fine_ciclo_quoziente_2_PMV

			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			xor edi,edi
			mov edi, eax
			mov edx,0
			mov eax, 16
			mul ecx
			add eax,edi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov edx, [ebp+8]
			movups xmm0, [edx+eax]	; ho caricato la prima riga

			mov edx,0
			mov eax, 16
			mul ecx
			; in eax ho ecx*16 = indice del vettore
			mov edx, [ebp+12]
			movups xmm4, [edx+eax]
			; ho caricato il vettore in xmm4

			mulps xmm0,xmm4

			haddps xmm0,xmm0
			haddps xmm0,xmm0
			; nei primi 32 bit ho la somma di 4 elementi moltiplcati per i 4 del vettore

			addss xmm7,xmm0
			; ho aggiunto la somma temporanea nel registro

			add ecx,1
			jmp ciclo_quoziente_2_PMV

		fine_ciclo_quoziente_2_PMV:
		; ho terminato le colonne prese a 4 ora conto le singole

		mov edx,0
		mov eax,4
		mul ecx
		mov ecx,eax
		mov ebx, [ebp+20]

		cmp ecx,ebx
		jge fine_ciclo_resto_2_PMV

		ciclo_resto_2_PMV:
			mov ebx,[ebp+20]
			cmp ecx,ebx
			jge fine_ciclo_resto_2_PMV

			mov edx,0
			mov eax,4
			mul ecx
			mov ebx, eax ; in ebx ho ecx*4 che è un indice di colonna	
			mov edx,0
			mov eax, [ebp+20]
			mul dword[quattro]
			mul esi
			add eax,ebx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov edx, [ebp+8]
			movss xmm0, [edx+eax]
			; nei primi 32 di xmm0 ho inserito l'elemento della matrice

			mov edx,0
			mov eax,4
			mul ecx
			mov edx,[ebp+12]
			movss xmm4, [edx+eax]
			; ho in xmm4 l'elemento del vettore
			
			mulss xmm0,xmm4

			addss xmm7,xmm0

			add ecx,1
			jmp ciclo_resto_2_PMV

		fine_ciclo_resto_2_PMV:
		
		; prendere i primi 32 bit di xmm7 ed aggiungerli in posizione del vettore ebp+28
		mov edx,0
		mov eax,4
		mul esi
		
		mov ecx,[ebp+28]
		movss [ecx+eax],xmm7

		pop edi

		add esi,1
		jmp ciclo_esterno_2_PMV

	uscita_prodMatrVett_ass_PMV:
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante






calcolaNorma_ass:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi
	; ------------------------------------------------------------
	; legge i parametri dal Record di Attivazione corrente
	; ------------------------------------------------------------

	; 8 puntatore, 12 la length, 16 puntatore alla norma

	; per calcolare la norma devo fare loop unrolling quoziente su un registro xmm per length/4 volte
	; e successivamente per la differenza devo eseguire un loop resto

	xorps xmm0,xmm0
	mov ecx, [ebp+16] ; prendo il puntatore alla norma
	movss [ecx],xmm0

	mov edx,0
	mov eax,[ebp+12]
	div dword[quattro]
	mov edi,eax

	xor esi, esi
	ciclo_quoz_CN:
		cmp esi, edi
		je continua_CN
		mov ebx,[ebp+8] ; prendo il puntatore all'array
		mov edx,0
		mov eax, 16 ; moltiplico esi per 16 in maniera tale da saltare 4 elementi
		mul esi
		movups xmm0, [ebx+eax] ; prendo 4 elementi
		mulps xmm0,xmm0
		haddps xmm0,xmm0 ; sommo pos 0 e 1 di xmm
		haddps xmm0,xmm0 ; stessa cosa
		xorps xmm1,xmm1 ; azzero xmm1
		movss xmm1, [ecx] ; inserisco in xmm1 la norma
		addss xmm1,xmm0 ; sommo norma a quello calcolato prima
		movss [ecx],xmm1 ; la rimetto in memoria
		add esi,1
		jmp ciclo_quoz_CN

	continua_CN:
		mov edx,0
		mov eax,esi
		mul dword[quattro]
		mov esi,eax
		mov edi, [ebp+12]

		cmp esi,edi
		je fine_calcolaNorma_ass_CN

	ciclo_resto_CN:
		cmp esi,edi
		je fine_calcolaNorma_ass_CN
		mov ebx,[ebp+8]
		movss xmm0,[ebx+4*esi]
		mulss xmm0,xmm0
		xorps xmm1,xmm1
		movss xmm1, [ecx]
		addss xmm1, xmm0
		movss [ecx],xmm1
		add esi,1
		jmp ciclo_resto_CN


	fine_calcolaNorma_ass_CN:
	movss xmm0, [ecx]
	sqrtss xmm0,xmm0
	movss [ecx],xmm0

	;mov edx,[ecx]
	;push edx
	;push esad
	;call printf
	;pop edx
	;pop edx
	;pushad
	;call getchar
	;popad

	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante






