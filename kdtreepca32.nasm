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

		; elaborazione
		
		; esempio: stampa input->n e di input->k
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
		; [EAX] contiene l'indirizzo della stringa con il nome del file
		; [EAX+4] contiene l'indirizzo di partenza del data set
		; [EAX+8] contiene l'indirizzo di partenza del query set
;		printi dword[EAX+12]	; a 12 byte dall'inizio della struct si trova n
;		printi dword[EAX+16]	; a 4 byte da n si trova k
		

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
	jne uscita_prodMatrVett_ass

	mov edx,0
	mov eax, [ebp+16]
	div dword[quattro]
	mov edi, eax ; ho ottenuto il nRighe/4 in EDI

	xor esi,esi
	ciclo_esterno_1:
		cmp esi,edi
		jge fine_ciclo_esterno_1
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
		ciclo_interno_1:
			cmp ecx,ebx
			jge fine_ciclo_interno_1

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
			div dword[quattro]
			mov esi,eax
			; riporto esi al valore del ciclo esterno 1

			jmp ciclo_interno_1

		fine_ciclo_interno_1:
		
		xor edi,edi
		; ci manca il ciclo interno 2 in cui cicliamo sulle singole colonne dopo aver
		; controllato che non abbiamo raggiunto il nColonne giusto (ecx*4 == nColonne)
		mov edx,0
		mov eax,4
		mul ecx
		mov ecx,eax
		mov ebx, [ebp+20]

		cmp ecx,ebx
		je concluso_4_righe

		ciclo_interno_2:
			cmp ecx,ebx
			jge concluso_4_righe

		mov edx, ecx
		push edx
		push intero
		call printf
		pop edx
		pop edx
		pushad
		call getchar
		popad

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

		movss [tmp],xmm0
		mov edx,[tmp]
		push edx
		push esad
		call printf
		pop edx
		pop edx
		pushad
		call getchar
		popad


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
			_00:movss xmm3, [edx+eax]

			; nei primi 32 bit di xmm0,xmm1,xmm2,xmm3 ho i valori delle colonne
			
			_10:mov edx,0
			_12:mov eax,4
			_13:mul ecx
			_14:mov edx,[ebp+12]
			_15:movss xmm4, [edx+eax]
			; ho in xmm4 l'elemento del vettore

			_20:mulss xmm0,xmm4
			_22:mulss xmm1,xmm4
			_21:mulss xmm2,xmm4
			_23:mulss xmm3,xmm4

			_30:shufps xmm0,xmm1,0
			_31:shufps xmm2,xmm3,0
			_32:shufps xmm0,xmm2,136 ; come prima

			_44:addps xmm7,xmm0

			_55:sub esi,3
			_51:mov edx,0
			_52:mov eax,esi
			_53:div dword[quattro]
			_54:mov esi,eax

			add ecx,1
			jmp ciclo_interno_2

		concluso_4_righe:
		pop edi

		mov ebx,[ebp+28] ; prendo il puntatore del vettore risultato
		mov edx,0
		mov eax,16
		mul edi
		movups [ebx+eax],xmm7

	fine_ciclo_esterno_1:
	mov edx,0
	mov eax,esi
	mul dword[quattro]
	mov esi,eax
	mov edi,[ebp+16]
	cmp esi,edi ; controllo se ho fatto tutte le righe con 4 alla volta
	je uscita_prodMatrVett_ass

	ciclo_esterno_2:	; svolgiamo una riga alla volta
		cmp esi,edi
		jge uscita_prodMatrVett_ass







	uscita_prodMatrVett_ass:
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

	mov ecx, [ebp+16] ; prendo il puntatore alla norma

	mov edx,0
	mov eax,[ebp+12]
	div dword[quattro]
	mov edi,eax

	xor esi, esi
	ciclo_quoz:
		cmp esi, edi
		je continua
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
		jmp ciclo_quoz

	continua:
		mov edx,0
		mov eax,esi
		mul dword[quattro]
		mov esi,eax
		mov edi, [ebp+12]

		cmp esi,edi
		je fine_calcolaNorma_ass

	ciclo_resto:
		cmp esi,edi
		je fine_calcolaNorma_ass
		mov ebx,[ebp+8]
		movss xmm0,[ebx+4*esi]
		mulss xmm0,xmm0
		xorps xmm1,xmm1
		movss xmm1, [ecx]
		addss xmm1, xmm0
		movss [ecx],xmm1
		add esi,1
		jmp ciclo_resto


	fine_calcolaNorma_ass:
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

