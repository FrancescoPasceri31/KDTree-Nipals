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

global tess


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

