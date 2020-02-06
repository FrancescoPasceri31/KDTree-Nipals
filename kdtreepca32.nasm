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
stringa: db 0,'Stampa',0
quattro: dd 4
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4
tmp: resd 1

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
		push 		eax
		push		ecx
		push		edx
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		mov	eax, [ebp+8] ;v
		mov ebx, [ebp+12] ;dim
		;mov ecx, [ebp+16]

		xor esi, esi
		xor edx, edx
ciclo_1:
;;		div ebx,4
		mov esi,[ebx]
		xor edi,edi
		cmp edi, esi
		jge fine_ciclo
;;		mov xmm0,[eax+4*edi]
		
		;problema?
		add edi,1
fine_ciclo:
;;		mul esi,4				; indice da cui partire (esi=16)
		mov edi, [ebx]
		sub edi,esi				; 19 -16 =3 
ciclo_2:
		cmp edi,0
		je fine_tutto_pt1
		mov ecx,[eax+esi]		;elementi
;;		mul ecx, ecx
		add edx, ecx			;somma += elementoi
		add esi,1
		sub edi,1
		jmp ciclo_2

		

fine_tutto_pt1:
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop edx
		pop ecx
		pop eax
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante


tess:
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

push stringa
call printf
pop ecx


	mov eax, [ebp+8] ; prendo punt array
	mov edi, [ebp+12] ; prendo length in edi
	; ho la norma in [ebp+16]

	; per calcolare la norma devo fare loop unrolling quoziente su un registro xmm per length/4 volte
	; e successivamente per la differenza devo eseguire un loop resto

	movss xmm0, [ebp+12]
	shufps xmm0,xmm0,0
	movss xmm1, [quattro]
	shufps xmm1,xmm1,0
	divps xmm0,xmm1	; metto in xmm0 la divisione length/4 
	movss [tmp], xmm0 ; salvo in edi la length del vettore / 4
	mov edi, [tmp]
	xorps xmm0,xmm0
	xorps xmm1,xmm1

	xor esi, esi
	ciclo_quoz:
		cmp esi, edi
		je ciclo_resto
		movups xmm0, [eax+4*esi]
		movss [tmp],xmm0
		mulps xmm0,xmm0
		haddps xmm0,xmm0
		haddps xmm0,xmm0
		xorps xmm1,xmm1
		movss xmm1, [ebp+16]
		addss xmm1,xmm0
		movss [ebp+16],xmm1
		add esi,1
		jmp ciclo_quoz


	xor esi,esi
	mov [tmp],edi
	movss xmm0,[tmp]
	mulss xmm0,[quattro]	; moltiplico il numero di cicli precedenti * 4 in maniera tale da trovare la posizione dove riempire
	movss [tmp],xmm0	; cella da riempire
	mov esi, [tmp]
	xor edi,edi
	mov edi, [ebp+12]	; carico in edi la length del vettore ed eseguo length-cella dove sono
	ciclo_resto:
		cmp esi,edi
		je fine_tess
		mov edx, [eax+4*esi]
		mov [tmp],edx
		movss xmm0,[tmp]
		mulss xmm0,xmm0
		xorps xmm1,xmm1
		movss xmm1, [ebp+16]
		addss xmm1, xmm0
		movss [ebp+16],xmm1
		add esi,1
		jmp ciclo_resto

	fine_tess:
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante

