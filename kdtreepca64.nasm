; ---------------------------------------------------------
; PageRank con istruzioni AVX a 64 bit
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
;     nasm -f elf64 pagerank64.nasm 
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
esad: db 'valore: %x',10,0
intero: db 'valore: %d',10,0
duevalori: db '%d %d',10,0
trevalori: db '%d %d %x',10,0
stampaXMM: db '%x %x %x %x',10,0
otto: dd 8
quattro: dd 4
;
;align 32
;vec1:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 32
;vec2:		resq	4
tmp: resd 2
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

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global prova
global calcolaNorma_ass_64
global prodScalare_ass_64
global distanzaEuclidea_ass_64
global divisioneVettoreScalare_ass_64

msg	db 'n:',0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		
		; rdi = indirizzo della struct input
		; esempio: stampa input->n e di input->k
		; rdi contiente l'indirizzo della struttura contenente i parametri
		; [rdi] contiene l'indirizzo della stringa con il nome del file
		; [rdi+8] contiene l'indirizzo di partenza del data set
		; [rdi+16] contiene l'indirizzo di partenza del query set
		; rdi -> puntatore ad elemento
		; rsi -> puntatore a scalare
		; rdx -> puntatore a risultato
		vbroadcastss ymm0, [rdi]
		vbroadcastss ymm1, [rsi]
		vpermilps ymm0,ymm0,0
		vmulss xmm0,xmm0,xmm1
		vmovss [rdx],xmm0
		;movsx rax, dword[rdi+28]		; a 4 byte da n si trova k
		;printi rax

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante

divisioneVettoreScalare_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------
	
	; v = +8
	; sca = +12 
	; nCol = +16
	; rs = +20 -> è un puntatore

	mov r8, rdi
	mov r9, rsi
	mov r10, rdx
	mov r11, rcx 

	mov rdx,0
	mov rax,r10
	div dword[otto]
	mov rdi,rax ; in rdi ho il numero di volte che posso ciclare con ymm0 = nColonne/8

	xor rsi,rsi
	ciclo_quoziente_DVS:
		cmp rsi,rdi
		jge fine_ciclo_quoziente_DVS

		mov rbx,r8 ; prendo il puntatore all'array
		mov rdx,0
		mov rax, 32 ; moltiplico esi per 32 in maniera tale da saltare 8 elementi
		mul rsi
		vmovups ymm0, [rbx+rax] ; prendo 8 elementi

		vbroadcastss ymm1, [r9]

		vdivps ymm0,ymm0,ymm1

		;call printf

		vmovups [r11+rax],ymm0

		add rsi, 1
		jmp ciclo_quoziente_DVS

	fine_ciclo_quoziente_DVS:
	mov rdx,0
	mov rax,8
	mul rsi
	mov rsi,rax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

	ciclo_resto_DVS:
		cmp rsi,r10
		jge fine_ciclo_resto_DVS

		mov rdx,0
		mov rax,4
		mul rsi
		mov rbx,r8
		vmovss xmm2, [rbx+rax]
		; ho caricato l'elemento del vettore

		vmovss xmm1, [r9]
		; ho caricato lo scalare

		vdivss xmm2,xmm1

		vmovss [r11+rax],xmm2

		add rsi,1
		jmp ciclo_resto_DVS

	fine_ciclo_resto_DVS:

	; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante







distanzaEuclidea_ass_64:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------

		; float* P +8 ,float* Q +12 , int dimen +16, float* dist +20
		
		mov r8, rdi
		mov r9, rsi
		mov r10, rdx
		mov r11, rcx

		vxorps xmm0,xmm0,xmm0
		vmovss [r11],xmm0

		mov rdx,0
		mov rax,r10
		div dword[otto]
		mov rdi, rax
		
		vxorps ymm7,ymm7,ymm7

		xor rsi,rsi
		ciclo_quoziente_DE:
			cmp rsi,rdi
			jge fine_quoziente_DE

			mov rdx,0
			mov rax,32
			mul rsi
			mov rdx,r8
			mov rbx, r9
			vmovups ymm0, [rdx+rax]
			vmovups ymm1, [rbx+rax]
			; ho caricato i 4 elementi dei due vettori

			vsubps ymm0,ymm0,ymm1
			vmulps ymm0,ymm0,ymm0

			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm0,ymm0,ymm0

			vaddss xmm7,xmm0

			add rsi,1
			jmp ciclo_quoziente_DE

		fine_quoziente_DE:
		mov rdx,0
		mov rax,8
		mul rsi
		mov rsi,rax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

		ciclo_resto_DE:
			cmp rsi,r10
			jge fine_ciclo_resto_DE

			mov rdx,0
			mov rax,4
			mul rsi
			mov rbx,r8
			mov rcx,r9
			vmovss xmm0, [ebx+eax]
			vmovss xmm1, [ecx+eax]
			; ho caricato i due elementi

			vsubss xmm0,xmm1
			vmulss xmm0,xmm0

			vaddss xmm7,xmm0

			add rsi,1
			jmp ciclo_resto_DE

		fine_ciclo_resto_DE:
		
		vsqrtss xmm7,xmm7,xmm7
		movss [r11],xmm7

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante



prodScalare_ass_64:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------

		; (float* v1,int dim1,float* v2,int dim2,float* rs)

		mov rax,r8
		mov r8,rdi	; v1
		mov r9, rsi	; dim1
		mov r10, rdx	; v2
		mov r11, rcx	; dim2
		mov r12,rax	; rs

		cmp r9,rcx
		jne uscita_prodottoScalare_ass_PS

		mov rdx,0
		mov rax,r9
		div dword[otto]
		mov rdi,rax ; in rdi ho il numero di volte che posso ciclare con ymm0 = nColonne/8

		vxorps ymm7,ymm7,ymm7

		xor rsi,rsi
		ciclo_quoziente_PS:
			cmp rsi,rdi
			jge fine_ciclo_quoziente_PS

			_1:mov rdx,0
			_11:mov rax,32
			_12:mul rsi
			_13:mov rdx,r8
			_14:mov rbx, r10

			_15:vmovups ymm0, [rdx+rax]
			_16:vmovups ymm1, [rbx+rax]
			; ho caricato gli 8 elementi dei due vettori

			_2:vmulps ymm0,ymm0,ymm1
			_22:vhaddps ymm0,ymm0,ymm0
			_23:vhaddps ymm0,ymm0,ymm0
			_24:vhaddps ymm0,ymm0,ymm0

			vaddss xmm7,xmm0

			add rsi, 1
			jmp ciclo_quoziente_PS

		fine_ciclo_quoziente_PS:
		mov rdx,0
		mov rax,8
		mul rsi
		mov rsi,rax ; in esi ho esi*4 = colonna corrente da cui iniziare l'analisi singola

		ciclo_resto_PS:
			cmp rsi,r9
			jge fine_ciclo_resto_PS

			mov rdx,0
			mov rax,4
			mul rsi
			mov rbx,r8
			mov rcx,r10

			vmovss xmm0, [ebx+eax]
			vmovss xmm1, [ecx+eax]
			; ho caricato i due elementi

			vmulss xmm0,xmm1

			vaddss xmm7,xmm0

			add rsi,1
			jmp ciclo_resto_PS

		fine_ciclo_resto_PS:
		
		vmovss [r12],xmm7

		uscita_prodottoScalare_ass_PS:

		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante




















calcolaNorma_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------

	; float* arr, int length, float* norma

	; per calcolare la norma devo fare loop unrolling quoziente su un registro ymm per length/8 volte
	; e successivamente per la differenza devo eseguire un loop resto

    mov r8, rdi ; puntatore array
    mov r9, rsi ; length array
    mov r10,rdx    ; puntatore all'elemento risultato

	vxorps xmm1,xmm1

	mov rdx,0
	mov rax,r9
	div qword[otto]
	mov rdi, rax

	xor rsi, rsi
	ciclo_quoz_CN:
		cmp rsi, rdi
		je continua_CN
		mov rbx,r8 ; prendo il puntatore all'array
		mov rdx,0
		mov rax, 32 ; moltiplico esi per 32 in maniera tale da saltare 8 elementi
		mul rsi
		vmovups ymm0, [rbx+rax] ; prendo 8 elementi
		vmulps ymm0,ymm0,ymm0
		
        vhaddps ymm0,ymm0,ymm0 	; sommo pos 0 e 1 di xmm
		vhaddps ymm0,ymm0,ymm0 	; stessa cosa
		vhaddps ymm0,ymm0,ymm0	; ho ridotto nei primi 32 bit la somma degli 8 elementi

		vaddss xmm1,xmm0 ; sommo norma a quello calcolato prima
		
		add rsi,1
		jmp ciclo_quoz_CN

	continua_CN:
		mov rdx,0
		mov rax,rsi
		mul dword[otto]
		mov rsi,rax

	ciclo_resto_CN:
		cmp rsi, r9
		jge fine_calcolaNorma_ass_CN
		
		mov rbx,r8
		vmovss xmm0,[rbx+4*rsi]
		vmulss xmm0,xmm0
		addss xmm1, xmm0
		
		add esi,1
		jmp ciclo_resto_CN

	fine_calcolaNorma_ass_CN:

	vsqrtss xmm1,xmm1
	vmovss [r10],xmm1

	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp					; ripristina il Base Pointer
	ret							; torna alla funzione C chiamante

