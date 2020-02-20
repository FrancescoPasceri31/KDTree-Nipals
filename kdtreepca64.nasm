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
somma_PMV: resd 8
ymmTMP: resd 8

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
global calcolaNorma_ass_64
global prodScalare_ass_64
global distanzaEuclidea_ass_64
global trasponi_ass_64
global divisioneVettoreScalare_ass_64
global prodMatrVett_ass_64
global sottrazioneMatrici_ass_64
global prodMatr_ass_64


msg	db 'n:',0
nl	db 10,0


prodMatr_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali
	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------
	
	; float* m1 +8, int nRighe1 +12, int nColonne1 +16, float* m2 +20, int nRighe2 +24, int nColonne2 +28, float* risultato +32

	mov r12, r8
	mov r13, r9
	mov r8, rdi
	mov r9, rsi
	mov r10, rdx
	mov r11, rcx
	mov r14, [rbp+16]

	cmp r10, r12
	jne fine_prodMatrici_ass_PM

	mov rdx,0
	mov rax, r9
	div dword[quattro]
	mov rcx,rax

	xor rsi,rsi
	ciclo_quoziente_righe1_PM:
		cmp rsi, rcx
		jge fine_ciclo_quoziente_righe1_PM

		push rcx 
		
		mov rdx,0
		mov rax,4
		mul rsi
		mov rsi,rax

		xor rdi,rdi
		ciclo_colonne_1_PM:
			cmp rdi, r13
			jge fine_ciclo_colonne_1_PM

			mov rdx,0
			mov rax, r10
			div dword[otto]
			mov rcx, rax ; rcx = nCol1/8

			vxorps ymm7,ymm7,ymm7

			xor rbx,rbx
			ciclo_quoziente_1_PM:
				cmp rbx,rcx
				jge fine_ciclo_quoziente_1_PM

				mov rdx,0
				mov rax,8
				mul rbx
				mov rbx,rax	; moltiplico per 8 per ottenere il numero del blocco da 8 della colonna

				; prendo elementi della prima colonna della seconda matrice
				push rcx
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP],xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+4], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+8], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+12], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+16], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+20], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+24], xmm1

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]
				vmovss [ymmTMP+28], xmm1

				vxorps ymm1,ymm1,ymm1
				vmovups ymm1,[ymmTMP]
				; in ymm1 ho la colonna della seconda matrice

				sub rbx,7
				mov rdx,0
				mov rax,rbx
				div dword[otto]
				mov rbx,rax

				push rcx
				mov rdx,0
				mov rax,4
				mul r10
				mul rsi
				mov rcx, rax
				mov rdx,0
				mov rax, 32
				mul rbx
				add rax,rcx
				pop rcx
				mov rdx, r8
				vmovups ymm0, [rdx+rax]

				push rcx
				add rsi,1
				mov rdx,0
				mov rax,4
				mul r10
				mul rsi
				mov rcx, rax
				mov rdx,0
				mov rax, 32
				mul rbx
				add rax,rcx
				pop rcx
				mov rdx,r8
				vmovups ymm2, [rdx+rax]

				push rcx
				add rsi,1
				mov rdx,0
				mov rax,4
				mul r10
				mul rsi
				mov rcx, rax
				mov rdx,0
				mov rax, 32
				mul rbx
				add rax,rcx
				pop rcx
				mov rdx,r8
				vmovups ymm3, [rdx+rax]

				push rcx
				add rsi,1
				mov rdx,0
				mov rax,4
				mul r10
				mul rsi
				mov rcx, rax
				mov rdx,0
				mov rax, 32
				mul rbx
				add rax,rcx
				pop rcx
				mov rdx,r8
				vmovups ymm4, [rdx+rax]

				sub rsi,3

				; ho ymm0,2,3,4 -> 8 elementi delle 4 righe di m1
				; ho ymm1 -> 8 elementi della colonna di m2

				vmulps ymm0,ymm0,ymm1
				vmulps ymm2,ymm2,ymm1
				vmulps ymm3,ymm3,ymm1
				vmulps ymm4,ymm4,ymm1

				vhaddps ymm0,ymm0,ymm0
				vhaddps ymm0,ymm0,ymm0

				vhaddps ymm2,ymm2,ymm2
				vhaddps ymm2,ymm2,ymm2


				vhaddps ymm3,ymm3,ymm3
				vhaddps ymm3,ymm3,ymm3

				vhaddps ymm4,ymm4,ymm4
				vhaddps ymm4,ymm4,ymm4

				vperm2f128 ymm8,ymm0,ymm0,00000001b
				vperm2f128 ymm9,ymm2,ymm2,00000001b
				vperm2f128 ymm10,ymm3,ymm3,00000001b
				vperm2f128 ymm11,ymm4,ymm4,00000001b
				; inserisco la parte alta di ymm0 in ymm8 e per gli altri

				vaddss xmm0,xmm0,xmm8
				vaddss xmm2,xmm2,xmm9
				vaddss xmm3,xmm3,xmm10
				vaddss xmm4,xmm4,xmm11

				vmovss [xmmTMP],xmm0
				vmovss [xmmTMP+4],xmm2
				vmovss [xmmTMP+8],xmm3
				vmovss [xmmTMP+12],xmm4

				vmovups xmm0,[xmmTMP]

				vaddps xmm7,xmm7,xmm0

				add rbx,1
				jmp ciclo_quoziente_1_PM

			fine_ciclo_quoziente_1_PM:

			mov rdx,0
			mov rax,8
			mul rbx
			mov rbx,rax

			ciclo_resto_1_PM:
				cmp rbx,r10
				jge inserisco_in_memoria_PM

				mov rdx,0
				mov rax,rbx
				mul r13
				add rax,rdi
				mul dword[quattro]
				; rax = (rbx*nCol2 + rdi)*4
				mov rdx, r11
				vmovss xmm1, [rdx+rax]
				; ho caricato in xmm1 l'elemento della seconda matrice

				mov rdx,0
				mov rax, r10
				mul rsi
				add rax,rbx
				mul dword[quattro]
				mov rdx, r8
				vmovss xmm0, [rdx+rax]

				add rsi,1
				mov rdx,0
				mov rax, r10
				mul rsi
				add rax,rbx
				mul dword[quattro]
				mov rdx,r8
				vmovss xmm2, [rdx+rax]

				add rsi,1
				mov rdx,0
				mov rax, r10
				mul rsi
				add rax,rbx
				mul dword[quattro]
				mov rdx, r8
				vmovss xmm3, [rdx+rax]

				add rsi,1
				mov rdx,0
				mov rax, r10
				mul rsi
				add rax,rbx
				mul dword[quattro]
				mov rdx,r8
				vmovss xmm4, [rdx+rax]
				; nei primi 32 bit di xmm0,2,3,4 ho caricato gli elementi delle righe

				sub rsi,3

				vmulss xmm0,xmm0,xmm1
				vmulss xmm2,xmm2,xmm1
				vmulss xmm3,xmm3,xmm1
				vmulss xmm4,xmm4,xmm1

				vshufps xmm0, xmm0,xmm2,00000000b ;[aabb]
				vshufps xmm3, xmm3, xmm4, 00000000b ;[ccdd]

				vshufps xmm1, xmm0, xmm3, 10001000b ;[abcd]

				vaddps xmm7,xmm7,xmm1

				add rbx,1
				jmp ciclo_resto_1_PM
				
			inserisco_in_memoria_PM:
			
			mov rdx,0
			mov rax, r13
			mul rsi
			add rax,rdi
			mov rdx,0
			mul dword[quattro]
			mov rdx,r14
			vmovss [rdx+rax],xmm7
			
			vmovups [xmmTMP],xmm7

			add rsi,1
			mov rdx,0
			mov rax,r13
			mul rsi
			add rax,rdi
			mov rdx,0
			mul dword[quattro]
			mov rdx,r14
			vmovss xmm0, [xmmTMP+4]
			vmovss [rdx+rax],xmm0

			add rsi,1
			mov rdx,0
			mov rax,r13
			mul rsi
			add rax,rdi
			mov rdx,0
			mul dword[quattro]
			mov rdx,r14
			vmovss xmm1, [xmmTMP+8]
			vmovss [rdx+rax],xmm1

			add rsi,1
			mov rdx,0
			mov rax,r13
			mul rsi
			add rax,rdi
			mul dword[quattro]
			mov rdx,r14
			vmovss xmm2, [xmmTMP+12]
			vmovss [rdx+rax],xmm2

			sub rsi,3

			add rdi,1
			jmp ciclo_colonne_1_PM

		fine_ciclo_colonne_1_PM:


		mov rdx,0
		mov rax,rsi
		div dword[quattro]
		mov rsi,rax

		pop rcx
		add rsi,1
		jmp ciclo_quoziente_righe1_PM

	fine_ciclo_quoziente_righe1_PM:
	
	; abbiamo analizzato le righe a 4 alla volta ora dobbiamo calcolare 1 alla volta

	mov rdx,0
	mov rax,4
	mul rsi
	mov rsi,rax ; ottengo in rsi la riga da analizzare

	ciclo_resto_righe1_PM:
		cmp rsi,r9
		jge fine_prodMatrici_ass_PM

		xor rdi,rdi
		ciclo_colonne_2_PM:
			cmp rdi, r13
			jge fine_ciclo_colonne_2_PM

			mov rdx,0
			mov rax, r10
			div dword[otto]
			mov rcx, rax ; rcx = nCol1/8

			vxorps ymm9,ymm9,ymm9

			xor rbx,rbx
			ciclo_quoziente_2_PM:
				cmp rbx,rcx
				jge fine_ciclo_quoziente_2_PM

				mov rdx,0
				mov rax,8
				mul rbx
				mov rbx,rax
				
				; prendo elementi della prima colonna della seconda matrice
				push rcx
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm1, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm2, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm3, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm4, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm5, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm6, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm7, [rdx+rax]

				push rcx
				add rbx,1
				mov rdx,0
				mov rax,4
				mul r13
				mul rbx
				mov rcx,rax
				mov rdx,0
				mov rax,4
				mul rdi
				add rax,rcx
				pop rcx
				mov rdx,r11
				vmovss xmm8, [rdx+rax]

				; xmmi [x x x i]

				vshufps xmm1,xmm1,xmm2,0	; b b a a
				vshufps xmm3,xmm3,xmm4,0	; d d c c
				vshufps xmm5,xmm5,xmm6,0	; f f e e
				vshufps xmm7,xmm7,xmm8,0	; h h g g

				vshufps xmm1,xmm1,xmm3,10001000b ; [d c b a]
				vshufps xmm5,xmm5,xmm7,10001000b ; [h g f e]

				vperm2f128 ymm1,ymm1,ymm5,00100000b	; [h g f e d c b a] 

				sub rbx,7
				mov rdx,0
				mov rax,rbx
				div dword[otto]
				mov rbx,rax

				push rcx
				mov rdx,0
				mov rax,4
				mul r10
				mul rsi
				mov rcx, rax	; rcx= indRiga1*nCol1*4
				mov rdx,0
				mov rax, 32
				mul rbx
				add rax,rcx	; rax = indRiga1*nCol1*4 + indCol*32
				pop rcx
				mov rdx, r8
				vmovups ymm0, [rdx+rax]
				; ho ymm0 -> 8 elementi delle 4 righe di m1
				; ho ymm1 -> 8 elementi della colonna di m2

				vmulps ymm0,ymm0,ymm1

				vhaddps ymm0,ymm0,ymm0
				vhaddps ymm0,ymm0,ymm0

				vperm2f128 ymm8,ymm0,ymm0,00000001b

				vaddss xmm0,xmm0,xmm8

				vaddss xmm9,xmm9,xmm0 ; sommo in xmm7 il contributo

				add rbx,1
				jmp ciclo_quoziente_2_PM

			fine_ciclo_quoziente_2_PM:

			mov rdx,0
			mov rax,8
			mul rbx
			mov rbx,rax

			ciclo_resto_2_PM:
				cmp rbx,r10
				jge inserisco_in_memoria_2_PM

				mov rdx,0
				mov rax,rbx
				mul r13
				add rax,rdi
				mul dword[quattro]
				; rax = (rbx*nCol2 + rdi)*4
				mov rdx, r11
				vmovss xmm1, [rdx+rax]
				; ho caricato in xmm1 l'elemento della seconda matrice


				mov rdx,0
				mov rax, r10
				mul rsi
				add rax,rbx
				mul dword[quattro]
				mov rdx, r8
				vmovss xmm0, [rdx+rax]
				; nei primi 32 bit di xmm0 ho caricato l' elemento della righa

				vmulss xmm0,xmm0,xmm1

				vaddss xmm9,xmm9,xmm0

				add rbx,1
				jmp ciclo_resto_2_PM
				
			inserisco_in_memoria_2_PM:
			
			mov rdx,0
			mov rax,r13
			mul rsi
			add rax,rdi
			mov rdx,0
			mul dword[quattro]
			mov rdx,r14
			vmovss [rdx+rax],xmm9
			
			add rdi,1
			jmp ciclo_colonne_2_PM

		fine_ciclo_colonne_2_PM:

		add rsi,1
		jmp ciclo_resto_righe1_PM

	
	fine_prodMatrici_ass_PM:
	
    ; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp					; ripristina il Base Pointer
	ret							; torna alla funzione C chiamante



sottrazioneMatrici_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------

	; float* m1 +8, float* m2 +12, int nRighe1 +16, int nColonne1 +20, float* risultato +32

	mov r14, r8
	mov r8,rdi
	mov r9, rsi
	mov r10, rdx
	mov r11, rcx

	mov rdx,0
	mov rax, r10
	div dword[quattro]
	mov rdi, rax ; ho ottenuto il nRighe/4 in RDI

	xor rsi,rsi
	ciclo_esterno_1_SM:
		cmp rsi,rdi
		jge fine_ciclo_esterno_1_SM
		; altrimenti prendo 4 righe alla volta

		mov rdx,0
		mov rax, r11
		div dword[otto]	
		mov rbx,rax	; in rbx ho nColonne/8

		mov r15, rdi

		mov rdx,0
		mov rax,4
		mul rsi
		mov rsi,rax
		; inserisco in RSI = RSI*4 per trovare il primo indice riga delle 4 da analizzare

		xor rcx,rcx
		ciclo_interno_1_SM:
			cmp rcx,rbx
			jge fine_ciclo_interno_1_SM

			; in ymm0-1-2-3 mantengo la prima matrice
			; in ymm4-5-6-7 mantengo la seconda matrice

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (32 byte*indiceColonna) [rsi*r11*4 + 32*ecx]
			mov rdx, r8
			vmovups ymm0, [rdx+rax]	; ho caricato la prima riga
			mov rdx, r9
			vmovups ymm4, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm1, [rdx+rax]
			mov rdx, r9
			vmovups ymm5, [rdx+rax]
			
			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm2, [rdx+rax]
			mov rdx, r9
			vmovups ymm6, [rdx+rax]
			
			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm3, [rdx+rax]
			mov rdx, r9
			vmovups ymm7, [rdx+rax]
			; ho caricato le 4 righe in ymm0,ymm1,ymm2,ymm3
			; e della seconda matrice in ymm4,5,6,7

			vsubps ymm0,ymm0,ymm4
			vsubps ymm1,ymm1,ymm5
			vsubps ymm2,ymm2,ymm6
			vsubps ymm3,ymm3,ymm7
			; ho sottratto gli elementi e devo inserirli in memoria

			sub rsi,3
			mov rdx,0
			mov rax,rsi
			mov rsi,rax			

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (16 byte*indiceColonna) [esi*[ebp+20]*4 + 16*ecx]
			mov rdx, r14
			vmovups [rdx+rax], ymm0	; ho caricato la prima riga

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi,rax
			mov rdx,0
			mov rax,32
			mul rcx
			add rax,rdi
			mov rdx, r14
			vmovups [rdx+rax],ymm1

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi,rax
			mov rdx,0
			mov rax,32
			mul rcx
			add rax,rdi
			mov rdx, r14
			vmovups [rdx+rax],ymm2

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi,rax
			mov rdx,0
			mov rax,32
			mul rcx
			add rax,rdi
			mov rdx, r14
			vmovups [rdx+rax],ymm3	
			; ho caricato le 4 righe da xmm0,xmm1,xmm2,xmm3
			; in memoria

			sub rsi,3
			; riporto esi al valore del ciclo esterno 1

			add rcx,1
			jmp ciclo_interno_1_SM

		fine_ciclo_interno_1_SM:
		
		xor rdi,rdi
		; ci manca il ciclo interno 2 in cui cicliamo sulle singole colonne dopo aver
		; controllato che non abbiamo raggiunto il nColonne giusto (ecx*4 == nColonne)
		mov rdx,0
		mov rax,8
		mul rcx
		mov rcx,rax

		ciclo_interno_2_SM:
			cmp rcx,r11
			jge concluso_4_righe_SM

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm0, [rdx+rax]
			mov rdx,r9
			vmovss xmm4, [rdx+rax]
 
			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm1, [rdx+rax]
			mov rdx, r9
			vmovss xmm5, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm2, [rdx+rax]
			mov rdx, r9
			vmovss xmm6, [rdx+rax]
			
			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm3, [rdx+rax]
			mov rdx, r9
			vmovss xmm7, [rdx+rax]
			; nei primi 32 bit di xmm0,xmm1,xmm2,xmm3 ho i valori delle colonne
			
			vsubss xmm0,xmm0,xmm4
			vsubss xmm1,xmm1,xmm5
			vsubss xmm2,xmm2,xmm6
			vsubss xmm3,xmm3,xmm7
			; ho diviso tutti gli elementi per lo scalare

			sub rsi,3
			mov rdx,0
			mov rax,rsi
			mov rsi,rax

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r14
			vmovss [rdx+rax],xmm0
 
			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in ebx ho ecx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r14
			vmovss [rdx+rax],xmm1

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in ebx ho ecx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r14
			vmovss [rdx+rax],xmm2

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in ebx ho ecx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r14
			vmovss [rdx+rax],xmm3

			sub rsi,3
			mov rdx,0
			mov rax,rsi
			mov rsi,rax

			add ecx,1
			jmp ciclo_interno_2_SM

		concluso_4_righe_SM:
		mov rdi, r15

		mov rdx,0
		mov rax,rsi
		div dword[quattro]
		mov rsi,rax

		add rsi,1
		jmp ciclo_esterno_1_SM

	fine_ciclo_esterno_1_SM:
	mov rdx,0
	mov rax,rsi
	mul dword[quattro]
	mov rsi,rax
	

	ciclo_esterno_2_SM:	; svolgiamo una riga alla volta se non ho fatto tutte le righe con 4 alla volta
		cmp rsi,r10
		jge uscita_sottrazioneMatrice_ass_SM
		
		mov r15,rdi

		; in rsi ho l'indice riga
		mov rdx,0
		mov rax, r11
		div dword[otto] ; eax = nColonne/4
		mov rbx, rax

		vxorps ymm7,ymm7,ymm7

		xor rcx,rcx
		ciclo_quoziente_2_SM:
			cmp rcx,rbx
			jge fine_ciclo_quoziente_2_SM

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (32 byte*indiceColonna) [rsi*r11*4 + 32*rcx]
			mov rdx, r8
			vmovups ymm0, [rdx+rax]	; ho caricato la prima riga
			mov rdx, r9
			vmovups ymm4, [rdx+rax]

			vsubps ymm0,ymm0,ymm4

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (32 byte*indiceColonna) [rsi*r11*4 + 32*rcx]
			mov rdx, r14
			vmovups [rdx+rax],ymm0	; ho caricato la prima riga

			add rcx,1
			jmp ciclo_quoziente_2_SM

		fine_ciclo_quoziente_2_SM:
		; ho terminato le colonne prese a 4 ora conto le singole

		mov rdx,0
		mov rax,8
		mul rcx
		mov rcx,rax

		ciclo_resto_2_SM:
			cmp rcx,r11
			jge fine_ciclo_resto_2_SM

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna	
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov rdx, r8
			vmovss xmm0, [rdx+rax]
			mov rdx, r9
			vmovss xmm4, [rdx+rax]
			; nei primi 32 di xmm0 ho inserito l'elemento della matrice

			vsubss xmm0,xmm0,xmm4

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in ebx ho ecx*4 che è un indice di colonna	
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [esi*eax + ebx]
			mov rdx, r14
			vmovss [rdx+rax],xmm0

			add rcx,1
			jmp ciclo_resto_2_SM

		fine_ciclo_resto_2_SM:
		
		mov rdi, r15

		add rsi,1
		jmp ciclo_esterno_2_SM

	uscita_sottrazioneMatrice_ass_SM:
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp					; ripristina il Base Pointer
	ret							; torna alla funzione C chiamante





prodMatrVett_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------


	; elaborazione

	; float* m = ebp+8 , 
	; float* v = ebp+12, 
	; int nRighe = ebp+16, 
	; int nColonne = ebp+20,
	; int lengthVettore = ebp+24,
	; float* risultato = ebp+28

	mov r12, r8
	mov r13, r9
	mov r8, rdi
	mov r9, rsi
	mov r10, rdx
	mov r11, rcx


	cmp r11,r12
	jne uscita_prodMatrVett_ass_PMV

	mov rdx,0
	mov rax, r10
	div dword[otto]
	mov rdi, rax ; ho ottenuto il nRighe/8 in RDI

	xor rsi,rsi
	ciclo_esterno_1_PMV:
		cmp rsi,rdi
		jge fine_ciclo_esterno_1_PMV
		; altrimenti prendo 4 righe alla volta

		mov rdx,0
		mov rax, r11
		div dword[otto]	
		mov rbx,rax	; in rbx ho nColonne/8

		vxorps ymm9,ymm9,ymm9

		vmovups [somma_PMV],ymm9 ; sum = 0.0

		mov r15, rdi

		mov rdx,0
		mov rax,8
		mul rsi
		mov rsi,rax
		; inserisco in RSI = RSI*8 per trovare il primo indice riga delle 8 da analizzare

		xor rcx,rcx
		ciclo_interno_1_PMV:
			cmp rcx,rbx
			jge fine_ciclo_interno_1_PMV

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (32 byte*indiceColonna) [rsi*r20*4 + 32*rcx]
			mov rdx, r8
			vmovups ymm0, [rdx+rax]	; ho caricato la prima riga

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm1, [rdx+rax]
	
			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm2, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm3, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm4, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm5, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm6, [rdx+rax]

			add rsi,1
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			mov rdx, r8
			vmovups ymm7, [rdx+rax]
			; ho caricato le 8 righe in ymm0,ymm1,ymm2,ymm3, ymm4, ymm5, ymm6, ymm7

			mov rdx,0
			mov rax, 32
			mul rcx
			; in rax ho rcx*32 = indice del vettore
			mov rdx, r9
			vmovups ymm8, [rdx+rax]
			; ho caricato il vettore in ymm4

			vmulps ymm0,ymm0,ymm8
			vmulps ymm1,ymm1,ymm8
			vmulps ymm2,ymm2,ymm8
			vmulps ymm3,ymm3,ymm8
			vmulps ymm4,ymm4,ymm8
			vmulps ymm5,ymm5,ymm8
			vmulps ymm6,ymm6,ymm8
			vmulps ymm7,ymm7,ymm8

			vxorps ymm8,ymm8,ymm8
			vmovups [ymmTMP],ymm8

			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm1,ymm1,ymm1
			vhaddps ymm2,ymm2,ymm2
			vhaddps ymm3,ymm3,ymm3
			vhaddps ymm4,ymm4,ymm4
			vhaddps ymm5,ymm5,ymm5
			vhaddps ymm6,ymm6,ymm6
			vhaddps ymm7,ymm7,ymm7

			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm1,ymm1,ymm1
			vhaddps ymm2,ymm2,ymm2
			vhaddps ymm3,ymm3,ymm3
			vhaddps ymm4,ymm4,ymm4
			vhaddps ymm5,ymm5,ymm5
			vhaddps ymm6,ymm6,ymm6
			vhaddps ymm7,ymm7,ymm7

			; abbiamo liberi i registri ymm 8->15
			vperm2f128 ymm8,ymm0,ymm0,00000001b
			vperm2f128 ymm9,ymm1,ymm1,00000001b
			vperm2f128 ymm10,ymm2,ymm2,00000001b
			vperm2f128 ymm11,ymm3,ymm3,00000001b
			vperm2f128 ymm12,ymm4,ymm4,00000001b
			vperm2f128 ymm13,ymm5,ymm5,00000001b
			vperm2f128 ymm14,ymm6,ymm6,00000001b
			vperm2f128 ymm15,ymm7,ymm7,00000001b
			
			vaddss xmm0,xmm8
			vaddss xmm1,xmm9
			vaddss xmm2,xmm10
			vaddss xmm3,xmm11
			vaddss xmm4,xmm12
			vaddss xmm5,xmm13
			vaddss xmm6,xmm14
			vaddss xmm7,xmm15
			; nei 32 bit di ognuno ho la somma degli elementi della riga ymm0,ymm1,....,ymm7

			vmovss [ymmTMP],xmm0
			vmovss [ymmTMP+4],xmm1
			vmovss [ymmTMP+8],xmm2
			vmovss [ymmTMP+12],xmm3
			vmovss [ymmTMP+16],xmm4
			vmovss [ymmTMP+20],xmm5
			vmovss [ymmTMP+24],xmm6
			vmovss [ymmTMP+28],xmm7
			; ho ottenuto il vettore in mem con le 8 somme [a,b,c,d,e,f,g,h]

			vxorps ymm1, ymm1, ymm1
			vxorps ymm2, ymm2, ymm2
			vmovups ymm1, [ymmTMP]
			vmovups ymm2, [somma_PMV]
			vaddps ymm1,ymm1,ymm2
			vmovups [somma_PMV],ymm1

			sub rsi,7
			mov rdx,0
			mov rax,rsi
			mov rsi,rax
			; riporto esi al valore del ciclo esterno 1


			add rcx,1
			jmp ciclo_interno_1_PMV

		fine_ciclo_interno_1_PMV:
		
		xor rdi,rdi
		; ci manca il ciclo interno 2 in cui cicliamo sulle singole colonne dopo aver
		; controllato che non abbiamo raggiunto il nColonne giusto (rcx*8 == nColonne)
		mov rdx,0
		mov rax,8
		mul rcx
		mov rcx,rax

		ciclo_interno_2_PMV:
			cmp rcx,r11
			jge concluso_4_righe_PMV

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm0, [rdx+rax]
 
			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm1, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm2, [rdx+rax]


			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm3, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm4, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm5, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm6, [rdx+rax]

			xor rax,rax
			add rsi,1
			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in rbx ho rcx*4 che è un indice di colonna		
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm7, [rdx+rax]
			; nei primi 32 bit di xmm0,xmm1,xmm2,xmm3 ho i valori delle colonne
			
			mov rdx,0
			mov rax,4
			mul rcx
			mov rdx,r9
			vmovss xmm8, [rdx+rax]
			; ho in xmm8 l'elemento del vettore

			vmulss xmm0,xmm8
			vmulss xmm1,xmm8
			vmulss xmm2,xmm8
			vmulss xmm3,xmm8
			vmulss xmm4,xmm8
			vmulss xmm5,xmm8
			vmulss xmm6,xmm8
			vmulss xmm7,xmm8
			; abbiamo tutto moltiplicato
			; ho nei 32 bit lower tutte le 8 moltiplicazioni dei registri

			vaddss xmm0, [somma_PMV]
			vaddss xmm1, [somma_PMV+4]
			vaddss xmm2, [somma_PMV+8]
			vaddss xmm3, [somma_PMV+12]
			vaddss xmm4, [somma_PMV+16]
			vaddss xmm5, [somma_PMV+20]
			vaddss xmm6, [somma_PMV+24]
			vaddss xmm7, [somma_PMV+28]

			vmovss [somma_PMV], xmm0
			vmovss [somma_PMV+4], xmm1
			vmovss [somma_PMV+8], xmm2
			vmovss [somma_PMV+12], xmm3
			vmovss [somma_PMV+16], xmm4
			vmovss [somma_PMV+20], xmm5
			vmovss [somma_PMV+24], xmm6
			vmovss [somma_PMV+28], xmm7

			sub rsi,7
			mov rdx,0
			mov rax,rsi
			mov rsi,rax
			; riporto rsi al valore iniziale del ciclo

			add rcx,1
			jmp ciclo_interno_2_PMV

		concluso_4_righe_PMV:
		mov rdi, r15

		mov rdx,0
		mov rax, rsi
		div dword[otto]
		mov rsi,rax

		mov rbx,r13 ; prendo il puntatore del vettore risultato
		mov rdx,0
		mov rax,32
		mul rsi

		vmovups ymm7, [somma_PMV]
		vmovups [rbx+rax],ymm7

		add rsi,1
		jmp ciclo_esterno_1_PMV

	fine_ciclo_esterno_1_PMV:
	mov rdx,0
	mov rax,rsi
	mul dword[otto]
	mov rsi,rax
	
	ciclo_esterno_2_PMV:	; svolgiamo una riga alla volta se non ho fatto tutte le righe considerandone 8 alla volta
		cmp rsi,r10
		jge uscita_prodMatrVett_ass_PMV
		
		mov r15,rdi

		; in rsi ho l'indice riga
		mov rdx,0
		mov rax, r11
		div dword[otto] ; eax = nColonne/8
		mov rbx, rax

		vxorps ymm9,ymm9, ymm9

		xor rcx,rcx
		ciclo_quoziente_2_PMV:
			cmp rcx,rbx
			jge fine_ciclo_quoziente_2_PMV

			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			xor rdi,rdi
			mov rdi, rax
			mov rdx,0
			mov rax, 32
			mul rcx
			add rax,rdi
			; ho ottenuto indiceRiga*(nColonne*4 byte) + (32 byte*indiceColonna) [rsi*r11*4 + 32*rcx]
			mov rdx, r8
			vmovups ymm0, [rdx+rax]	; ho caricato la prima riga

			mov rdx,0
			mov rax, 32
			mul rcx
			; in rax ho rcx*32 = indice del vettore
			mov rdx, r9
			vmovups ymm4, [rdx+rax]
			; ho caricato il vettore in ymm4

			vmulps ymm0,ymm0,ymm4

			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm0,ymm0,ymm0

			vperm2f128 ymm1,ymm0,ymm0,00000001b
			vaddss xmm0,xmm1
			; nei primi 32 bit ho la somma degli 8 elementi moltiplcati per gli 8 del vettore

			vaddss xmm9,xmm9,xmm0

			add rcx,1
			jmp ciclo_quoziente_2_PMV

		fine_ciclo_quoziente_2_PMV:
		; ho terminato le colonne prese a 8, ora conto le singole

		mov rdx,0
		mov rax,8
		mul rcx
		mov rcx,rax

		ciclo_resto_2_PMV:
			cmp rcx,r11
			jge fine_ciclo_resto_2_PMV

			mov rdx,0
			mov rax,4
			mul rcx
			mov rbx, rax ; in ebx ho ecx*4 che è un indice di colonna	
			mov rdx,0
			mov rax, r11
			mul dword[quattro]
			mul rsi
			add rax,rbx
			; ho calcolato indiceRiga*(nColonne*4)+indiceColonna [rsi*rax + rbx]
			mov rdx, r8
			vmovss xmm0, [rdx+rax]
			; nei primi 32 di xmm0 ho inserito l'elemento della matrice

			mov rdx,0
			mov rax,4
			mul rcx
			mov rdx,r9
			vmovss xmm4, [rdx+rax]
			; ho in xmm4 l'elemento del vettore
			
			vmulss xmm0,xmm0,xmm4

			vaddss xmm9,xmm9,xmm0

			add rcx,1
			jmp ciclo_resto_2_PMV

		fine_ciclo_resto_2_PMV:
		
		; prendere i primi 32 bit di xmm9 ed aggiungerli in posizione del vettore risultato
		mov rdx,0
		mov rax,4
		mul rsi
		
		mov rcx,r13
		vmovss [rcx+rax],xmm9

		mov rdi,r15

		add rsi,1
		jmp ciclo_esterno_2_PMV

	uscita_prodMatrVett_ass_PMV:
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante






trasponi_ass_64:
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; ------------------------------------------------------------
	; I parametri sono passati nei registri
	; ------------------------------------------------------------

	; float* m +8, int nRighe +12, int nColonne +16, float* risultato +20

	mov r8, rdi
	mov r9, rsi
	mov r10, rdx
	mov r11, rcx

	xor rsi,rsi
	ciclo_righe_TRA:
		cmp rsi,r9
		jge fine_ciclo_righe_TRA

		mov rdx,0
		mov rbx,16
		mov rax,r10
		div rbx
		mov rbx,rax ; in rbx = nCol/16

		xor rdi,rdi
		ciclo_quoziente_TRA:
			cmp rdi,rbx
			jge fine_ciclo_quoziente_TRA

			mov rdx,0
			mov rax,2
			mul rdi
			mov rdi,rax

			mov rdx,0
			mov rax,r10
			mul dword[quattro]
			mul rsi
			mov rcx,rax ; rcx = rax = nCol*4*rsi
			mov rdx,0
			mov rax,32
			mul rdi
			add rax,rcx ; rax = nCol*4*rsi + rdi*32 [esi indice riga, edi indice colonna]
			mov rdx,r8
			vmovups ymm0,[rdx+rax]

			add rdi,1
			mov rdx,0
			mov rax,32
			mul rdi
			add rax,rcx
			mov rdx,r8
			vmovups ymm1,[rdx+rax]
			; in ymm0 e ymm1 ho i primi 16 elementi della riga

			sub rdi,1

			mov rdx,0
			mov rax,8
			mul rdi
			mov rdi,rax ; moltiplico rdi * 8 per andare a riempire la riga giusta

			vmovups ymm2,ymm0
			vmovups ymm3,ymm0
			vmovups ymm4,ymm0
			vmovups ymm5,ymm0
			vmovups ymm6,ymm0
			vmovups ymm7,ymm0
			vmovups ymm8, ymm0

			vmovups ymm9,ymm1
			vmovups ymm10,ymm1
			vmovups ymm11,ymm1
			vmovups ymm12,ymm1
			vmovups ymm13,ymm1
			vmovups ymm14,ymm1
			vmovups ymm15, ymm1

			vperm2f128 ymm5,ymm5,ymm5,00010001b	
			vperm2f128 ymm6,ymm6,ymm6,00010001b
			vperm2f128 ymm7,ymm7,ymm7,00010001b
			vperm2f128 ymm8,ymm8,ymm8,00010001b

			vshufps ymm2,ymm2,ymm0,1
			vshufps ymm3,ymm3,ymm0,2
			vshufps ymm4,ymm4,ymm0,3
			vshufps ymm5,ymm5,ymm5,4
			
			vshufps ymm6,ymm6,ymm6,5
			vshufps ymm7,ymm7,ymm7,6
			vshufps ymm8,ymm8,ymm8,7

			vshufps ymm9,ymm9,ymm1,1
			vshufps ymm10,ymm10,ymm1,2
			vshufps ymm11,ymm11,ymm1,3

			vperm2f128 ymm12,ymm12,ymm12,00010001b
			vperm2f128 ymm13,ymm13,ymm13,00010001b
			vperm2f128 ymm14,ymm14,ymm14,00010001b
			vperm2f128 ymm15,ymm15,ymm15,00010001b
			vshufps ymm12,ymm12,ymm12,4
			vshufps ymm13,ymm13,ymm13,5
			vshufps ymm14,ymm14,ymm14,6
			vshufps ymm15,ymm15,ymm15,7
			; in ogni registro nei primi 32 bit ho gli elementi della colonna della trasposta

			mov rdx,0
			mov rax,4
			mul rsi
			mov rcx,rax ; rcx = rsi*8

			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm0			

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm2
			
			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm3

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm4

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm5

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm6

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm7

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm8

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm1

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm9

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm10

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm11

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm12

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm13

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm14

			add rdi,1
			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm15

			sub rdi,15

			mov rdx,0
			mov rcx,16
			mov rax,rdi
			div rcx
			mov rdi,rax

			add rdi,1
			jmp ciclo_quoziente_TRA
		
		fine_ciclo_quoziente_TRA:
		mov rdx,0
		mov rax,16
		mul rdi
		mov rdi,rax

		ciclo_resto_TRA:
			cmp rdi,r10
			jge fine_ciclo_TRA
			
			mov rdx,0
			mov rax,r10
			mul dword[quattro]
			mul rsi
			mov rcx,rax ; ecx = eax = nCol*4*esi
			mov rdx,0
			mov rax,4
			mul rdi
			add rax,rcx ; eax = nCol*4*esi + edi*4 [esi indice riga, edi indice colonna]
			mov rdx,r8
			vmovss xmm0,[rdx+rax]

			mov rdx,0
			mov rax,4
			mul rsi
			mov rcx,rax ; ecx = esi*4

			mov rdx,0
			mov rax,r9
			mul rdi
			mul dword[quattro]
			add rax,rcx
			mov rdx,r11
			vmovss [rdx+rax],xmm0

			add rdi,1
			jmp ciclo_resto_TRA

		fine_ciclo_TRA:
		
		add rsi,1
		jmp ciclo_righe_TRA

	fine_ciclo_righe_TRA:
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

		vdivss xmm2,xmm2,xmm1

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
			;vhaddps ymm0,ymm0,ymm0

			vperm2f128 ymm8,ymm0,ymm0,00000001b

			vaddss xmm0,xmm0,xmm8

			vaddss xmm7,xmm7,xmm0

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
		vmovss [r11],xmm7

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

		mov r12,r8	; rs
		mov r8,rdi	; v1
		mov r9, rsi	; dim1
		mov r10, rdx	; v2
		mov r11, rcx	; dim2

		cmp r9,r11
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

			mov rdx,0
			mov rax,32
			mul rsi
			mov rdx,r8
			mov rbx, r10

			vmovups ymm0, [rdx+rax]
			vmovups ymm1, [rbx+rax]
			; ho caricato gli 8 elementi dei due vettori

			vmulps ymm0,ymm0,ymm1

			vhaddps ymm0,ymm0,ymm0
			vhaddps ymm0,ymm0,ymm0
			;vhaddps ymm0,ymm0,ymm0

			vperm2f128 ymm8,ymm0,ymm0,00000001b

			vaddss xmm0,xmm0,xmm8

			vaddss xmm7,xmm7,xmm0

			add rsi, 1
			jmp ciclo_quoziente_PS

		fine_ciclo_quoziente_PS:
		mov rdx,0
		mov rax,8
		mul rsi
		mov rsi,rax ; in rsi ho rsi*8 = colonna corrente da cui iniziare l'analisi singola

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

			vmulss xmm0,xmm0,xmm1

			vaddss xmm7,xmm7,xmm0

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
	div dword[otto]
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
		;vhaddps ymm0,ymm0,ymm0	; ho ridotto nei primi 32 bit la somma degli 8 elementi

		vperm2f128 ymm8, ymm0, ymm0,00000001b

		vaddss xmm0,xmm0,xmm8

		vaddss xmm1,xmm1,xmm0 ; sommo norma a quello calcolato prima
		
		add rsi,1
		jmp ciclo_quoz_CN

	continua_CN:
		mov rdx,0
		mov rax,8
		mul rsi
		mov rsi,rax

	ciclo_resto_CN:
		cmp rsi, r9
		jge fine_calcolaNorma_ass_CN
		
		mov rbx,r8
		vmovss xmm0,[rbx+4*rsi]
		vmulss xmm0,xmm0,xmm0
		vaddss xmm1, xmm1,xmm0
		
		add esi,1
		jmp ciclo_resto_CN

	fine_calcolaNorma_ass_CN:

	vsqrtss xmm1,xmm1,xmm1
	vmovss [r10],xmm1

	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
		
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp					; ripristina il Base Pointer
	ret							; torna alla funzione C chiamante


