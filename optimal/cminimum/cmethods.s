	.file	"cmethods.c"
	.intel_syntax noprefix
	.text
	.p2align 2,,3
	.globl	_derivative
	.def	_derivative;	.scl	2;	.type	32;	.endef
_derivative:
	sub	esp, 28
	fld	DWORD PTR LC0
	fstp	QWORD PTR [esp]
	call	[DWORD PTR [esp+32]]
	add	esp, 28
	ret
	.p2align 2,,3
	.globl	_derivative1
	.def	_derivative1;	.scl	2;	.type	32;	.endef
_derivative1:
	push	ebx
	sub	esp, 72
	mov	ebx, DWORD PTR [esp+80]
	fld	QWORD PTR [esp+84]
	fld	QWORD PTR [esp+92]
	fld	st(1)
	fadd	st, st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [esp+16]
	fstp	QWORD PTR [esp+32]
	call	ebx
	fstp	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+32]
	fstp	QWORD PTR [esp]
	call	ebx
	fsubr	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+16]
	fdivp	st(1), st
	add	esp, 72
	pop	ebx
	ret
	.p2align 2,,3
	.globl	_derivative2
	.def	_derivative2;	.scl	2;	.type	32;	.endef
_derivative2:
	push	ebx
	sub	esp, 56
	mov	ebx, DWORD PTR [esp+64]
	fld	QWORD PTR [esp+68]
	fld	QWORD PTR [esp+76]
	fstp	QWORD PTR [esp+32]
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [esp+16]
	call	ebx
	fstp	QWORD PTR [esp+40]
	fld	QWORD PTR [esp+16]
	fsub	QWORD PTR [esp+32]
	fstp	QWORD PTR [esp]
	call	ebx
	fsubr	QWORD PTR [esp+40]
	fdiv	QWORD PTR [esp+32]
	add	esp, 56
	pop	ebx
	ret
	.p2align 2,,3
	.globl	_derivative3
	.def	_derivative3;	.scl	2;	.type	32;	.endef
_derivative3:
	push	ebx
	sub	esp, 72
	mov	ebx, DWORD PTR [esp+80]
	fld	QWORD PTR [esp+84]
	fld	QWORD PTR [esp+92]
	fld	st(1)
	fadd	st, st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [esp+16]
	fstp	QWORD PTR [esp+32]
	call	ebx
	fstp	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+16]
	fld	QWORD PTR [esp+32]
	fsubrp	st(1), st
	fstp	QWORD PTR [esp]
	call	ebx
	fsubr	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+16]
	fadd	st, st(0)
	fdivp	st(1), st
	add	esp, 72
	pop	ebx
	ret
	.p2align 2,,3
	.globl	_gradient1
	.def	_gradient1;	.scl	2;	.type	32;	.endef
_gradient1:
	push	ebp
	push	edi
	push	esi
	push	ebx
	sub	esp, 60
	mov	ebp, DWORD PTR [esp+80]
	mov	esi, DWORD PTR [esp+84]
	mov	edi, DWORD PTR [esp+92]
	fld	QWORD PTR [esp+96]
	fstp	QWORD PTR [esp+32]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	call	ebp
	fstp	QWORD PTR [esp+40]
	test	edi, edi
	je	L5
	xor	ebx, ebx
	.p2align 2,,3
L7:
	fld	QWORD PTR [esi+ebx*8]
	fld	QWORD PTR [esp+32]
	fadd	st, st(1)
	fstp	QWORD PTR [esi+ebx*8]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	fstp	QWORD PTR [esp+16]
	call	ebp
	fsub	QWORD PTR [esp+40]
	fdiv	QWORD PTR [esp+32]
	mov	eax, DWORD PTR [esp+88]
	fstp	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [esp+16]
	fstp	QWORD PTR [esi+ebx*8]
	inc	ebx
	cmp	ebx, edi
	jne	L7
L5:
	add	esp, 60
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 2,,3
	.globl	_gradient2
	.def	_gradient2;	.scl	2;	.type	32;	.endef
_gradient2:
	push	ebp
	push	edi
	push	esi
	push	ebx
	sub	esp, 60
	mov	ebp, DWORD PTR [esp+80]
	mov	esi, DWORD PTR [esp+84]
	mov	edi, DWORD PTR [esp+92]
	fld	QWORD PTR [esp+96]
	fstp	QWORD PTR [esp+32]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	call	ebp
	fstp	QWORD PTR [esp+40]
	test	edi, edi
	je	L10
	xor	ebx, ebx
	.p2align 2,,3
L12:
	fld	QWORD PTR [esi+ebx*8]
	fld	QWORD PTR [esp+32]
	fsubr	st, st(1)
	fstp	QWORD PTR [esi+ebx*8]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	fstp	QWORD PTR [esp+16]
	call	ebp
	fsub	QWORD PTR [esp+40]
	fdiv	QWORD PTR [esp+32]
	mov	eax, DWORD PTR [esp+88]
	fstp	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [esp+16]
	fstp	QWORD PTR [esi+ebx*8]
	inc	ebx
	cmp	ebx, edi
	jne	L12
L10:
	add	esp, 60
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 2,,3
	.globl	_gradient3
	.def	_gradient3;	.scl	2;	.type	32;	.endef
_gradient3:
	push	ebp
	push	edi
	push	esi
	push	ebx
	sub	esp, 60
	mov	ebp, DWORD PTR [esp+80]
	mov	esi, DWORD PTR [esp+84]
	mov	edi, DWORD PTR [esp+92]
	fld	QWORD PTR [esp+96]
	fst	QWORD PTR [esp+24]
	test	edi, edi
	je	L18
	xor	ebx, ebx
	fadd	st, st(0)
	fstp	QWORD PTR [esp+40]
	.p2align 2,,3
L16:
	fld	QWORD PTR [esi+ebx*8]
	fst	QWORD PTR [esp+16]
	fadd	QWORD PTR [esp+24]
	fstp	QWORD PTR [esi+ebx*8]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	call	ebp
	fstp	QWORD PTR [esp+32]
	fld	QWORD PTR [esp+16]
	fsub	QWORD PTR [esp+24]
	fstp	QWORD PTR [esi+ebx*8]
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	call	ebp
	fsub	QWORD PTR [esp+32]
	fdiv	QWORD PTR [esp+40]
	mov	eax, DWORD PTR [esp+88]
	fstp	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [esp+16]
	fstp	QWORD PTR [esi+ebx*8]
	inc	ebx
	cmp	ebx, edi
	jne	L16
	jmp	L14
L18:
	fstp	st(0)
	.p2align 2,,3
L14:
	add	esp, 60
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 2,,3
	.globl	_trapesium1
	.def	_trapesium1;	.scl	2;	.type	32;	.endef
_trapesium1:
	push	edi
	push	esi
	push	ebx
	sub	esp, 80
	mov	esi, DWORD PTR [esp+96]
	mov	edi, DWORD PTR [esp+100]
	fld	QWORD PTR [esp+104]
	fst	QWORD PTR [esp+64]
	fsub	QWORD PTR [esp+112]
	fstp	QWORD PTR [esp+48]
	xor	edx, edx
	mov	DWORD PTR [esp+72], edi
	mov	DWORD PTR [esp+76], edx
	fild	QWORD PTR [esp+72]
	fdivr	QWORD PTR [esp+48]
	fstp	QWORD PTR [esp+48]
	dec	edi
	xor	ebx, ebx
	fldz
	.p2align 2,,3
L20:
	xor	edx, edx
	mov	DWORD PTR [esp+72], ebx
	mov	DWORD PTR [esp+76], edx
	fild	QWORD PTR [esp+72]
	fmul	QWORD PTR [esp+48]
	fadd	QWORD PTR [esp+64]
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [esp+32]
	fstp	QWORD PTR [esp+16]
	call	esi
	fstp	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+32]
	fadd	QWORD PTR [esp+48]
	fstp	QWORD PTR [esp]
	call	esi
	fadd	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+16]
	faddp	st(1), st
	inc	ebx
	cmp	ebx, edi
	jbe	L20
	fld	QWORD PTR [esp+48]
	fmul	DWORD PTR LC3
	fmulp	st(1), st
	add	esp, 80
	pop	ebx
	pop	esi
	pop	edi
	ret
	.p2align 2,,3
	.globl	_trapesium2
	.def	_trapesium2;	.scl	2;	.type	32;	.endef
_trapesium2:
	push	ebp
	push	edi
	push	esi
	push	ebx
	sub	esp, 108
	mov	esi, DWORD PTR [esp+128]
	fld	QWORD PTR [esp+132]
	fstp	QWORD PTR [esp+56]
	fld	QWORD PTR [esp+140]
	fst	QWORD PTR [esp+72]
	fsubr	QWORD PTR [esp+148]
	fdiv	QWORD PTR [esp+56]
	fstp	QWORD PTR [esp]
	call	_round
	fnstcw	WORD PTR [esp+86]
	mov	ax, WORD PTR [esp+86]
	mov	ah, 12
	mov	WORD PTR [esp+84], ax
	fldcw	WORD PTR [esp+84]
	fistp	QWORD PTR [esp+88]
	fldcw	WORD PTR [esp+86]
	mov	edi, DWORD PTR [esp+88]
	dec	edi
	xor	ebx, ebx
	fldz
	.p2align 2,,3
L23:
	xor	edx, edx
	mov	DWORD PTR [esp+88], ebx
	mov	DWORD PTR [esp+92], edx
	fild	QWORD PTR [esp+88]
	fmul	QWORD PTR [esp+56]
	fadd	QWORD PTR [esp+72]
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [esp+32]
	fstp	QWORD PTR [esp+16]
	call	esi
	fstp	QWORD PTR [esp+64]
	fld	QWORD PTR [esp+32]
	fadd	QWORD PTR [esp+56]
	fstp	QWORD PTR [esp]
	call	esi
	fadd	QWORD PTR [esp+64]
	fld	QWORD PTR [esp+16]
	faddp	st(1), st
	inc	ebx
	cmp	ebx, edi
	jbe	L23
	fld	QWORD PTR [esp+56]
	fmul	DWORD PTR LC3
	fmulp	st(1), st
	add	esp, 108
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 2,,3
	.globl	_tomasAlgorithm
	.def	_tomasAlgorithm;	.scl	2;	.type	32;	.endef
_tomasAlgorithm:
	push	ebp
	push	edi
	push	esi
	push	ebx
	sub	esp, 92
	mov	edx, DWORD PTR [esp+112]
	mov	DWORD PTR [esp+48], edx
	mov	esi, DWORD PTR [esp+116]
	mov	ecx, DWORD PTR [esp+120]
	mov	DWORD PTR [esp+28], ecx
	mov	ebp, DWORD PTR [esp+124]
	mov	edi, DWORD PTR [esp+128]
	mov	DWORD PTR [esp+52], edi
	mov	edi, DWORD PTR [esp+132]
	lea	edx, [0+edi*8]
	mov	DWORD PTR [esp+40], edx
	mov	DWORD PTR [esp], edx
	call	_malloc
	mov	ebx, eax
	mov	ecx, DWORD PTR [esp+40]
	mov	DWORD PTR [esp], ecx
	call	_malloc
	test	edi, edi
	je	L26
	lea	edx, [edi-1]
	mov	DWORD PTR [esp+32], edx
	sal	edx, 3
	mov	DWORD PTR [esp+36], edx
	add	edx, ebx
	mov	DWORD PTR [esp+44], edx
	mov	ecx, DWORD PTR [esp+36]
	add	ecx, ebp
	mov	DWORD PTR [esp+56], ecx
	mov	edx, DWORD PTR [esp+48]
	add	edx, DWORD PTR [esp+36]
	mov	DWORD PTR [esp+60], edx
	lea	ecx, [-16+edi*8]
	lea	edx, [ebx+ecx]
	mov	DWORD PTR [esp+64], edx
	mov	edx, DWORD PTR [esp+36]
	add	edx, esi
	mov	DWORD PTR [esp+68], edx
	add	ecx, eax
	mov	DWORD PTR [esp+72], ecx
	mov	ecx, DWORD PTR [esp+36]
	add	ecx, eax
	mov	DWORD PTR [esp+76], ecx
	xor	ecx, ecx
	mov	edx, DWORD PTR [esp+28]
	mov	DWORD PTR [esp+28], edi
	jmp	L30
	.p2align 2,,3
L37:
	fld	QWORD PTR [ebp+0]
	fdiv	QWORD PTR [esi]
	fstp	QWORD PTR [ebx]
	fld	QWORD PTR [edx]
	fchs
	fdiv	QWORD PTR [esi]
	fstp	QWORD PTR [eax]
L28:
	inc	ecx
	cmp	ecx, DWORD PTR [esp+28]
	je	L36
L30:
	test	ecx, ecx
	je	L37
	cmp	DWORD PTR [esp+32], ecx
	je	L38
	mov	edi, DWORD PTR [esp+48]
	fld	QWORD PTR [edi+ecx*8]
	fld	QWORD PTR [ebx-8+ecx*8]
	fmul	st, st(1)
	fsubr	QWORD PTR [ebp+0+ecx*8]
	fxch	st(1)
	fmul	QWORD PTR [eax-8+ecx*8]
	fadd	QWORD PTR [esi+ecx*8]
	fdiv	st(1), st
	fxch	st(1)
	fstp	QWORD PTR [ebx+ecx*8]
	fld	QWORD PTR [edx+ecx*8]
	fchs
	fdivrp	st(1), st
	fstp	QWORD PTR [eax+ecx*8]
	inc	ecx
	cmp	ecx, DWORD PTR [esp+28]
	jne	L30
L36:
	cmp	DWORD PTR [esp+32], -1
	je	L26
	mov	esi, DWORD PTR [esp+52]
	add	esi, DWORD PTR [esp+36]
	mov	ecx, DWORD PTR [esp+52]
	add	ecx, DWORD PTR [esp+40]
	mov	edx, DWORD PTR [esp+32]
	mov	ebp, edx
	mov	edi, DWORD PTR [esp+44]
	jmp	L33
	.p2align 2,,3
L31:
	fld	QWORD PTR [eax+edx*8]
	fmul	QWORD PTR [ecx]
	fadd	QWORD PTR [ebx+edx*8]
	fstp	QWORD PTR [ecx-8]
	dec	edx
	sub	ecx, 8
	cmp	edx, -1
	je	L26
L33:
	cmp	ebp, edx
	jne	L31
	fld	QWORD PTR [edi]
	fstp	QWORD PTR [esi]
	dec	edx
	sub	ecx, 8
	cmp	edx, -1
	jne	L33
L26:
	mov	DWORD PTR [esp], ebx
	mov	DWORD PTR [esp+24], eax
	call	_free
	mov	eax, DWORD PTR [esp+24]
	mov	DWORD PTR [esp+112], eax
	add	esp, 92
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	jmp	_free
	.p2align 2,,3
L38:
	mov	edi, DWORD PTR [esp+60]
	fld	QWORD PTR [edi]
	fld	st(0)
	mov	edi, DWORD PTR [esp+64]
	fmul	QWORD PTR [edi]
	mov	edi, DWORD PTR [esp+56]
	fsubr	QWORD PTR [edi]
	fxch	st(1)
	mov	edi, DWORD PTR [esp+72]
	fmul	QWORD PTR [edi]
	mov	edi, DWORD PTR [esp+68]
	fadd	QWORD PTR [edi]
	fdivp	st(1), st
	mov	edi, DWORD PTR [esp+44]
	fstp	QWORD PTR [edi]
	mov	edi, DWORD PTR [esp+76]
	mov	DWORD PTR [edi], 0
	mov	DWORD PTR [edi+4], 0
	jmp	L28
	.section .rdata,"dr"
	.align 4
LC0:
	.long	1084227584
	.align 4
LC3:
	.long	1056964608
	.def	_round;	.scl	2;	.type	32;	.endef
	.def	_malloc;	.scl	2;	.type	32;	.endef
	.def	_free;	.scl	2;	.type	32;	.endef
