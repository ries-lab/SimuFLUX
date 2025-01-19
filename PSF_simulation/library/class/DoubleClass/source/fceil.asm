;Optimized double class functions for MATLAB on x86 computers.
;Copyright � Marcel Leutenegger, 2003-2008, �cole Polytechnique F�d�rale de Lausanne (EPFL),
;Laboratoire d'Optique Biom�dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
;
;    This library is free software; you can redistribute it and/or modify it under
;    the terms of the GNU Lesser General Public License as published by the Free
;    Software Foundation; version 2.1 of the License.
;
;    This library is distributed in the hope that it will be useful, but WITHOUT ANY
;    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
;    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
;
;    You should have received a copy of the GNU Lesser General Public License along
;    with this library; if not, write to the Free Software Foundation, Inc.,
;    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

[segment .data align=16]
fpcw:	dw	0x1B7F


[segment .text align=16]
global	_fceil


;void fceil(double* or, double* oi, const double* sr, const double* si, int n)
;
_fceil:
	push ebx
	fldcw [fpcw]
	mov ecx,[esp+24]	; n
	mov ebx,[esp+16]	; sr
	sub ecx,byte 1
	mov edx,[esp+8]		; or
.1:	fld qword [ebx+ecx*8]
	fxam
	fnstsw ax
	sahf
       jc	.2
	frndint
.2:	fstp qword [edx+ecx*8]
	sub ecx,byte 1
       jns	.1
	mov ecx,[esp+24]	; n
	mov ebx,[esp+20]	; si
	sub ecx,byte 1
	test ebx,ebx
	mov edx,[esp+12]	; oi
       jz	.5
.3:	fld qword [ebx+ecx*8]
	fxam
	fnstsw ax
	sahf
       jc	.4
	frndint
.4:	fstp qword [edx+ecx*8]
	sub ecx,byte 1
       jns	.3
.5:	pop ebx
       retn
