' *****************************************************************************
'Subject: FreeBasic include file for the GMP functions called by cfr_lib
'Code   : GMP library 6.2.1 with FreeBasic 1.10.0

#pragma once

#ifdef __FB_64BIT__
  #inclib "gmp-64"
#else
  #inclib "gmp-32"
#endif

#include once "crt/long.bi"
#include once "crt/stddef.bi"

' *****************************************************************************
extern "C"

type mp_limb_t as culong
type mp_bitcnt_t as culong

type __mpz_struct
   _mp_alloc as long      ' the number of limbs currently allocated at _mp_d
   _mp_size as long       ' the number of limbs used
   _mp_d as mp_limb_t ptr ' a pointer to an array of limbs
end type

type mpz_srcptr as const __mpz_struct ptr
type mpz_ptr as __mpz_struct ptr
type mpz_t as __mpz_struct

#if defined(__FB_WIN32__) and defined(LIBGMP_DLL)
   const __GMP_LIBGMP_DLL = 1
   extern import mp_bits_per_limb alias "__gmp_bits_per_limb" as const long
   extern import gmp_errno alias "__gmp_errno" as long
   extern import gmp_version alias "__gmp_version" as const zstring const ptr
#else
   const __GMP_LIBGMP_DLL = 0
   extern mp_bits_per_limb alias "__gmp_bits_per_limb" as const long
   extern gmp_errno alias "__gmp_errno" as long
   extern gmp_version alias "__gmp_version" as const zstring const ptr
#endif

'initialization
declare sub mpz_init alias "__gmpz_init"(byval as mpz_ptr)
declare sub mpz_init2 alias "__gmpz_init2"(byval as mpz_ptr, byval as mp_bitcnt_t)
declare sub mpz_clear alias "__gmpz_clear"(byval as mpz_ptr)

'assignment
declare sub mpz_set alias "__gmpz_set"(byval as mpz_ptr, byval as mpz_srcptr)
declare sub mpz_set_ui alias "__gmpz_set_ui"(byval as mpz_ptr, byval as culong)
declare sub mpz_set_si alias "__gmpz_set_si"(byval as mpz_ptr, byval as clong)
declare function mpz_set_str alias "__gmpz_set_str"(byval as mpz_ptr, byval as const zstring ptr, byval as long) as long

'conversion
declare function mpz_get_d alias "__gmpz_get_d"(byval as mpz_srcptr) as double
declare function mpz_get_d_2exp alias "__gmpz_get_d_2exp"(byval as clong ptr, byval as mpz_srcptr) as double
declare function mpz_get_str alias "__gmpz_get_str"(byval as zstring ptr, byval as long, byval as mpz_srcptr) as zstring ptr

'arithmetic
declare sub mpz_add alias "__gmpz_add"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)
declare sub mpz_add_ui alias "__gmpz_add_ui"(byval as mpz_ptr, byval as mpz_srcptr, byval as culong)

declare sub mpz_sub alias "__gmpz_sub"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)

declare sub mpz_mul alias "__gmpz_mul"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)
declare sub mpz_mul_si alias "__gmpz_mul_si"(byval as mpz_ptr, byval as mpz_srcptr, byval as clong)
declare sub mpz_mul_ui alias "__gmpz_mul_ui"(byval as mpz_ptr, byval as mpz_srcptr, byval as culong)

declare sub mpz_addmul alias "__gmpz_addmul"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)
declare sub mpz_addmul_ui alias "__gmpz_addmul_ui"(byval as mpz_ptr, byval as mpz_srcptr, byval as culong)

declare sub mpz_submul alias "__gmpz_submul"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)
declare sub mpz_submul_ui alias "__gmpz_submul_ui"(byval as mpz_ptr, byval as mpz_srcptr, byval as culong)

declare sub mpz_mul_2exp alias "__gmpz_mul_2exp"(byval as mpz_ptr, byval as mpz_srcptr, byval as mp_bitcnt_t)

declare sub mpz_neg alias "__gmpz_neg"(byval as mpz_ptr, byval as mpz_srcptr)
declare sub mpz_abs alias "__gmpz_abs"(byval as mpz_ptr, byval as mpz_srcptr)

'division
declare sub mpz_fdiv_q alias "__gmpz_fdiv_q"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)

declare sub mpz_tdiv_q alias "__gmpz_tdiv_q"(byval as mpz_ptr, byval as mpz_srcptr, byval as mpz_srcptr)
declare function mpz_tdiv_q_ui alias "__gmpz_tdiv_q_ui"(byval as mpz_ptr, byval as mpz_srcptr, byval as culong) as culong
declare sub mpz_tdiv_q_2exp alias "__gmpz_tdiv_q_2exp"(byval as mpz_ptr, byval as mpz_srcptr, byval as mp_bitcnt_t)

'exponentiation
declare sub mpz_ui_pow_ui alias "__gmpz_ui_pow_ui"(byval as mpz_ptr, byval as culong, byval as culong)

'fast testing for negative, zero, and positive:
#define mpz_sgn(z) iif((z)->_mp_size < 0, -1, -((z)->_mp_size > 0))
#define mpz_odd_p(z) ((-((z)->_mp_size <> 0)) and clng((z)->_mp_d[0]))

'size
declare function mpz_sizeinbase alias "__gmpz_sizeinbase"(byval as mpz_srcptr, byval as long) as uinteger
declare function mpz_size alias "__gmpz_size"(byval as mpz_srcptr) as uinteger

end extern
