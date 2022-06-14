#define _CRT_SECURE_NO_WARNINGS // ms spam

#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "crude-xa.h"

// shorten types and defs

#define D XA_DIG
#define V XA_VAL

#define D_BYTES   XA_DIG_BYTES
#define D_BITS    XA_DIG_BITS
#define D_DEC_VAL XA_DIG_DEC_VAL
#define D_DEC_POW XA_DIG_DEC_POW
#define D_MASK(u) XA_DIG_MASK(u) 

#ifdef XA_VAL_LEAKS
#define V_LEAKS
#endif

#ifdef XA_AUTO_TEST

#undef xa_load
#undef xa_add
#undef xa_sub
#undef xa_mul

XA_VAL* xa_load(long double f);
XA_VAL* xa_add(const XA_VAL* a, const XA_VAL* b);
XA_VAL* xa_sub(const XA_VAL* a, const XA_VAL* b);
XA_VAL* xa_mul(const XA_VAL* a, const XA_VAL* b);

static void check(const char *desc, long double xa, long double ld)
{
	if (fabs(xa-ld) > fabs(ld) * pow(2.0, -52))
	{
		printf("%s\n",desc);
		// exit(0);
	}
}

V* xa_load_check(long double f)
{
	V* v = xa_load(f);
	long double vf = xa_extr(v);
	assert(f == vf);
	return v;
}

V* xa_add_check(const V* a, const V* b)
{
	long double fa = xa_extr(a);
	long double fb = xa_extr(b);

	V* va = xa_load(fa);
	V* vb = xa_load(fb);

	long double fc = fa+fb;
	V* vc = xa_add(va,vb);

	long double vcf = xa_extr(vc);
	//assert(vcf == fc);
	check("ADD",vcf, fc);

	V* v = xa_add(a, b);

	xa_free(va);
	xa_free(vb);
	xa_free(vc);

	return v;
}

V* xa_sub_check(const V* a, const V* b)
{
	long double fa = xa_extr(a);
	long double fb = xa_extr(b);

	V* va = xa_load(fa);
	V* vb = xa_load(fb);

	long double fc = fa-fb;
	V* vc = xa_sub(va,vb);

	long double vcf = xa_extr(vc);
	//assert(vcf == fc);
	check("SUB",vcf, fc);

	V* v = xa_sub(a, b);

	xa_free(va);
	xa_free(vb);
	xa_free(vc);

	return v;
}

V* xa_mul_check(const V* a, const V* b)
{
	long double fa = xa_extr(a);
	long double fb = xa_extr(b);

	V* va = xa_load(fa);
	V* vb = xa_load(fb);

	long double fc = fa*fb;
	V* vc = xa_mul(va,vb);

	long double vcf = xa_extr(vc);
	//assert(vcf == fc);
	check("MUL",vcf, fc);

	V* v = xa_mul(a, b);

	xa_free(va);
	xa_free(vb);
	xa_free(vc);

	return v;
}

#endif

// TO CONSIDER:
// refactor V::data from BigEndian to LittleEndian
// that will make possible to ignore need for V reallocations:
// - carry_overflows in xa_abs_add, but with inital alloc bigger by 1 digit!
// - borrow underflows in xa_sub_abs
// - reallocs in xa_mul

// TO CONSIDER:
// allow multi-chunk values, for results of add/sub with arguments
// having mantissas far disjoint like xa_add(0x123p+1000,0x456p-1000)
// so instead of >2000 bits of mantissa we'd have 2 short ones (~32 bits each).
// any subsequent operations must take care of merging chunks 
// probably (re)splitting wouldn't be worth such trouble
// alternatively instead of multi chunk we could wrap such values
// as tree of operators with referenced arguments

/*
typedef uint16_t D;
static const int D_BYTES = sizeof(D);
static const int D_BITS = D_BYTES * 8;
static const int D_DEC_VAL = 10000; // greatest 10^N smaller than 2^D_BITS
static const int D_DEC_POW = 4;	  // log10(D_DEC_VAL)
*/

/*
typedef uint8_t D;
static const int D_BYTES = sizeof(D);
static const int D_BITS = D_BYTES * 8;
static const int D_DEC_VAL = 100; // greatest 10^N smaller than 2^D_BITS
static const int D_DEC_POW = 2;	  // log10(D_DEC_VAL)
*/

typedef int _bool;
#define _false 0
#define _true 1

#ifdef V_LEAKS
static int alloc_bytes = 0;
static int num_allocs = 0;

struct Head
{
	Head* prev;
	Head* next;
	int id;
};

static const int head_bytes = ((sizeof(Head) -1) | 15) +1;

static Head* leak_head = 0;
static Head* leak_tail = 0;

static int break_id = 0;

void xa_break(int id)
{
	break_id = id;
}

int xa_leaks(int* bytes)
{
	if (bytes)
		*bytes = alloc_bytes;
	return num_allocs;
}

#endif

static int xa_pool_size = 0;
static V** xa_pool;

void xa_pool_alloc(int s)
{
	for (int i = s; i < xa_pool_size; i++)
	{
		V* v = xa_pool[i];
		while (v)
		{
			V* n = *(V**)v;
			free(v);
			v = n;
		}
	}

	xa_pool = (V**)realloc(xa_pool,sizeof(V*)*s);

	if (s > xa_pool_size)
		memset(xa_pool + xa_pool_size, 0, sizeof(V*)*(s - xa_pool_size));

	xa_pool_size = s;
}

void xa_pool_free()
{
	for (int i = 0; i < xa_pool_size; i++)
	{
		V* v = xa_pool[i];
		while (v)
		{
			V* n = *(V**)v;
			free(v);
			v = n;
		}
	}

	free(xa_pool);
	xa_pool_size = 0;
	xa_pool = 0;
}

V* xa_alloc(int digs)
{
	int s = sizeof(V) + digs * D_BYTES - D_BYTES;

	#ifdef V_LEAKS
	alloc_bytes += s;
	num_allocs++;
	#endif

	#ifdef V_LEAKS
	static int id = 1;
	if (id == break_id)
	{
		printf("break alloc %d\n", break_id);
	}

	Head* q = (Head*)malloc(s+head_bytes);
	q->prev = 0;
	q->next = leak_head;
	if (leak_head)
		leak_head->prev = q;
	else
		leak_tail = q;
	leak_head = q;
	q->id = id++;

	V* v = (V*)((intptr_t)q + head_bytes);
	#else

	V* v;
	if (digs < xa_pool_size && xa_pool[digs])
	{
		v = xa_pool[digs];
		xa_pool[digs] = *(V**)v;
	}
	else
		v = (V*)malloc(s);
	#endif

	v->refs = 1;
	v->digs = digs;
	return v;
}

void xa_free(V* v)
{
	if (v)
    {
        v->refs--;
        if (!v->refs)
		{
			#ifdef V_LEAKS
			int s = sizeof(V) + v->digs * D_BYTES - D_BYTES;
			alloc_bytes -= s;
			num_allocs--;
			Head* q = (Head*)((intptr_t)v - head_bytes);
			if (q->prev)
				q->prev->next = q->next;
			else
				leak_head = q->next;
			if (q->next)
				q->next->prev = q->prev;
			else
				leak_tail = q->prev;
			free(q);
			#else
			int digs = v->digs;
			if (digs < xa_pool_size)
			{
				*(V**)v = xa_pool[digs];
				xa_pool[digs] = v;
			}
			else
				free(v);
			#endif
		}
    }
}

void xa_grab(V* v)
{
	if (v)
		v->refs++;
}

int xa_extr_dec(const V* v, char **str)
{
	if (!v)
	{
		if (str)
		{
			char *out = (char *)malloc(2);
			out[0] = '0';
			out[1] = 0;
			*str = out;
		}

		return 1;
	}

	// max string length is:
	// numof quotient digits * log10(max_dig_val+1) + numof fraction digits * dig_bits
	// possibly trimmed to given precision

	//////////////////////////////////////////////////////////////////
	// prepare temp data buffer of length max(quotient_bits,fraction_bits)

	int quots = v->quot >= 0 ? v->quot + 1 : 0;
	int fracs = (int)v->digs > v->quot ? (int)v->digs - v->quot - 1 : 0;

	int digs = quots > fracs ? quots : fracs;
	D* buf = (D*)malloc(digs*D_BYTES);

	int decs = 4;				  // sign,period, round up and terminator
	decs += quots * D_BITS / 3; // bits * log10(2)
	decs += fracs * D_BITS;	  // wow, that can be big!

	char *out = (char *)malloc(decs);
	int dec = 0;

	//////////////////////////////////////////////////////////////////
	// quotient part: (less to more significant dec digits)

	// todo:
	// do not insert lo significant zeros into out, just count'em!
	// so if we have no fraction, we can replace zeros with exponent

	if (quots == 1)
	{
		//
		D rem = v->data[0];
		while (rem)
		{
			D d = rem % (D)10;
			rem /= (D)10;
			out[dec++] = '0' + (char)d;
		}
	}
	else if (quots > 0)
	{
		// first division pass:
		// v->data / 10^N -> buf
		D rem = 0;
		_bool scan = _true;
		int from = 0;
		int to = quots < (int)v->digs ? quots : (int)v->digs;
		
		// TODO: KEEP TRACK on last non zero buf[i], use it as 'to' in next iter
		for (int i = from; i < to; i++)
		{
			uint64_t x = ((uint64_t)rem << D_BITS) + v->data[i];
			D y = (D)(x / D_DEC_VAL);
			buf[i] = y;
			rem = (D)(x % D_DEC_VAL);

			if (scan && y)
			{
				scan = _false;
				from = i;
			}
		}

		// continue with out of data range
		for (int i = to; i<quots; i++)
		{
			if (!rem)
			{
				memset(buf + i, 0, (quots-i) * D_BYTES);
				break;
			}

			uint64_t x = (uint64_t)rem << D_BITS;
			D y = (D)(x / D_DEC_VAL);
			buf[i] = y;
			rem = (D)(x % D_DEC_VAL);

			if (scan && y)
			{
				scan = _false;
				from = i;
			}
		}

		if (scan)
		{
			// final dump, exclude leading zeros
			while (rem)
			{
				D d = rem % (D)10;
				rem /= (D)10;
				out[dec++] = '0' + (char)d;
			}
		}
		else
		{
			// full dump
			for (int i = 0; i < D_DEC_POW; i++)
			{
				D d = rem % (D)10;
				rem /= (D)10;
				out[dec++] = '0' + (char)d;
			}
		}

		while (!scan)
		{
			// repeat divisions buf / 10^N -> buf
			// until result > 0

			rem = 0;
			scan = _true;

			for (int i = from; i < quots; i++)
			{
				uint64_t x = ((uint64_t)rem << D_BITS) + buf[i];
				D y = (D)(x / D_DEC_VAL);
				buf[i] = y;
				rem = (D)(x % D_DEC_VAL);

				if (scan && y)
				{
					scan = _false;
					from = i;
				}
			}

			if (scan)
			{
				// final dump, exclude leading zeros
				while (rem)
				{
					D d = rem % (D)10;
					rem /= (D)10;
					out[dec++] = '0' + (char)d;
				}
			}
			else
			{
				// full dump
				for (int i = 0; i < D_DEC_POW; i++)
				{
					D d = rem % (D)10;
					rem /= (D)10;
					out[dec++] = '0' + (char)d;
				}
			}
		}
	}

	if (v->sign)
		out[dec++] = '-';

	// strrev(out);
	for (int i = 0; i < dec / 2; i++)
	{
		int j = dec - 1 - i;
		char s = out[i];
		out[i] = out[j];
		out[j] = s;
	}

	if (!dec)
		out[dec++] = '0';

	//////////////////////////////////////////////////////////////////
	// fractional part: (more to less significant dec digits)

	static const uint64_t mask = (~(uint64_t)0) >> (64 - D_BITS);

	if (fracs > 0)
	{
		out[dec++] = '.';

		// first pass multiplication:
		// v->data * 10^N -> buf
		int from = fracs-1;//v->digs - 1;
		int to = fracs > (int)v->digs ? fracs - (int)v->digs : 0;
		int ofs = v->digs - fracs;
		_bool scan = _true;
		uint64_t mul = 0;
		// TODO: KEEP TRACK on last non zero buf[i], use it as 'to' in next iter
		for (int i = from; i >= to; i--)
		{
			mul += (uint64_t)v->data[i+ofs] * D_DEC_VAL;
			buf[i] = D_MASK(mul); //(D)(mul & mask);
			mul >>= D_BITS;

			if (scan && buf[i])
			{
				scan = _false;
				from = i;
			}
		}

		// continue over out of data range
		int i = to-1;
		if (i>=0)
		{
			if (mul)
			{
				buf[i] = D_MASK(mul); //(D)(mul & mask);
				mul >>= D_BITS;
				if (scan && buf[i])
				{
					scan = _false;
					from = i;
				}

				i--;
			}

			if (i>=0)
				memset(buf,0,(i+1) * D_BYTES);
		}

		D rem = (D)mul;

		if (scan)
		{
			// final dump (reversed), exclude trailing zeros
			int z = 0;
			for (int i = 0; i < D_DEC_POW; i++)
			{
				D d = rem % (D)10;
				rem /= (D)10;

				if (i == z && !d)
				{
					z++;
					continue;
				}

				out[dec + D_DEC_POW - 1 - i] = '0' + (char)d;
			}
			dec += D_DEC_POW - z;
		}
		else
		{
			// full dump (reversed)
			for (int i = 0; i < D_DEC_POW; i++)
			{
				D d = rem % (D)10;
				rem /= (D)10;
				out[dec + D_DEC_POW - 1 - i] = '0' + (char)d;
			}
			dec += D_DEC_POW;
		}

		while (!scan)
		{
			// repeat multiplications buf * 10^N -> buf
			// until result contains fraction != 0
			scan = _true;
			mul = 0;
			// TODO: KEEP TRACK on last non zero buf[i], use it as 'to' in next iter
			for (int i = from; i >= 0; i--)
			{
				mul += (uint64_t)buf[i] * D_DEC_VAL;
				buf[i] = D_MASK(mul); //(D)(mul & mask);
				mul >>= D_BITS;

				if (scan && buf[i])
				{
					scan = _false;
					from = i;
				}
			}

			D rem = (D)mul;

			if (scan)
			{
				// final dump (reversed), exclude trailing zeros
				int z = 0;
				for (int i = 0; i < D_DEC_POW; i++)
				{
					D d = rem % (D)10;
					rem /= (D)10;

					if (i == z && !d)
					{
						z++;
						continue;
					}

					out[dec + D_DEC_POW - 1 - i] = '0' + (char)d;
				}
				dec += D_DEC_POW - z;
			}
			else
			{
				// full dump (reversed)
				for (int i = 0; i < D_DEC_POW; i++)
				{
					D d = rem % (D)10;
					rem /= (D)10;
					out[dec + D_DEC_POW - 1 - i] = '0' + (char)d;
				}
				dec += D_DEC_POW;
			}
		}
	}

	free(buf);

	if (str)
	{
		*str = out;
		out[dec] = 0;
	}
	else
		free(out);

	return dec;
}

int xa_extr_hex(const V* v, char **str)
{
	if (!v)
	{
		if (str)
		{
			*str = (char*)malloc(7);
			strcpy(*str, "0x0p+0");
		}

		return 6;
	}

	char* out = 0;

	if (str)
	{
		int len = 17; // -0x0p+2000000000#
		len += (v->digs * D_BITS + 3) / 4;
		out = (char*)malloc(len);
	}

	int len = 0;
	if (out)
	{
		if (v->sign)
			out[len++] = '-';
		out[len++] = '0';
		out[len++] = 'x';
	}
	else
		len += v->sign ? 3 : 2;

	int i = 0;
	int right = D_BITS; // how many bits are still available in data[i]
	D d = v->data[i];
	while (!(d & (1 << (right - 1))))
		right--;

	int e = v->quot * D_BITS + right - 4;

	int period = _false;
	uint64_t acc = 0;
	int left = 64; // how many bits left free in accumulator
	while (1)
	{
		if (left < right)
		{
			// we're full
			acc <<= left;
			acc |= v->data[i] >> (right - left);
			right -= left;

			// check tail condition
			_bool tail = _false;

			if (i == v->digs - 1)
			{
				if ((v->data[i] & ((1 << right) - 1)) == 0)
					tail = _true;
			}

			if (tail)
			{
				// dump until acc is not 0
				if (out)
				{
					while (acc)
					{
						int d = (int)(acc >> 60);
						out[len++] = d <= 9 ? '0' + d : 'a' - 10 + d;
						acc <<= 4;
						if (!period && acc)
						{
							period = len;
							out[len++] = '.';
						}
					}
				}
				else
				{
					while (acc) // possible bit-scan inoutuction
					{
						len++;
						acc <<= 4;
						if (!period && acc)
							period = len;
					}
				}

				// we're done
				break;
			}
			else
			{
				// dump all 16 digits
				if (out)
				{
					for (int h = 0; h < 16; h++)
					{
						int d = (int)(acc >> 60);
						out[len++] = d <= 9 ? '0' + d : 'a' - 10 + d;
						if (!period)
						{
							period = len;
							out[len++] = '.';
						}
						acc <<= 4;
					}
				}
				else
				{
					if (period)
						len += 16;
					else
					{
						period = len + 1;
						len += 17;
					}
				}
			}

			acc = 0;
			left = 64;
		}
		else
		{
			acc <<= right;
			acc |= v->data[i];
			left -= right;

			i++;

			if (i == v->digs)
			{
				// final, partial dump
				acc <<= left;
				if (out)
				{
					while (acc)
					{
						int d = (int)(acc >> 60);
						out[len++] = d <= 9 ? '0' + d : 'a' - 10 + d;
						acc <<= 4;
						if (!period && acc)
						{
							period = len;
							out[len++] = '.';
						}
					}
				}
				else
				{
					while (acc) // possible bit-scan inoutuction
					{
						len++;
						acc <<= 4;
						if (!period && acc)
						{
							period = len + 1;
							len++;
						}
					}
				}
				break;
			}

			right = D_BITS;
		}
	}

	char exp_buf[16];
	int exp_len = sprintf(exp_buf, "p%+d", e);

	if (out)
		strcpy(out + len, exp_buf);
	len += exp_len;

	if (str)
		*str = out;

	return len;
}


long double xa_extr(const V* v)
{
	// TODO:
	// rounding modes
	// #include <fenv.h>
	// fegetround() ->
	//  FE_DOWNWARD   Round downward.
	//  FE_TONEAREST  Round to nearest.  <- implemented
	//  FE_TOWARDZERO Round toward zero. <- no rounding
	//  FE_UPWARD     Round upward.

	if (!v)
		return 0;

	// scan hi zeros
	D d = v->data[0];
	D m = ((D)1) << (D_BITS - 1);
	int z = 0;
	// possible use of bit-scan instruction
	while (!(d & m))
	{
		m >>= 1;
		z++;
	}

	uint64_t u = 0;
	for (int i = 0; i < (int)v->digs; i++)
	{
		int shl = 64 - (D_BITS - z) - i * D_BITS;
		if (shl <= -(int)D_BITS)
			break;

		uint64_t x;
		if (shl >= 0)
			x = ((uint64_t)v->data[i]) << shl;
		else
			x = ((uint64_t)v->data[i]) >> -shl;

		u |= x;
	}

	int e = v->quot * D_BITS - z - (64 - D_BITS);

	// ROUNDING: to nearest, tie to even

	int rdig = (64 + z) / D_BITS;
	if (rdig < (int)v->digs)
	{
		// check 1 bit after least significant bit of destinatiom
		int rbit = D_BITS - 1 - (64 + z - rdig * D_BITS);

		if (v->data[rdig] & (((D)1) << rbit)) // if it is '1' we'll increment u
		{
			// hold on a sec!
			// do not increment if u&1 is '0' and rbit is the very last '1' bit in data
			if ((u & 1) == 0 && rdig == v->digs - 1 && (v->data[rdig] & ((1 << rbit) - 1)) == 0)
			{
				// don't round ties to odd
			}
			else
			{
				u++;	// round up
				if (!u) // if u overflows
				{
					// set u=1 and increment e by 64
					e += 64;
					u = 1;
				}
			}
		}
	}
	return v->sign ? -ldexp((long double)u, e) : ldexp((long double)u, e);
}

V* xa_load(long double v)
{
	if (!v)
	{
		// null is treated as exact 0.0
		return 0;
	}

	assert(D_BITS >= 4 && D_BITS <= 63 && "digit misconfiguration");
	assert(FLT_RADIX == 2 && "other representations would not be exact");
	assert(LDBL_MANT_DIG <= 64 && "currently mantissa extraction is limited to 64 bits");

	// static const long double hi = ldexp(1.0, 64);
	static const long double hi = 0x1.0p+64; // cuz ms

	int sign;

	int e;
	uint64_t u;
	if (v > 0)
	{
		u = (uint64_t)(frexp(v, &e) * hi);
		sign = 0;
	}
	else
	{
		u = (uint64_t)(frexp(-v, &e) * hi);
		sign = 1;
	}

	// shift out lo significant zeros
	// possible use of bit-scan instruction
	int insignificant = 0;
	while (!((u >> insignificant) & 1))
		insignificant++;

	// align
	int shl = 0;
	if (e > 0)
	{
		int f = e / D_BITS;
		shl = e - f * D_BITS;
		e = f;
	}
	else if (e < 0)
	{
		int f = (e - D_BITS + 1) / D_BITS;
		shl = e - f * D_BITS;
		e = f;
	}

	int digs;
	if (shl)
	{
		digs = 1 + (64 - shl - insignificant + D_BITS - 1) / D_BITS;
	}
	else
	{
		e--;
		shl += D_BITS;
		digs = (64 - insignificant + D_BITS - 1) / D_BITS;
	}

	// allows using bigger dig than dig_bits
	static const uint64_t mask = (~(uint64_t)0) >> (64 - D_BITS);

	V *ret = xa_alloc(digs);
	ret->sign = sign;
	ret->quot = e;

	for (int dig = 0; dig < digs; dig++)
	{
		int shr = 64 - shl - dig * D_BITS;
		if (shr >= 0)
			ret->data[dig] = D_MASK((u >> (shr & 63)) & mask); //(D)((u >> (shr & 63)) & mask);
		else
			ret->data[dig] = D_MASK((u << ((-shr) & 63)) & mask); //(D)((u << ((-shr) & 63)) & mask);
	}

	return ret;
}

V* xa_cpy(const V* v)
{
	if (!v)
		return 0;
	V *ret = xa_alloc(v->digs);
	memcpy(ret->data, v->data, v->digs * D_BYTES);
	ret->sign = v->sign;
	ret->quot = v->quot;

	return ret;
}

/*
int xa_cmp(const V* v)
{
	// compares with 0
	if (!v)
		return 0;
	if (v->sign)
		return -1;
	return 1;
}
*/

int xa_cmp(const V *a, const V *b)
{
	if (a == b)
		return 0;

	if (b == 0)
		return a->sign ? -1 : 1;

	if (a == 0)
		return b->sign ? 1 : -1;

	if (a->sign)
	{
		if (!b->sign)
			return -1;
	}
	else
	{
		if (b->sign)
			return 1;
	}

	// signs are equal
	// compare exponents
	if (a->quot < b->quot)
		return a->sign ? 1 : -1;

	if (a->quot > b->quot)
		return a->sign ? -1 : 1;

	// signs and exponents are equal
	// compare data
	const D *a_d = a->data;
	const D *b_d = b->data;
	int digs = a->digs < b->digs ? a->digs : b->digs;

	for (int dig = 0; dig < digs; dig++)
	{
		if (a->data[dig] < b->data[dig])
			return a->sign ? 1 : -1;

		if (a->data[dig] > b->data[dig])
			return a->sign ? -1 : 1;
	}

	// signs exponents and common data range are same
	// compare dig lengths!!!

	// longer one means bigger abs val
	if (a->digs == b->digs)
		return 0;

	if (a->sign)
		return a->digs < b->digs ? 1 : -1;

	return a->digs < b->digs ? -1 : 1;
}

const V* xa_min(const V *a, const V *b)
{
	return xa_cmp(a, b) <= 0 ? a : b; // reference maybe?
}

const V* xa_max(const V *a, const V *b)
{
	return xa_cmp(a, b) >= 0 ? a : b; // reference maybe?
}

V* xa_abs(const V* v)
{
	if (!v)
		return 0;
	V *r = xa_cpy(v); // reference with abs bit ?
	r->sign = 0;
	return r;
}

V* xa_neg(const V* v)
{
	if (!v)
		return 0;
	V *r = xa_cpy(v); // reference with neg bit ?
	r->sign ^= 1;
	return r;
}

static int carry_overflows = 0;
static V* xa_add_abs(const V *a, const V *b)
{
	// internal, abs(a)+abs(b)
	// used by xa_add & xa_sub
	// always returns positive sign
	// cases where a==0 || b==0 are not handled here

	int head = a->quot > b->quot ? a->quot : b->quot;
	int tail_a = a->quot - a->digs;
	int tail_b = b->quot - b->digs;
	int digs = head - (tail_a < tail_b ? tail_a : tail_b);

	if (tail_a >= b->quot || tail_b >= a->quot)
	{
		// if a and b are not overlapping then
		// we can simply copy each source to its corresponding position in destination
		V* v = xa_alloc(digs);
		v->sign = 0;
		v->quot = head;

		int a_ofs = head - a->quot;
		memcpy(v->data + a_ofs, a->data, a->digs * D_BYTES);

		int b_ofs = head - b->quot;
		memcpy(v->data + b_ofs, b->data, b->digs * D_BYTES);

		// and fill the gap if needed
		if (tail_a > b->quot)
			memset(v->data + a->digs + a_ofs, 0, (tail_a - b->quot) * D_BYTES);
		else if (tail_b > a->quot)
			memset(v->data + b->digs + b_ofs, 0, (tail_b - a->quot) * D_BYTES);

		return v;
	}

	// static const uint64_t mask = (((uint64_t)1)<<D_BITS) - 1;
	static const uint64_t mask = (~(uint64_t)0) >> (64 - D_BITS);

	uint64_t sum = 0;

	int adig = a->digs - 1;
	int bdig = b->digs - 1;
	int zero = 1;

	if (tail_a == tail_b)
	{
		// perform addition from tail
		// to see how many zeros we can get there.
		// we keep result of lowest significant non-zero digit and carry flag
		// se we will be able to continue after allocating destination
		while (adig >= 0 && bdig >= 0)
		{
			sum += a->data[adig--];
			sum += b->data[bdig--];

			if (!D_MASK(sum)/*(sum & mask)*/)
			{
				digs--;
			}
			else
			{
				zero = 0;
				break;
			}

			sum >>= D_BITS;
		}

		if (zero)
		{
			while (adig >= 0)
			{
				sum += a->data[adig--];

				if (!D_MASK(sum)/*(sum & mask)*/)
				{
					digs--;
				}
				else
				{
					zero = 0;
					break;
				}

				sum >>= D_BITS;
			}

			while (bdig >= 0)
			{
				sum += b->data[bdig--];

				if (!D_MASK(sum)/*(sum & mask)*/)
				{
					digs--;
				}
				else
				{
					zero = 0;
					break;
				}

				sum >>= D_BITS;
			}
		}

		/*
		// todo: fixme to avoid 0 digs allocs & reallocs to 1
		if (!digs)
		{
			head++;
			digs = 1;
		}
		*/
	}

	// we alloc V optimistically
	// assuming there will be no need for extra digit at head
	// but if that happens we will have to re-alloc

	// now we can alloc V
	V* v = xa_alloc(digs);
	v->sign = 0;
	v->quot = head;

	int dig = digs - 1;

	if (tail_a < tail_b)
	{
		// note: non overlapping case already handled
		if (a->quot <= b->quot)
		{
			// [  [aaaaa]
			// [bbbbb]

			// memcpy aa]
			int num = tail_b - tail_a;
			memcpy(v->data + digs - num, a->data + a->digs - num, num * D_BYTES);
			adig -= num;
			dig -= num;
			// sum aa+bb
			num = a->quot - tail_b;
			for (int i = 0; i < num; i++)
			{
				sum += a->data[adig--];
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
			// continue with [bb (possibly 0 len)
			// TODO: once carry disapear we can finish with memcpy b
			num = b->quot - a->quot;
			for (int i = 0; i < num; i++)
			{
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
		}
		else
		{
			// [aaaaaaaa]
			//   [bbbb]

			// memcpy a]
			int num = tail_b - tail_a;
			memcpy(v->data + digs - num, a->data + a->digs - num, num * D_BYTES);
			adig -= num;
			dig -= num;
			// sum aaaa+bbbb
			num = b->quot - tail_b;
			for (int i = 0; i < num; i++)
			{
				sum += a->data[adig--];
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
			// continue with [a
			num = a->quot - b->quot;
			for (int i = 0; i < num; i++)
			{
				if (!sum)
				{
					// memcpy !!! 2
					// break
					memcpy(v->data, a->data, (num-i) * D_BYTES);
					break;
				}
				sum += a->data[adig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
		}
	}
	else if (tail_a > tail_b)
	{
		// note: non overlapping case already handled
		if (b->quot <= a->quot)
		{
			// [aaaaa]
			// [  [bbbbb]

			// memcpy bb]
			int num = tail_a - tail_b;
			memcpy(v->data + digs - num, b->data + b->digs - num, num * D_BYTES);
			bdig -= num;
			dig -= num;
			// sum aa+bb
			num = b->quot - tail_a;
			for (int i = 0; i < num; i++)
			{
				sum += a->data[adig--];
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
			// continue wth [aa (possibly 0 len)
			num = a->quot - b->quot;
			for (int i = 0; i < num; i++)
			{
				if (!sum)
				{
					// memcpy !!! 1
					// break
					memcpy(v->data, a->data, (num-i) * D_BYTES);
					break;
				}
				sum += a->data[adig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
		}
		else
		{
			//   [aaaa]
			// [bbbbbbbb]

			// memcpy b]
			int num = tail_a - tail_b;
			memcpy(v->data + digs - num, b->data + b->digs - num, num * D_BYTES);
			bdig -= num;
			dig -= num;
			// sum aaaa+bbbb
			num = a->quot - tail_a;
			for (int i = 0; i < num; i++)
			{
				sum += a->data[adig--];
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
			// continue with [b
			num = b->quot - a->quot;
			for (int i = 0; i < num; i++)
			{
				if (!sum)
				{
					// memcpy !!! 3
					// break
					memcpy(v->data, b->data, (num-i) * D_BYTES);
					break;
				}
				sum += b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
			}
		}
	}
	else
	{
		// we've already started this case
		// store current values and continue

		if (!zero && dig>=0)
		{
			v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
			sum >>= D_BITS;
		}

		while (adig >= 0 && bdig >= 0)
		{
			sum += a->data[adig--];
			sum += b->data[bdig--];
			v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
			sum >>= D_BITS;
		}

		while (adig >= 0)
		{
			sum += a->data[adig--];
			v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
			sum >>= D_BITS;
		}

		while (bdig >= 0)
		{
			sum += b->data[bdig--];
			v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
			sum >>= D_BITS;
		}
	}

	if (sum) // ugh! carry overflow
	{
		carry_overflows++;

		digs++;
		head++;

		V *r = xa_alloc(digs);
		r->sign = 0;
		r->quot = head;

		memcpy(r->data + 1, v->data, (digs - 1) * D_BYTES);
		r->data[0] = D_MASK(sum); //sum & mask;

		xa_free(v);
		v = r;
	}

	return v;
}

static int borrow_underflows = 0;
static V* xa_sub_abs(const V *a, const V *b)
{
	if (a == b)
		return 0;

	// internal, abs(a)-abs(b)
	// used by xa_add & xa_sub
	// returns negative sign if abs(a) < abs(b)
	// return 0 if abs(a)==abs(b)

	// IF quots are same
	// 	march from head to tail until digs are equal
	// 	keep going even if reached one tail but another digit is 0
	// 	if reached both tails, return 0 !!!
	// 	remember which number (by comparing last digits) was greater
	// 	remember num of digits excluded
	// ELSE
	//  bigger number is the one with bigger exponent
	//  no skipped digits

	int head = a->quot > b->quot ? a->quot : b->quot;
	int tail_a = a->quot - a->digs;
	int tail_b = b->quot - b->digs;
	int digs = head - (tail_a < tail_b ? tail_a : tail_b);

	_bool swapped = _false;
	if (a->quot == b->quot)
	{
		int dig = 0;
		int num = a->digs < b->digs ? a->digs : b->digs;
		for (; dig < num; dig++)
		{
			if (a->data[dig] > b->data[dig])
				break;
			if (a->data[dig] < b->data[dig])
			{
				const V *s = a;
				a = b;
				b = s;

				swapped = _true;
				tail_a = a->quot - a->digs;
				tail_b = b->quot - b->digs;
				break;
			}
		}

		if (dig == digs)
		{
			// a & b are equal!
			return 0;
		}

		if (dig == num)
		{
			// entire shorter arg was identical to the longer
			// clip additional 0s and return remaining part of longer arg
			if (a->digs < b->digs)
			{
				const V *s = a;
				a = b;
				b = s;
				swapped = _true;
			}

			// while a (longer) contins 0s
			while (dig < digs && !a->data[dig])
				dig++;

			// now we should simply alloc V for ramaining part of A and do memcpy
			V* v = xa_alloc(digs - dig);
			v->sign = swapped ? 1 : 0;
			v->quot = head - dig;
			memcpy(v->data, a->data + dig, v->digs * D_BYTES);

			return v;
		}

		// trim
		head -= dig;
		digs -= dig;
	}
	else if (a->quot < b->quot)
	{
		const V *s = a;
		a = b;
		b = s;
		swapped = _true;
		tail_a = a->quot - a->digs;
		tail_b = b->quot - b->digs;
	}

	// if tails are same
	// march from tail to head until digs are equal (sub==0)
	// remember num of digits excluded
	int adig = a->digs - 1;
	int bdig = b->digs - 1;

	if (tail_a == tail_b)
	{
		int dig = 0;
		int num = a->digs < b->digs ? a->digs : b->digs;

		for (; dig < num; dig++)
		{
			if (a->data[adig] != b->data[bdig])
				break;
			adig--;
			bdig--;
		}

		if (dig == num)
		{
			// continue with longer till contains 0s
			while (dig < digs && !a->data[adig])
			{
				adig--;
				dig++;
			}
		}

		// trim
		tail_a += dig;
		tail_b += dig;
		digs -= dig;
	}

	int dig = digs - 1;

	// ALLOC V with digs = head - tail - head_excluded - tail_excluded
	// IF bigger was value B, swap A&B and store negative sign in V
	V* v = xa_alloc(digs);
	v->sign = swapped ? 1 : 0;
	v->quot = head;

	static const uint64_t mask = (~(uint64_t)0) >> (64 - D_BITS);
	
	//static const D dmsk = (D)mask;
	static const D dmsk = (D)((~(uint64_t)0) >> (64 - D_BITS)); // cuz ms

	// run borrowing algorithm on trimmed range (A is always bigger here but not neccessarily longer)
	uint64_t sum = 0; // 0 - not borrowed, 1 - borrowed

	if (tail_a >= b->quot && bdig>=0)
	{
		// not trimmed
		// [aa] [bbbb]

		v->data[dig--] = ((D)0 - b->data[bdig--]) & dmsk; // first digit with 1<<bits borrow (use 0)
		int num = digs - (a->quot - b->quot) - 1;
		for (int i = 0; i < num; i++)
			v->data[dig--] = ((D)-1 - b->data[bdig--]) & dmsk; // next digits with (1<<bits)-1 borrow (use -1)

		// if D_BITS == D_BYTES*8
		// we could make memset here
		num = tail_a - b->quot;
		D fill = (D)mask;
		for (int i = 0; i < num; i++)
			v->data[dig--] = fill;

		v->data[dig--] = (a->data[adig--] - (D)1) & dmsk; // handle borrow, here a->data can't be 0

		num = a->digs - 1;
		if (num)
			memcpy(v->data, a->data, num * D_BYTES);
	}
	else if (a->quot == b->quot)
	{
		if (tail_a > tail_b)
		{
			// possibly head trimmed!!!!
			// [aaaa]
			// [bbbbbbbb]

			v->data[dig--] = ((D)0 - b->data[bdig--]) & dmsk;
			int num = tail_a - tail_b - 1;
			for (int i = 0; i < num; i++)
				v->data[dig--] = ((D)-1 - b->data[bdig--]) & dmsk;

			sum = (uint64_t)~0;
			while (dig >= 0)
			{
				sum += a->data[adig--];
				sum -= b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}
		}
		else
		{
			if (tail_b > tail_a)
			{
				// possibly head trimmed!!!!
				// [aaaaaaaa]
				// [bbbb]

				int num = tail_b - tail_a; // must be >0
				adig -= num;
				dig -= num;
				memcpy(v->data + dig + 1, a->data + adig + 1, num * D_BYTES);
			}

			// possibly head and / or tail trimmed!!!!
			// [aaaaaa]
			// [bbbbbb]

			while (dig >= 0)
			{
				sum += a->data[adig--];
				sum -= b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}
		}
	}
	else
	{
		if (tail_a > tail_b)
		{
			// not trimmed
			// [aaaa]
			//    [bbbbbbb]

			v->data[dig--] = ((D)0 - b->data[bdig--]) & dmsk;
			int num = tail_a - tail_b - 1;
			for (int i = 0; i < num; i++)
				v->data[dig--] = ((D)-1 - b->data[bdig--]) & dmsk;

			sum = ~(uint64_t)0;
			num = b->quot - tail_a;
			for (int i = 0; i < num; i++)
			{
				sum += a->data[adig--];
				sum -= b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}

			while (dig >= 0)
			{
				sum += a->data[adig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}
		}
		else
		{
			if (tail_a < tail_b)
			{
				// not trimmed
				// [aaaaaaa]
				//   [bbb]

				int num = tail_b - tail_a; // must be >0
				adig -= num;
				dig -= num;
				memcpy(v->data + dig + 1, a->data + adig + 1, num * D_BYTES);
			}

			// possibly tail trimmed!!!!
			// [aaaaaaaa]
			//     [bbbb]

			while (bdig >= 0)
			{
				sum += a->data[adig--];
				sum -= b->data[bdig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}

			while (dig >= 0)
			{
				sum += a->data[adig--];
				v->data[dig--] = D_MASK(sum); //(D)(sum & mask);
				sum >>= D_BITS;
				if (sum)
					sum = (uint64_t)~0;
			}
		}
	}

	// it may happen that due to borrowing one or more head digits has gone (equals 0)
	// in such case realloc V and trim zeros out
	int clip = 0;
	while (v->data[clip] == 0)
		clip++;
	if (clip)
	{
		borrow_underflows += clip;

		head-=clip;
		digs-=clip;

		V *r = xa_alloc(digs);
		r->sign = swapped ? 1 : 0;
		r->quot = head;

		memcpy(r->data, v->data + clip, digs * D_BYTES);

		xa_free(v);
		v = r;
	}

	return v;
}

V* xa_add(const V *a, const V *b)
{
	if (!a)
	{
		if (!b)
			return 0;
		return xa_cpy(b);
	}

	if (!b)
		return xa_cpy(a);

	if (a->sign == b->sign)
	{
		V *r = xa_add_abs(a, b);
		r->sign = a->sign;
		return r;
	}

	return b->sign ? xa_sub_abs(a, b) : xa_sub_abs(b, a);
}

V* xa_sub(const V *a, const V *b)
{
	if (!a)
	{
		if (!b)
			return 0;
		V* v = xa_cpy(b);
		v->sign ^= 1;
		return v;
	}


	if (!b)
		return xa_cpy(a);

	if (a->sign == b->sign)
	{
		V *r = xa_sub_abs(a, b);
		if (r)
			r->sign ^= a->sign;
		return r;
	}

	if (b->sign)
		return xa_add_abs(a, b);

	V *r = xa_add_abs(b, a);
	r->sign = 1;

	return r;
}

V* xa_mul(const V *a, const V *b)
{
	if (!a || !b)
		return 0;

	static const uint64_t mask = (~(uint64_t)0) >> (64 - D_BITS);

	int digs = a->digs + b->digs;
	int quot = a->quot + b->quot + 1;

	// optimize size by checking himul
	_bool predicted = _false;
	uint64_t himul = (uint64_t)a->data[0] * (uint64_t)b->data[0];
	if (himul <= mask*245/256)
	{
		/*
		digs--;
		quot--;
		*/
		predicted = _true;
	}

	// optimize size by checking lomul
	int tail_clip = 0;
	uint64_t lomul = (uint64_t)a->data[a->digs - 1] * (uint64_t)b->data[b->digs - 1];
	if (!(lomul & mask))
	{
		tail_clip = 1;
		digs--;
	}

	V* v = xa_alloc(digs);
	v->sign = a->sign ^ b->sign;
	v->quot = quot;

	// make b (nested loop) longer than a 
	if (a->digs > b->digs)
	{
		const V *s = a;
		a = b;
		b = s;
	}

	// first pass a*b[0]->v
	int dig = digs - 1;
	int adig = a->digs - 1;
	D aval = a->data[adig];
	uint64_t mul = 0;

	if (tail_clip)
	{
		mul = lomul >> D_BITS;
		for (int bdig = b->digs - 2; bdig >= 0; bdig--)
		{
			mul += (uint64_t)aval * (uint64_t)b->data[bdig];
			v->data[dig--] = D_MASK(mul); //(D)(mul & mask);
			mul >>= D_BITS;
		}
	}
	else
	{
		for (int bdig = b->digs - 1; bdig >= 0; bdig--)
		{
			mul += (uint64_t)aval * (uint64_t)b->data[bdig];
			v->data[dig--] = D_MASK(mul); //(D)(mul & mask);
			mul >>= D_BITS;
		}
	}
	v->data[dig] = (D)(mul); // must fit, no need to mask
	adig--;

	while (adig >= 0)
	{
		v->data[dig-1]=0; // uninitialized yet, we are going to access it on last bdig iter
		dig = digs + tail_clip - (a->digs - adig);
		D aval = a->data[adig];
		uint64_t mul = 0;
		for (int bdig = b->digs - 1; bdig >= 0; bdig--)
		{
			mul += (uint64_t)aval * (uint64_t)b->data[bdig] + v->data[dig];
			v->data[dig--] = D_MASK(mul); //(D)(mul & mask);
			mul >>= D_BITS;
		}
		v->data[dig] = (D)(mul + v->data[dig]); // must fit, no need to mask
		adig--;
	}

	if (!v->data[0])
	{
		digs--;
		V *r = xa_alloc(digs);
		r->quot = v->quot - 1;
		r->sign = v->sign;
		memcpy(r->data, v->data + 1, digs * D_BYTES);

		xa_free(v);
		v = r;
	}

	return v;
}

///////////////////////////////////////////////////////////////////////////////////
