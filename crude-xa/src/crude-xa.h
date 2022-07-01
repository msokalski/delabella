#ifndef CRUDE_XA_H
#define CRUDE_XA_H

#include <stdint.h>

// must be enabled for xa_leaks
// #define XA_VAL_LEAKS

// crude-xa internal tests, this is super slow,
// don't use unless you suspect a bug in crude-xa
//#define XA_AUTO_TEST

typedef uint32_t XA_DIG;
#define XA_DIG_BYTES ((int)sizeof(XA_DIG))
#define XA_DIG_BITS (XA_DIG_BYTES * 8)
#define XA_DIG_DEC_VAL 1000000000 // greatest 10^N smaller than 2^XA_DIG_BITS
#define XA_DIG_DEC_POW 9		   // log10(DIG_DEC_VAL) == N

// use this if (XA_DIG_BITS < DIG_BYTES * 8)
// #define XA_DIG_MASK(u) ((XA_DIG)((u) & mask)) 

// this one otherwise (faster)
#define XA_DIG_MASK(u) ((XA_DIG)(u))

typedef struct
{
	unsigned int sign : 1;  // 0 positive, 1 negative
	unsigned int digs : 30; // data size in digits
	int quot;				// first digit multiplier = 2^(XA_DIG_BITS*quot)
	int refs;               // grab counter
	XA_DIG data[1];			// big endian (digits' bytes are in machine order)
} XA_VAL;

#ifdef __cplusplus
extern "C" {
#endif

// LO LEVEL C API

void xa_pool_alloc(int s);
void xa_pool_free();
int xa_pool_stat();

XA_VAL* xa_alloc(int digs);
void xa_free(XA_VAL* v);
void xa_grab(XA_VAL* v);

#ifdef XA_VAL_LEAKS
int xa_leaks(int* bytes);
void xa_break(int id);
#endif

int xa_extr_dec(const XA_VAL* v, char** str);
int xa_extr_hex(const XA_VAL* v, char** str);
long double xa_extr(const XA_VAL* v);
XA_VAL* xa_load(long double v);
XA_VAL* xa_cpy(const XA_VAL* v);
int xa_cmp(const XA_VAL* a, const XA_VAL* b);
const XA_VAL* xa_min(const XA_VAL* a, const XA_VAL* b);
const XA_VAL* xa_max(const XA_VAL* a, const XA_VAL* b);
XA_VAL* xa_abs(const XA_VAL* v);
XA_VAL* xa_neg(const XA_VAL* v);

#ifdef XA_AUTO_TEST
// only these 4 are auto tested
// assuming xa_extr works fine
#define xa_load xa_load_check
#define xa_add xa_add_check
#define xa_sub xa_sub_check
#define xa_mul xa_mul_check
#endif

XA_VAL* xa_load(long double v);
XA_VAL* xa_add(const XA_VAL* a, const XA_VAL* b);
XA_VAL* xa_sub(const XA_VAL* a, const XA_VAL* b);
XA_VAL* xa_mul(const XA_VAL* a, const XA_VAL* b);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <iostream>

struct XA_REF
{
    ~XA_REF()
    {
        xa_free(v);
    }

	XA_REF()
    {
        v=0;
    }

	int cmp(const XA_REF& w) const
	{
		return xa_cmp(v, w.v);
	}

	operator long double () const
	{
		return xa_extr(v);
	}

	XA_REF(const long double& f)
    {
        v=xa_load(f);
    }

	XA_REF(const XA_REF& w)
    {
        xa_grab(w.v);
        v=w.v;
    }

	XA_REF& operator = (const XA_REF& w)
    {
        xa_grab(w.v);
		xa_free(v);
		v=w.v;
        return *this;
    }

	XA_REF& operator = (const long double& f)
    {
		xa_free(v);
        v=xa_load(f);
        return *this;
    }

	XA_REF(XA_REF&& w)
    {
        v = w.v;
        w.v = 0;
    }

	XA_REF& operator = (XA_REF&& w)
    {
		xa_free(v);
		v = w.v;
		w.v = 0;
		return *this;
	}

    bool operator == (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) == 0;
		xa_free(j);
		return r;
    }

    bool operator == (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) == 0;
    }

    bool operator != (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) != 0;
		xa_free(j);
		return r;
    }

    bool operator != (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) != 0;
    }

    bool operator < (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) < 0;
		xa_free(j);
		return r;
    }

    bool operator < (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) < 0;
    }

    bool operator > (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) > 0;
		xa_free(j);
		return r;
    }

    bool operator > (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) > 0;
    }

    bool operator <= (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) <= 0;
		xa_free(j);
		return r;
    }

    bool operator <= (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) <= 0;
    }

    bool operator >= (int i) const
    {
		XA_VAL* j = xa_load(i);
        bool r = xa_cmp(v,j) >= 0;
		xa_free(j);
		return r;
    }

    bool operator >= (const XA_REF& w) const
    {
        return xa_cmp(v,w.v) >= 0;
    }

	XA_REF operator + (const XA_REF& w) const
    {
		XA_REF r;
        r.v = xa_add(v,w.v);
        return r;
    }

	XA_REF& operator += (const XA_REF& w)
    {
		XA_VAL* r = xa_add(v,w.v);
        xa_free(v);
        v = r;
        return *this;
    }

	XA_REF operator - (const XA_REF& w) const
    {
		XA_REF r;
        r.v = xa_sub(v,w.v);
        return r;
    }

	XA_REF& operator -= (const XA_REF& w)
    {
		XA_VAL* r = xa_sub(v,w.v);
        xa_free(v);
        v = r;
        return *this;
    }    

	XA_REF operator * (const XA_REF& w) const
    {
		XA_REF r;
        r.v = xa_mul(v,w.v);
        return r;
    }

	XA_REF& operator *= (const XA_REF& w)
    {
		XA_VAL* r = xa_mul(v,w.v);
        xa_free(v);
        v = r;
        return *this;
    }
 
	XA_REF operator << (int i) const
    {
		if (v && i>0)
		{
			// no field names cuz ms
			XA_VAL m = {/*sign:*/0,/*digs:*/1,/*quot:*/i/(int)XA_DIG_BITS,/*refs:*/1,
			       /*data:*/{(XA_DIG)1<<(i%XA_DIG_BITS)}};
			XA_REF r;
			r.v = xa_mul(&m,v);
			return r;
		}
        return *this;
    }

	XA_REF& operator <<= (int i)
    {
		if (v && i>0)
		{
			// no field names cuz ms
			XA_VAL m = {/*sign:*/0,/*digs:*/1,/*quot:*/i/(int)XA_DIG_BITS,/*refs:*/1,
			       /*data:*/{(XA_DIG)1<<(i%XA_DIG_BITS)}};
			XA_VAL* r = xa_mul(v,&m);
			xa_free(v);
			v=r;
		}
        return *this;
    }

	XA_REF operator >> (int i) const
    {
		if (v && i>0)
		{
			// no field names cuz ms
			XA_VAL m = {/*sign:*/0,/*digs:*/1,/*quot:*/-1-i/(int)XA_DIG_BITS,/*refs:*/1,
			       /*data:*/{(XA_DIG)1<<(XA_DIG_BITS-(i%XA_DIG_BITS))}};
			XA_REF r;
			r.v = xa_mul(&m,v);
			return r;
		}
        return *this;
    }

	XA_REF& operator >>= (int i)
    {
		if (v && i>0)
		{
			// no field names cuz ms
			XA_VAL m = {/*sign:*/0,/*digs:*/1,/*quot:*/-1-i/(int)XA_DIG_BITS,/*refs:*/1,
				   /*data:*/{(XA_DIG)1<<(XA_DIG_BITS-(i%XA_DIG_BITS))}};
			XA_VAL* r = xa_mul(v,&m);
			xa_free(v);
			v=r;
		}
        return *this;
    }

	XA_REF& operator + ()
	{
		return *this;
	}

	XA_REF& operator - ()
	{
		if (v)
		{
			if (v->refs == 1)
				v->sign ^= 1;
			else
			{
				XA_VAL* n = xa_neg(v);
				xa_free(v);
				v = n;
			}
		}
		return *this;
	}

	XA_REF& operator -- ()
	{
		static const XA_VAL one = {/*sign:*/0,/*digs:*/1,/*quot:*/0,/*refs:*/1,/*data:*/{(XA_DIG)1}};
		XA_VAL* r = xa_sub(v,&one);
		xa_free(v);
		v = r;
		return *this;
	}

	XA_REF& operator ++ ()
	{
		static const XA_VAL one = {/*sign:*/0,/*digs:*/1,/*quot:*/0,/*refs:*/1,/*data:*/{(XA_DIG)1}};
		XA_VAL* r = xa_add(v,&one);
		xa_free(v);
		v = r;
		return *this;
	}

	XA_REF operator -- (int)
	{
		static const XA_VAL one{/*sign:*/0,/*digs:*/1,/*quot:*/0,/*refs:*/1,/*data:*/{(XA_DIG)1}};
		XA_REF r;
		r.v = v;
		v = xa_sub(v,&one);
		return r;
	}

	XA_REF operator ++ (int)
	{
		static const XA_VAL one{/*sign:*/0,/*digs:*/1,/*quot:*/0,/*refs:*/1,/*data:*/{(XA_DIG)1}};
		XA_REF r;
		r.v = v;
		v = xa_add(v,&one);
		return r;
	}

	XA_VAL* operator -> () const
	{
		static XA_VAL null{/*sign:*/0,/*digs:*/0,/*quot:*/0,/*refs:*/1,/*data:*/{(XA_DIG)0}};

		if (v)
			return v;
		else
		{
			null.quot = 0;
			return &null;
		}
	}

    // todo
    // friend std::istream& operator >> (std::istream& is, XA_REF& w);

	friend std::ostream& operator << (std::ostream& os, const XA_REF& w);

    private: XA_VAL* v;
};

/*
// todo
inline std::istream& operator >> (std::istream& is, W& w)
{
    std::istream::sentry s(is);
    if (s) 
	{
		while (is.good()) 
		{
			char c = is.get();
			if (std::isspace(c,is.getloc())) 
				break;
			if (std::isdigit(c,is.getloc())) 
				w.digits+=c;
		}
	}
    return is;
}
*/

inline std::ostream& operator << (std::ostream& os, const XA_REF& w)
{
	char* str;
	int len;

	if (os.flags() & os.dec)
		len = xa_extr_dec(w.v, &str);
	else
	if (os.flags() & os.hex)
		len = xa_extr_hex(w.v, &str);
	else
		return os; // unk format

	// prints garbage at the end! why??
	// os.write(str,len); 

	os << str; // works ok!

	free(str);
	return os;
}

#define IA_FAST // assumes rounding mode is set to DOWNWARD!
//#define XA_AUTO_TEST

//#include <math.h> // (nexttoward is incredibly slow)
#include <float.h>

struct IA_VAL
{
	IA_VAL() : lo(0), hi(0) {}
	IA_VAL(const IA_VAL& v) : lo(v.lo), hi(v.hi) {}

#ifndef IA_FAST
	IA_VAL(double d) : lo(d), hi(d) {}
#else
	IA_VAL(double d) : lo(d), hi(-d) {}
#endif

#ifndef IA_FAST
	static double inflate_neg_lo(double lo)
	{
		//return nexttoward(lo, -INFINITY);
		return lo * (1.0 + 1.0 * DBL_EPSILON) - DBL_MIN;
	}

	static double inflate_pos_lo(double lo)
	{
		//return nexttoward(lo, -INFINITY);
		return lo * (1.0 - 0.5 * DBL_EPSILON) - DBL_MIN;
	}

	static double inflate_neg_hi(double hi)
	{
		//return nexttoward(hi, +INFINITY);
		return hi * (1.0 - 0.5 * DBL_EPSILON) + DBL_MIN;
	}

	static double inflate_pos_hi(double hi)
	{
		//return nexttoward(hi, +INFINITY);
		return hi * (1.0 + 1.0 * DBL_EPSILON) + DBL_MIN;
	}

	inline void inflate()
	{
		#ifdef XA_AUTO_TEST
		double l = lo >= 0 ? inflate_pos_lo(lo) : inflate_neg_lo(lo);
		double h = hi >= 0 ? inflate_pos_hi(hi) : inflate_neg_hi(hi);
		assert(l<lo&& h>hi);
		lo = l;
		hi = h;
		#else

		if (lo >= 0)
		{
			lo = inflate_pos_lo(lo);
			hi = inflate_pos_hi(hi);
		}
		else
		if (hi < 0)
		{
			lo = inflate_neg_lo(lo);
			hi = inflate_neg_hi(hi);
		}
		else
		{
			lo = inflate_neg_lo(lo);
			hi = inflate_pos_hi(hi);
		}

		#endif
	}
	#endif // !IA_FAST

	operator double() const
	{
		#ifndef IA_FAST
		return (lo + hi) / 2;
		#else
		return (lo - hi) / 2;
		#endif
	}

	IA_VAL operator + (const IA_VAL& v) const
	{
		IA_VAL r;
		#ifndef IA_FAST
		r.lo = lo + v.lo;
		r.hi = hi + v.hi;
		r.inflate();
		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		#else
		r.lo = lo + v.lo;
		r.hi = hi + v.hi;
		#ifdef XA_AUTO_TEST
		assert(r.lo <= -r.hi);
		#endif
		#endif
		return r;
	}

	IA_VAL operator - (const IA_VAL& v) const
	{
		IA_VAL r;
		#ifndef IA_FAST
		r.lo = lo - v.hi;
		r.hi = hi - v.lo;
		r.inflate();
		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		#else
		r.lo = lo + v.hi;
		r.hi = hi + v.lo;
		#ifdef XA_AUTO_TEST
		assert(r.lo <= -r.hi);
		#endif
		#endif
		return r;
	}

	IA_VAL operator * (const IA_VAL& v) const
	{
		IA_VAL r;

		#ifndef IA_FAST
		if (lo >= 0)
		{
			if (v.lo >= 0)	// ++ * ++ = ++
			{
				r.lo = inflate_pos_lo(lo * v.lo);
				r.hi = inflate_pos_hi(hi * v.hi);
			}
			else
			if (v.hi < 0)	// ++ * -- = --
			{
				r.lo = inflate_neg_lo(hi * v.lo);
				r.hi = inflate_neg_hi(lo * v.hi);
			}
			else			// ++ * -+ = -+
			{
				r.lo = inflate_neg_lo(hi * v.lo);
				r.hi = inflate_pos_hi(hi * v.hi);
			}
		}
		else
		if (hi < 0)
		{
			if (v.lo >= 0)	// -- * ++ = --
			{
				r.lo = inflate_neg_lo(lo * v.hi);
				r.hi = inflate_neg_hi(hi * v.lo);
			}
			else
			if (v.hi < 0)	// -- * -- = ++
			{
				r.lo = inflate_pos_lo(hi * v.hi);
				r.hi = inflate_pos_hi(lo * v.lo);
			}
			else			// -- * -+ = -+
			{
				r.lo = inflate_neg_lo(lo * v.hi);
				r.hi = inflate_pos_hi(lo * v.lo);
			}
		}
		else
		{
			if (v.lo >= 0)	// -+ * ++ = -+
			{
				r.lo = inflate_neg_lo(lo * v.hi);
				r.hi = inflate_pos_hi(hi * v.hi);
			}
			else
			if (v.hi < 0)	// -+ * -- = -+
			{
				r.lo = inflate_neg_lo(hi * v.lo);
				r.hi = inflate_pos_hi(lo * v.lo);
			}
			else			// -+ * -+ = -+
			{
				r.lo = inflate_neg_lo(std::min(lo * v.hi, hi * v.lo));
				r.hi = inflate_pos_hi(std::max(lo * v.lo, hi * v.hi));
			}
		}

		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		
		#else

		// 3x -hi * v.hi
		// 3x -lo * v.lo
		// 3x -lo * v.hi + 1x lo * -v.hi
		// 3x -hi * v.lo + 1x hi * -v.lo

		// 1x lo * v.lo
		// 1x lo * v.hi
		// 1x hi * v.lo
		// 1x hi * v.hi

		// 2-4 conditional checks

		//int C = ((*(uint64_t*)&lo) >> 63) + ((*(uint64_t*)&hi) >> 63); // 0,1,2
		//C += 3 * ( ((*(uint64_t*)&v.lo) >> 63) + ((*(uint64_t*)&v.hi) >> 63) ); // + 3*(0,1,2)

		
		if (lo >= 0)
		{
			// negating hi

			if (v.lo >= 0)	// ++ * ++ = ++
			{
				//static int c = C; // 4
				//assert(c == C);
				r.lo = lo * v.lo;
				r.hi = -hi * v.hi;
			}
			else
			if (v.hi > 0)	// ++ * -- = --
			{
				//static int c = C; // 4
				//assert(c == C);
				r.lo = -hi * v.lo;
				r.hi = lo * v.hi;
			}
			else			// ++ * -+ = -+
			{
				//static int c = C;
				//assert(c == C);
				r.lo = -hi * v.lo;
				r.hi = -hi * v.hi;
			}
		}
		else
		if (hi >= 0)
		{
			// negating lo

			if (v.lo >= 0)	// -- * ++ = --
			{
				//static int c = C; // 4
				//assert(c == C);
				r.lo = -lo * v.hi;
				r.hi = hi * v.lo;
			}
			else
			if (v.hi >= 0)	// -- * -- = ++
			{
				//static int c = C; // 4
				//assert(c == C);
				r.lo = hi * v.hi;
				r.hi = -lo * v.lo;
			}
			else			// -- * -+ = -+
			{
				//static int c = C;
				//assert(c == C);
				r.lo = -lo * v.hi;
				r.hi = -lo * v.lo;
			}
		}
		else
		{
			// negating v.lo and v.hi

			if (v.lo >= 0)	// -+ * ++ = -+
			{
				//static int c = C;
				//assert(c == C);
				r.lo = lo * -v.hi;
				r.hi = hi * -v.hi;
			}
			else
			if (v.hi >= 0)	// -+ * -- = -+
			{
				//static int c = C;
				//assert(c == C);
				r.lo = hi * -v.lo;
				r.hi = lo * -v.lo;
			}
			else			// -+ * -+ = -+
			{
				//static int c = C;
				//assert(c == C);
				r.lo = std::min(lo * -v.hi, hi * -v.lo);
				//r.hi = -std::max(lo * v.lo, hi * v.hi);
				r.hi = std::min(lo * -v.lo, hi * -v.hi);
			}
		}

		#ifdef XA_AUTO_TEST
		assert(r.lo <= -r.hi);
		#endif

#endif

		return r;
	}

	bool operator < (int i) const
	{
		#ifndef IA_FAST
		return hi < i;
		#else
		//return -hi < i;
		return hi > -i;
		#endif
	}

	bool operator < (double d) const
	{
		#ifndef IA_FAST
		return hi < d;
		#else
		return -hi < d;
		#endif
	}

	bool operator < (const IA_VAL& v) const
	{
		#ifndef IA_FAST
		return hi < v.lo;
		#else
		return -hi < v.lo;
		#endif
	}

	bool operator > (int i) const
	{
		return lo > i;
	}

	bool operator > (double d) const
	{
		return lo > d;
	}

	bool operator > (const IA_VAL& v) const
	{
		#ifndef IA_FAST
		return lo > v.hi;
		#else
		return lo > -v.hi;
		#endif
	}

	double lo, hi;
};

#endif // c++
#endif // CRUDE_XA_H
