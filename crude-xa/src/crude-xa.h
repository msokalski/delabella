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

#include <math.h>
struct IA_VAL
{
	IA_VAL() : lo(0), hi(0) {}
	IA_VAL(double d) : lo(d), hi(d) {}
	IA_VAL(const IA_VAL& v) : lo(v.lo), hi(v.hi) {}

	#define lower (1.0 - 0.5 * DBL_EPSILON)
	#define upper (1.0 + 1.0 * DBL_EPSILON)

	inline void explode()
	{
		/*
		lo = lo >= 0 ? lo * lower - DBL_MIN : lo * upper - DBL_MIN;
		hi = hi >= 0 ? hi * upper + DBL_MIN : hi * lower + DBL_MIN;
		*/
		// test
		double l = lo >= 0 ? lo * lower - DBL_MIN : lo * upper - DBL_MIN;
		double h = hi >= 0 ? hi * upper + DBL_MIN : hi * lower + DBL_MIN;
		assert(l<lo && h>hi);
		lo = l;
		hi = h;
	}

	operator double() const
	{
		return (lo + hi) / 2;
	}

	IA_VAL operator + (const IA_VAL& v) const
	{
		IA_VAL r;
		r.lo = lo + v.lo;
		r.hi = hi + v.hi;
		r.explode();
		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		return r;
	}

	IA_VAL operator - (const IA_VAL& v) const
	{
		IA_VAL r;
		r.lo = lo - v.hi;
		r.hi = hi - v.lo;
		r.explode();
		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		return r;
	}

	IA_VAL operator * (const IA_VAL& v) const
	{
		IA_VAL r;

		if (lo >= 0)
		{
			if (v.lo >= 0)	// +    +    +    +  ->  +    +
			{
				r.lo = lo * v.lo * lower - DBL_MIN;
				r.hi = hi * v.hi * upper + DBL_MIN;
			}
			else
			if (v.hi < 0)	// +    +    -    -  ->  -    -
			{
				//r.lo = hi * v.hi * upper - DBL_MIN;
				r.lo = hi * v.lo * upper - DBL_MIN;
				r.hi = lo * v.hi * lower + DBL_MIN;
			}
			else			// +    +    -    +  ->  -    +
			{
				r.lo = hi * v.lo * upper - DBL_MIN;
				r.hi = hi * v.hi * upper + DBL_MIN;
			}
		}
		else
		if (hi < 0)
		{
			if (v.lo >= 0)	// -    -    +    +  ->  -    -
			{
				r.lo = lo * v.hi * upper - DBL_MIN;
				r.hi = hi * v.lo * lower + DBL_MIN;
			}
			else
			if (v.hi < 0)	// -    -    -    -  ->  +    +
			{
				r.lo = hi * v.hi * lower - DBL_MIN;
				r.hi = lo * v.lo * upper + DBL_MIN;
			}
			else			// -    -    -    +  ->  -    +
			{
				r.lo = lo * v.hi * upper - DBL_MIN;
				//r.hi = hi * v.lo * upper + DBL_MIN;
				r.hi = lo * v.lo * upper + DBL_MIN;
			}
		}
		else
		{
			if (v.lo >= 0)	// -    +    +    +  ->  -    +
			{
				r.lo = lo * v.hi * upper - DBL_MIN;
				r.hi = hi * v.hi * upper + DBL_MIN;
			}
			else
			if (v.hi < 0)	// -    +    -    -  ->  -    +
			{
				r.lo = hi * v.lo * upper - DBL_MIN;
				r.hi = lo * v.lo * upper + DBL_MIN;
			}
			else			// -    +    -    +  ->  -    +
			{
				r.lo = std::min(lo * v.hi, hi * v.lo) * upper - DBL_MIN;
				r.hi = std::max(lo * v.lo, hi * v.hi) * upper + DBL_MIN;
			}
		}

		#ifdef XA_AUTO_TEST
		assert(r.lo < r.hi);
		#endif
		return r;
	}

	bool operator < (int i) const
	{
		return hi < i;
	}

	bool operator < (double d) const
	{
		return hi < d;
	}

	bool operator < (const IA_VAL& v) const
	{
		return hi < v.lo;
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
		return lo > v.hi;
	}

	double lo, hi;
};

#endif