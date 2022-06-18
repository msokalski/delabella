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

	const XA_REF eval() const
	{
		return *this;
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

template <typename A, typename B>
struct XA_ADD
{
	XA_ADD(const A& a, const B& b) : a(a), b(b) {}
	XA_REF eval() { return a.eval() + b.eval(); }
	A a;
	B b;
};

template <typename A, typename B>
struct XA_SUB
{
	XA_SUB(const A& a, const B& b) : a(a), b(b) {}
	XA_REF eval() { return a.eval() - b.eval(); }
	A a;
	B b;
};

template <typename A, typename B>
struct XA_MUL
{
	XA_MUL(const A& a, const B& b) : a(a), b(b) {}
	XA_REF eval() { return a.eval() * b.eval(); }
	A a;
	B b;
};

#endif