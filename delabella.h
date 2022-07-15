/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

//#define DELABELLA_AUTOTEST
// in case of troubles, allows to see if any assert pops up.
// define it globally (like with -DDELABELLA_AUTOTEST)

//#define DELABELLA_LEGACY double
// define it to enable compatibility with older api
// you can use: float, double or long double
// can be defined just prior to including this header

template<typename T = double>
struct IDelaBella3
{
	struct Vertex;
	struct Simplex;
	struct Iterator;

	struct Vertex
	{
		Vertex* next; // next in internal / boundary set of vertices
		Simplex* sew; // one of triangles sharing this vertex
		T x, y; // coordinates (input copy)
		int i; // index of original point

		inline const Simplex* StartIterator(Iterator* it/*not_null*/) const;
	};

	struct Simplex
	{
		Vertex* v[3]; // 3 vertices spanning this triangle
		Simplex* f[3]; // 3 adjacent faces, f[i] is at the edge opposite to vertex v[i]
		Simplex* next; // next triangle (of delaunay set or hull set)

		int index; // list index, if negative it is ~index'th in hull set 

		inline const Simplex* StartIterator(Iterator* it/*not_null*/, int around/*0,1,2*/) const;
	};

	struct Iterator
	{
		const Simplex* current;
		int around;

		const Simplex* Next()
		{
			int pivot = around + 1;
			if (pivot == 3)
				pivot = 0;

			Simplex* next = current->f[pivot];
			Vertex* v = current->v[around];

			if (next->v[0] == v)
				around = 0;
			else
			if (next->v[1] == v)
				around = 1;
			else
				around = 2;

			current = next;
			return current;
		}

		const Simplex* Prev()
		{
			int pivot = around - 1;
			if (pivot == -1)
				pivot = 2;

			Simplex* prev = current->f[pivot];
			Vertex* v = current->v[around];

			if (prev->v[0] == v)
				around = 0;
			else
			if (prev->v[1] == v)
				around = 1;
			else
				around = 2;

			current = prev;
			return current;
		}
	};

	static IDelaBella3<T>* Create();
	virtual ~IDelaBella3() = 0;

	virtual void Destroy() = 0;

	virtual void SetErrLog(int(*proc)(void* stream, const char* fmt, ...), void* stream) = 0;

	// return 0: no output 
	// negative: all points are colinear, output hull vertices form colinear segment list, no triangles on output
	// positive: output hull vertices form counter-clockwise ordered segment contour, delaunay and hull triangles are available
	// if 'y' pointer is null, y coords are treated to be located immediately after every x
	// if advance_bytes is less than 2*sizeof coordinate type, it is treated as 2*sizeof coordinate type  
	virtual int Triangulate(int points, const T* x, const T* y = 0, int advance_bytes = 0) = 0;

	// num of points passed to last call to Triangulate()
	virtual int GetNumInputPoints() const = 0;

	// num of indices returned from last call to Triangulate()
	virtual int GetNumOutputIndices() const = 0;

	// num of hull faces (non delaunay triangles)
	virtual int GetNumOutputHullFaces() const = 0;

	// num of boundary vertices
	virtual int GetNumBoundaryVerts() const = 0;

	// num of internal vertices
	virtual int GetNumInternalVerts() const = 0;

	virtual const Simplex* GetFirstDelaunaySimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Simplex* GetFirstHullSimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Vertex*  GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const Vertex*  GetFirstInternalVertex() const = 0;
	virtual const Vertex*  GetVertexByIndex(int i) const = 0;

	virtual int Constrain(int num, const int* pa, const int* pb, int advance_bytes) = 0;

	virtual int Polygonize(const Simplex* poly[/*GetNumOutputIndices()/3*/] = 0) = 0; // valid only if Triangulate() > 0

	#ifdef DELABELLA_LEGACY
	inline const Simplex* GetFirstDelaunayTriangle() const
	{
		return GetFirstDelaunaySimplex();
	}
	inline const Simplex* GetFirstHullTriangle() const
	{
		return GetFirstHullSimplex();
	}
	#endif
};

template <typename T>
inline const typename IDelaBella3<T>::Simplex* IDelaBella3<T>::Simplex::StartIterator(IDelaBella3<T>::Iterator* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

template <typename T>
inline const typename IDelaBella3<T>::Simplex* IDelaBella3<T>::Vertex::StartIterator(IDelaBella3<T>::Iterator* it/*not_null*/) const
{
	it->current = sew;
	if (sew->v[0] == this)
		it->around = 0;
	else
	if (sew->v[1] == this)
		it->around = 1;
	else
		it->around = 2;
	return sew;
}

#ifdef DELABELLA_LEGACY
typedef IDelaBella3<DELABELLA_LEGACY>  IDelaBella;
typedef IDelaBella::Simplex  DelaBella_Triangle;
typedef IDelaBella::Vertex   DelaBella_Vertex;
typedef IDelaBella::Iterator DelaBella_Iterator;
#endif

#endif
