/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

//#define DELABELLA_LEGACY double
// define it to enable compatibility with older api
// you can use: float, double or long double
// can be defined just prior to including this header

template<typename T = double>
struct IDelaBella2
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

	static IDelaBella2<T>* Create();
	virtual ~IDelaBella2() = 0;

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

	// called right after Triangulate() and/or Constrain() it returns number of triangles,
	// but if called after Polygonize() it returns number of polygons
	virtual int GetNumPolygons() const = 0;

	virtual const Simplex* GetFirstDelaunaySimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Simplex* GetFirstHullSimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Vertex*  GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const Vertex*  GetFirstInternalVertex() const = 0;
	virtual const Vertex*  GetVertexByIndex(int i) const = 0;

	virtual int Constrain(int num, const int* pa, const int* pb, int advance_bytes) = 0;

	virtual int Polygonize(const Simplex* poly[/*GetNumOutputIndices()/3*/] = 0) = 0; // valid only if Triangulate() > 0

	// GenVoronoiDiagramVerts(), valid only if Triangulate() > 0
	// generates VD vertices (for use with VD edges or VD polys indices)
	// <x>,<y> and/or <indices> can be null if their values are not needed
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   V = P + N
	// <x>,<y> array must be (at least) V elements long
	// first P <x>,<y> elements will contain internal points (divisor W=1)
	// next N <x>,<y> elements will contain edge normals (divisor W=0)
	// function returns number vertices filled (V) on success otherwise -1
	virtual int GenVoronoiDiagramVerts(T* x, T* y, int advance_bytes = 0) const = 0;


	// GenVoronoiDiagramEdges(), valid only if Triangulate() > 0
	// generates unidirected VD edges (without ones in opposite direction)
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   I = 2 * (N + M + P - 1)
	// <indices> must be (at least) I elements long
	// every pair of consecutive values in <indices> represent VD edge
	// there is no guaranteed correspondence between edges order and other data
	// function returns number of indices filled (I) on success otherwise -1
	virtual int GenVoronoiDiagramEdges(int* indices, int advance_bytes = 0) const = 0;

	// GenVoronoiDiagramPolys() valid only if Triangulate() > 0
	// generates VD polygons
	// assuming:
	//   M = GetNumBoundaryVerts()
	//   N = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   I = 3 * (N + M) + 2 * (P - 1) + M
	// <indices> must be (at least) I elements long
	// if <poly_ending> value != 0 it will be inserted into <indices> after every poly
	//    otherwise every poly will be prefixed by number of its indices
	// first N polys in <indices> represent closed VD cells, thay are in order
	// and corresponding to vertices from GetFirstInternalVertex() -> Vertex::next list
	// next M polys in <indices> represent open VD cells, thay are in order
	// and corresponding to vertices from GetFirstBoundaryVertex() -> Vertex::next list
	// function returns number of indices filled (I) on success otherwise -1
	virtual int GenVoronoiDiagramPolys(int* indices, int advance_bytes=0, int poly_ending=~0) const = 0;

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
};

template <typename T>
inline const typename IDelaBella2<T>::Simplex* IDelaBella2<T>::Simplex::StartIterator(IDelaBella2<T>::Iterator* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

template <typename T>
inline const typename IDelaBella2<T>::Simplex* IDelaBella2<T>::Vertex::StartIterator(IDelaBella2<T>::Iterator* it/*not_null*/) const
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
typedef IDelaBella2<DELABELLA_LEGACY>  IDelaBella;
typedef IDelaBella::Simplex  DelaBella_Triangle;
typedef IDelaBella::Vertex   DelaBella_Vertex;
typedef IDelaBella::Iterator DelaBella_Iterator;
#endif

#endif
