/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

#include <intrin.h>
#include <xmmintrin.h>

//#define DELABELLA_LEGACY double
// define it to enable compatibility with older api
// you can use: float, double or long double
// can be defined just prior to including this header

// choose I wisely!
// int8_t  -> max points = 19
// int16_t -> max points = 4682
// int32_t -> max points = 306783379
// int64_t -> max points > 306783379 

template<typename T = double, typename I = int>
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
		I i; // index of original point

		inline const Simplex* StartIterator(Iterator* it/*not_null*/) const;
	};

	struct Simplex
	{
		Vertex* v[3];  // 3 vertices spanning this triangle
		Simplex* f[3]; // 3 adjacent faces, f[i] is at the edge opposite to vertex v[i]
		Simplex* next; // next triangle (of delaunay set or hull set)
		I index;       // list index
		
		uint8_t flags; 
		// 0x01 : edge shared with f[0] is fixed at least once
		// 0x02 : edge shared with f[1] is fixed at least once
		// 0x04 : edge shared with f[2] is fixed at least once
		// 0x08 : edge shared with f[0] is fixed odd numof times
		// 0x10 : edge shared with f[1] is fixed odd numof times
		// 0x20 : edge shared with f[2] is fixed odd numof times
		// 0x40 : face belongs to interior (valid after classify)
		// 0x80 : face belongs to hull (all other bits are meaningless)

		inline const Simplex* StartIterator(Iterator* it/*not_null*/, int around/*0,1,2*/) const;
	};

	static IDelaBella2<T,I>* Create();
	virtual ~IDelaBella2() = 0;

	virtual void Destroy() = 0;

	virtual void SetErrLog(int(*proc)(void* stream, const char* fmt, ...), void* stream) = 0;

	// return 0: no output 
	// negative: all points are colinear, output hull vertices form colinear segment list, no triangles on output
	// positive: output hull vertices form counter-clockwise ordered segment contour, delaunay and hull triangles are available
	// if 'y' pointer is null, y coords are treated to be located immediately after every x
	// if advance_bytes is less than 2*sizeof coordinate type, it is treated as 2*sizeof coordinate type  
	virtual I Triangulate(I points, const T* x, const T* y = 0, size_t advance_bytes = 0) = 0;

	// num of points passed to last call to Triangulate()
	virtual I GetNumInputPoints() const = 0;

	// num of indices returned from last call to Triangulate()
	virtual I GetNumOutputIndices() const = 0;

	// num of hull faces (non delaunay triangles)
	virtual I GetNumOutputHullFaces() const = 0;

	// num of boundary vertices
	virtual I GetNumBoundaryVerts() const = 0;

	// num of internal vertices
	virtual I GetNumInternalVerts() const = 0;

	// called right after Triangulate() and/or Constrain() it returns number of triangles,
	// but if called after Polygonize() it returns number of polygons
	virtual I GetNumPolygons() const = 0;

	virtual const Simplex* GetFirstDelaunaySimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Simplex* GetFirstHullSimplex() const = 0; // valid only if Triangulate() > 0
	virtual const Vertex*  GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const Vertex*  GetFirstInternalVertex() const = 0;
	virtual const Vertex*  GetVertexByIndex(I i) const = 0;

	// if classify=true return number of interior faces 
	// and links these faces in front of delaunay simplex list,
	// if classify=false all faces are classified as interior
	virtual I ConstrainEdges(I edges, const I* pa, const I* pb, size_t advance_bytes) = 0;

	virtual I Polygonize(const Simplex* poly[/*GetNumOutputIndices()/3*/] = 0) = 0; // valid only if Triangulate() > 0

	// GenVoronoiDiagramVerts(), valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates VD vertices (for use with VD edges or VD polys indices)
	// <x>,<y> can be null if only number of vertices is needed
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   V = P + N
	// <x>,<y> array must be (at least) V elements long
	// first P <x>,<y> elements will contain internal points (divisor W=1)
	// next N <x>,<y> elements will contain edge normals (divisor W=0)
	// function returns number vertices filled (V) on success, otherwise 0
	virtual I GenVoronoiDiagramVerts(T* x, T* y, size_t advance_bytes = 0) const = 0;


	// GenVoronoiDiagramEdges(), valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates unidirected VD edges (without ones in opposite direction)
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   I = 2 * (N + M + P - 1)
	// <indices> must be (at least) I elements long or must be null
	// every pair of consecutive values in <indices> represent VD edge
	// there is no guaranteed correspondence between edges order and other data
	// function returns number of indices filled (I) on success, otherwise 0
	virtual I GenVoronoiDiagramEdges(I* indices, size_t advance_bytes = 0) const = 0;

	// GenVoronoiDiagramPolys() valid only if Triangulate() > 0
	// it makes sense to call it prior to constraining only
	// generates VD polygons
	// assuming:
	//   N = GetNumBoundaryVerts()
	//   M = GetNumInternalVerts()
	//   P = GetNumPolygons()
	//   I = 3 * (N + M) + 2 * (P - 1) + N
	// <indices> must be (at least) I elements long or must be null
	// first M polys in <indices> represent closed VD cells, thay are in order
	// and corresponding to vertices from GetFirstInternalVertex() -> Vertex::next list
	// next N polys in <indices> represent open VD cells, thay are in order
	// and corresponding to vertices from GetFirstBoundaryVertex() -> Vertex::next list
	// every poly written to <indices> is terminated with ~0 value
	// if both <indices> and <closed_indices> are not null, 
	// number of closed VD cells indices is written to <closed_indices>
	// function returns number of indices filled (I) on success, otherwise 0
	virtual I GenVoronoiDiagramPolys(I* indices, size_t advance_bytes=0, I* closed_indices=0) const = 0;

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

template <typename T, typename I>
inline const typename IDelaBella2<T,I>::Simplex* IDelaBella2<T,I>::Simplex::StartIterator(IDelaBella2<T,I>::Iterator* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

template <typename T, typename I>
inline const typename IDelaBella2<T,I>::Simplex* IDelaBella2<T,I>::Vertex::StartIterator(IDelaBella2<T,I>::Iterator* it/*not_null*/) const
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
typedef IDelaBella2<DELABELLA_LEGACY,int>  IDelaBella;
typedef IDelaBella::Simplex  DelaBella_Triangle;
typedef IDelaBella::Vertex   DelaBella_Vertex;
typedef IDelaBella::Iterator DelaBella_Iterator;
#endif

#endif
