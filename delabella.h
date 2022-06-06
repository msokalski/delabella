/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

#ifdef CRUDE_XA
// exact arithmetic floating point configuration
// its about 50x slower than regular floating point setup
// but if you need it it's priceless!
#include "crude-xa/src/crude-xa.h"
typedef XA_REF Signed14;   // BITS			xy coords
typedef XA_REF Signed15;   // BITS + 1		vect::xy
typedef XA_REF Unsigned28; // 2xBITS		z coord
typedef XA_REF Signed29;   // 2xBITS + 1	vect::z
typedef XA_REF Signed31;   // 2xBITS + 3	norm::z
typedef XA_REF Signed45;   // 3xBITS + 3	norm::xy
typedef XA_REF Signed62;   // 4xBITS + 6	dot(vect,norm)
#else
// regular floating point setup
/*
typedef double Signed14;		// BITS			xy coords
typedef double Signed15;		// BITS + 1		vect::xy
typedef long double Unsigned28; // 2xBITS		z coord
typedef long double Signed29;   // 2xBITS + 1	vect::z
typedef long double Signed31;   // 2xBITS + 3	norm::z
typedef long double Signed45;   // 3xBITS + 3	norm::xy
typedef long double Signed62;   // 4xBITS + 6	dot(vect,norm)
*/

typedef double Signed14;		// BITS			xy coords
typedef double Signed15;		// BITS + 1		vect::xy
typedef double Unsigned28; // 2xBITS		z coord
typedef double Signed29;   // 2xBITS + 1	vect::z
typedef double Signed31;   // 2xBITS + 3	norm::z
typedef double Signed45;   // 3xBITS + 3	norm::xy
typedef double Signed62;   // 4xBITS + 6	dot(vect,norm)
#endif

/*
// exact arithmetic using integers 
// can handle only a tiny range of input coords +/-8191
typedef int16_t  Signed14;		// BITS			xy coords
typedef int16_t  Signed15;		// BITS + 1		vect::xy
typedef uint32_t Unsigned28;	// 2xBITS		z coord
typedef int32_t  Signed29;		// 2xBITS + 1	vect::z
typedef int32_t  Signed31;		// 2xBITS + 3	norm::z
typedef int64_t  Signed45;		// 3xBITS + 3	norm::xy
typedef int64_t  Signed62;		// 4xBITS + 6	dot(vect,norm)
*/

// returns: positive value: number of triangle indices, negative: number of line segment indices (degenerated input)
//          triangle indices in abc array are always returned in clockwise order
// DEPRECIATED. move to new API either extern "C" or IDelaBella (C++)
int DelaBella(int points, const double* xy/*[points][2]*/, int* abc/*[2*points-5][3]*/, int (*errlog)(const char* fmt,...) = printf);

struct DelaBella_Triangle;
struct DelaBella_Iterator;

struct DelaBella_Vertex
{
	int i; // index of original point
	Signed14 x,y; // coordinates (input copy)
	DelaBella_Vertex* next; // next in internal / boundary set of vertices
	DelaBella_Triangle* sew; // one of triangles sharing this vertex

	inline const DelaBella_Triangle* StartIterator(DelaBella_Iterator* it/*not_null*/) const;
};

struct DelaBella_Triangle
{
	DelaBella_Vertex* v[3]; // 3 vertices spanning this triangle
	DelaBella_Triangle* f[3]; // 3 adjacent faces, f[i] is at the edge opposite to vertex v[i]
	DelaBella_Triangle* next; // next triangle (of delaunay set or hull set)

	struct Circumcenter
	{
		~Circumcenter() {} // prevents leaks from crude-xa
		Signed45 x; // voronoi vertex = {x/(-2*z),y/(-2*z)}
		Signed45 y; // keeping it as a ratio avoids any rounding
		Signed31 z; // additionally, if z<0 this is delaunay triangle
	} n; 

	int index; // list index

/*
               v[1]
    +-----------+-----------+
	 \         /.\         /
      \ f[0]  /. .\  f[2] /
 	   \     /. . .\  _--/--> iterator.around(0)->prev
        \   /. this.\/  /
	     \ /. . . . /\ / 
	 v[2] +--------\--+ v[0]
	       \        \/ 
            \  f[1] /'--> iterator.around(0)->next
             \     /
              \   /
			   \ /
                +
*/


	inline const DelaBella_Triangle* StartIterator(DelaBella_Iterator* it/*not_null*/, int around/*0,1,2*/) const;
};

struct DelaBella_Iterator
{
	const DelaBella_Triangle* current;
	int around;

	const DelaBella_Triangle* Next()
	{
		int pivot = around+1;
		if (pivot == 3)
			pivot = 0;

		DelaBella_Triangle* next = current->f[pivot];

		// we should not pass silhouette!
		// ...

		DelaBella_Vertex* v = current->v[around];

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

	const DelaBella_Triangle* Prev()
	{
		int pivot = around-1;
		if (pivot == -1)
			pivot = 2;

		DelaBella_Triangle* prev = current->f[pivot];

		// we should not pass silhouette!
		// ...

		DelaBella_Vertex* v = current->v[around];

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

inline const DelaBella_Triangle* DelaBella_Triangle::StartIterator(DelaBella_Iterator* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

inline const DelaBella_Triangle* DelaBella_Vertex::StartIterator(DelaBella_Iterator* it/*not_null*/) const
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

#ifdef __cplusplus
struct IDelaBella
{
	static IDelaBella* Create();
	virtual ~IDelaBella();

	virtual void Destroy() = 0;

	virtual void SetErrLog(int(*proc)(void* stream, const char* fmt, ...), void* stream) = 0;

	// return 0: no output 
	// negative: all points are colinear, output hull vertices form colinear segment list, no triangles on output
	// positive: output hull vertices form counter-clockwise ordered segment contour, delaunay and hull triangles are available
	// if 'y' pointer is null, y coords are treated to be located immediately after every x
	// if advance_bytes is less than 2*sizeof coordinate type, it is treated as 2*sizeof coordinate type  
	virtual int Triangulate(int points, const float* x, const float* y = 0, int advance_bytes = 0) = 0;
	virtual int Triangulate(int points, const double* x, const double* y = 0, int advance_bytes = 0) = 0;
	virtual int Triangulate(int points, const long double* x, const long double* y = 0, int advance_bytes = 0) = 0;

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

	virtual const DelaBella_Triangle* GetFirstDelaunayTriangle() const = 0; // valid only if Triangulate() > 0
	virtual const DelaBella_Triangle* GetFirstHullTriangle() const = 0; // valid only if Triangulate() > 0
	virtual const DelaBella_Vertex*   GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const DelaBella_Vertex*   GetFirstInternalVertex() const = 0;

	// search api
	// virtual const DelaBella_Query* CreateHashedVoronoiQuadTree();

	// Voronoi sub-api

	// legend:
	// dela_verts: number of input points, optionally excluding duplicated ones (num_verts from GetVertexArray())
	// dela_tris: number of returned triangles (verts/3) from Triangulate() or GetNumOutputIndices()
	// boundary_verts: GetNumBoundaryVerts()
	
	virtual int GetNumVoronoiVertices(bool center) const { return 0; }
	// = dela_tris + boundary_verts + center * dela_verts
	
	virtual int GetNumVoronoiIndices(bool center, bool loop_open, bool loop_closed, bool suffix) const { return 0; } 
	// = 2 * (dela_verts + dela_tris - 1) + boundary_verts 
	//   + (center + suffix) * dela_verts 
	//   + loop_closed * (dela_verts - boundary_verts) 
	//   + loop_open * boundary_verts

	virtual void BuildVoronoiDiagram(
		float* x, float* y,    // [GetNumVoronoiVertices()]
		int advance_bytes,     // x,y array advance bytes
		int* index_buffer,     // [GetNumVoronoiIndices()]
		int* num_site_indices, // [dela_verts], fills number of indices for every cell, (can be null)
		bool center,           // add dela vertex (voro site center), into xy and index_buffer as each site's first indice 
		bool loop_open,        // repeat looping indice after every open (boundary)site
		bool loop_closed,      // repeat looping indice after every closed (internal) site
		bool suffix            // add '-1' suffix after every site
	) const {}			
};
#else
void* DelaBella_Create();
void  DelaBella_Destroy(void* db);
void  DelaBella_SetErrLog(void* db, int(*proc)(void* stream, const char* fmt, ...), void* stream);
int   DelaBella_TriangulateFloat(void* db, int points, float* x, float* y = 0, int advance_bytes = 0);
int   DelaBella_TriangulateDouble(void* db, int points, double* x, double* y = 0, int advance_bytes = 0);
int   DelaBella_TriangulateLongDouble(void* db, int points, long double* x, long double* y = 0, int advance_bytes = 0);
int   DelaBella_GetNumInputPoints(void* db);
int   DelaBella_GetNumOutputIndices(void* db);
int   Delabella_GetNumOutputHullFaces(void* db);
int   Delabella_GetNumOutputHullVerts(void* db);
const DelaBella_Triangle* GetFirstDelaunayTriangle(void* db);
const DelaBella_Triangle* GetFirstHullTriangle(void* db);
const DelaBella_Vertex*   GetFirstHullVertex(void* db);
#endif

#endif
