/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#ifndef DELABELLA_H
#define DELABELLA_H

//#define DB_AUTO_TEST

template <typename> struct DelaBella_Triangle;
template <typename> struct DelaBella_Iterator;

template <typename T = double>
struct DelaBella_Vertex
{
	DelaBella_Vertex* next; // next in internal / boundary set of vertices
	DelaBella_Triangle<T>* sew; // one of triangles sharing this vertex
	T x, y; // coordinates (input copy)
	int i; // index of original point

	inline const DelaBella_Triangle<T>* StartIterator(DelaBella_Iterator<T>* it/*not_null*/) const;
};

template <typename T = double>
struct DelaBella_Triangle
{
	DelaBella_Vertex<T>* v[3]; // 3 vertices spanning this triangle
	DelaBella_Triangle* f[3]; // 3 adjacent faces, f[i] is at the edge opposite to vertex v[i]
	DelaBella_Triangle* next; // next triangle (of delaunay set or hull set)

	int index; // list index, if negative it is ~index'th in hull set 

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

	inline const DelaBella_Triangle* StartIterator(DelaBella_Iterator<T>* it/*not_null*/, int around/*0,1,2*/) const;
};

template <typename T = double>
struct DelaBella_Iterator
{
	const DelaBella_Triangle<T>* current;
	int around;

	const DelaBella_Triangle<T>* Next()
	{
		int pivot = around+1;
		if (pivot == 3)
			pivot = 0;

		DelaBella_Triangle<T>* next = current->f[pivot];
		DelaBella_Vertex<T>* v = current->v[around];

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

	const DelaBella_Triangle<T>* Prev()
	{
		int pivot = around-1;
		if (pivot == -1)
			pivot = 2;

		DelaBella_Triangle<T>* prev = current->f[pivot];
		DelaBella_Vertex<T>* v = current->v[around];

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

template <typename T = double>
inline const DelaBella_Triangle<T>* DelaBella_Triangle<T>::StartIterator(DelaBella_Iterator<T>* it/*not_null*/, int around/*0,1,2*/) const
{
	it->current = this;
	it->around = around;
	return this;
}

template <typename T = double>
inline const DelaBella_Triangle<T>* DelaBella_Vertex<T>::StartIterator(DelaBella_Iterator<T>* it/*not_null*/) const
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

template<typename T = double>
struct IDelaBella
{
	static IDelaBella<T>* Create();
	virtual ~IDelaBella() = 0;

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

	virtual const DelaBella_Triangle<T>* GetFirstDelaunayTriangle() const = 0; // valid only if Triangulate() > 0
	virtual const DelaBella_Triangle<T>* GetFirstHullTriangle() const = 0; // valid only if Triangulate() > 0
	virtual const DelaBella_Vertex<T>*   GetFirstBoundaryVertex() const = 0; // if Triangulate() < 0 it is list, otherwise closed contour! 
	virtual const DelaBella_Vertex<T>*   GetFirstInternalVertex() const = 0;
	virtual const DelaBella_Vertex<T>*   GetVertexByIndex(int i) const = 0;

	virtual int Constrain(int num, const int* pa, const int* pb, int advance_bytes) = 0;

	virtual int Polygonize(const DelaBella_Triangle<T>* poly[/*GetNumOutputIndices()/3*/] = 0) = 0; // valid only if Triangulate() > 0
};

#endif