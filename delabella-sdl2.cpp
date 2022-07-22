/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#define _CRT_SECURE_NO_WARNINGS

#define DELABELLA_LEGACY double

#define VORONOI
//#define VORONOI_POLYS
// otherwise EDGES

// override build define
//#undef WITH_DELAUNATOR 
//#define WITH_DELAUNATOR

struct int2048
{
    /*
    int2048()
    {
        memset(bits, 0, sizeof(bits));
    }

    int2048& operator += (double d)
    {
        // extract sign, exp, mantissa

    }

    void extract(uint8_t doubles, double* data)
    {
        // bitscan to locate exponent
    }

    uint32_t bits[64];
    */
};

// override build define
#undef WITH_CDT
#define WITH_CDT

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include <vector>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <random>

#include "delabella.h"

#undef GL_GLEXT_VERSION // silence MAC vs SDL2 fight
#include <GL/gl.h>

#ifndef _WIN32
#include <GL/glext.h>
#else
#pragma comment(lib,"OpenGL32.lib")
#pragma comment(lib,"SDL2.lib")
#undef main // on windows SDL does this weird thing
#endif

// competitors
#ifdef WITH_DELAUNATOR
#include "delaunator/delaunator-header-only.hpp"
#endif

#ifdef WITH_CDT
#include "CDT/include/CDT.h"

namespace CDT
{
typedef std::vector<TriInd> Poly;
template <typename T, typename L = LocatorKDTree<T> >
std::vector<Poly> Polygonize(Triangulation<T, L>& cdt)
{
    typedef TriInd PolyInd;
    // for every tri, here we store index of the poly the tri belongs to
    auto map = std::vector<PolyInd>(cdt.triangles.size());
    // vector of Polys, each as a vector of tri indices belonging to the Poly
    auto polys = std::vector<Poly>();

    for (TriInd iT = 0; iT < (TriInd)cdt.triangles.size(); iT++)
    {
        const auto& tri = cdt.triangles[iT];
        auto merged = false;

        // compare i'th tri with its adjacent tris
        for (const auto adj : tri.neighbors)
        {
            // but only if there's adjacent tri and it is processed already
            if (adj == noNeighbor || adj > iT)
                continue;
            // locate reflex vert in adj triangle
            const auto& vr =
                cdt.vertices[opposedVertex(cdt.triangles[adj], iT)];
            const auto& v1 = cdt.vertices[tri.vertices[0]];
            const auto& v2 = cdt.vertices[tri.vertices[1]];
            const auto& v3 = cdt.vertices[tri.vertices[2]];
            using predicates::adaptive::incircle;
            if (!incircle(vr.x, vr.y, v1.x, v1.y, v2.x, v2.y, v3.x, v3.y))
            {
                if (!merged)
                {
                    // append tri to already existing poly
                    merged = true;
                    const auto append_to = map[adj];
                    map[iT] = append_to;
                    polys[append_to].push_back(iT);
                }
                else
                {
                    const auto merge_to = map[iT];
                    const auto merge_from = map[adj];

                    if (merge_to == merge_from)
                        continue;

                    // funny case, tri is a bridge between 2 polys merge'em all
                    // together
                    for (const auto i : polys[merge_from])
                    {
                        map[i] = merge_to;            // remap
                        polys[merge_to].push_back(i); // merge
                    }
                    if (merge_from != (PolyInd)polys.size() - 1)
                    {
                        // replace merge_from poly with last poly in polys
                        polys[merge_from] = polys.back();
                        for (const auto i : polys[merge_from])
                        {
                            map[i] = merge_from; // remap
                        }
                    }
                    polys.pop_back();
                }
            }
        }

        if (!merged)
        {
            // at the moment, just alone tri
            // make a new poly for it
            map[iT] = (PolyInd)polys.size();
            polys.push_back({ iT });
        }
    }

    // post proc

    struct Stepper
    {
        const Triangulation<T, L>& cdt;
        TriInd current;
        int around;

        Stepper(const Triangulation<T, L>& cdt, TriInd t, int a) : cdt(cdt), current(t), around(a) {}
            
        void StepOver(int a) 
        { 
            around = a; 
        }

        TriInd Clockwise()
        {
            const auto& prev = cdt.triangles[current];

            current = prev.neighbors[around];
            const auto& next = cdt.triangles[current];

            VertInd v = prev.vertices[around];
            if (next.vertices[0] == v)
                around = 0;
            else
            if (next.vertices[1] == v)
                around = 1;
            else
                around = 2;

            return current;
        }
    };

    const PolyInd mask = (~(PolyInd)0) >> 1;

    for (PolyInd p = 0; p < polys.size(); p++)
    {
        const PolyInd q = p | ~mask; // unmarked

        // single triangle polys are ok already
        if (polys[p].size() == 1)
            continue;

        // find good starting triangle,
        // one with exeactly 1 inner edge
        TriInd first = noNeighbor;
        for (const auto t : polys[p])
        {
            if (first == noNeighbor)
            {
                const auto& tri = cdt.triangles[t];
                int inner_edges =
                    (tri.neighbors[0] != noNeighbor && (map[tri.neighbors[0]] & mask) == p) +
                    (tri.neighbors[1] != noNeighbor && (map[tri.neighbors[1]] & mask) == p) +
                    (tri.neighbors[2] != noNeighbor && (map[tri.neighbors[2]] & mask) == p);

                if (inner_edges == 1)
                    first = t;
            }

            // mark all tris as not inserted
            map[t] = q;
        }

        // we can clear current poly now, 
        // as we depend only on map and adjacency
        polys[p].clear();

        TriInd f = first; // current face
        bool step_on = false; // is current vertex inserted
        int insert = 2; // first triangle should end with 2

        Stepper it(cdt, f, 
            cdt.triangles[f].neighbors[0] != noNeighbor && map[cdt.triangles[f].neighbors[0]] == q ? 0 :
            cdt.triangles[f].neighbors[1] != noNeighbor && map[cdt.triangles[f].neighbors[1]] == q ? 1 : 2);

        while (1)
        {
            if (!step_on && map[f] == q)
            {
                step_on = true;
                map[f] = p; // mark as inserted

                if (it.around != insert)
                {
                    auto& tri = cdt.triangles[f];
                    static const int rot[3][3] = { {0,1,2},{2,0,1},{1,2,0} };
                    const int r = rot[it.around][insert];
                    const auto v = tri.vertices[r];
                    const auto n = tri.neighbors[r];
                    switch (r)
                    {
                    case 1:
                        tri.vertices[1] = tri.vertices[0];
                        tri.vertices[0] = tri.vertices[2];
                        tri.vertices[2] = v;
                        tri.neighbors[1] = tri.neighbors[0];
                        tri.neighbors[0] = tri.neighbors[2];
                        tri.neighbors[2] = n;
                        break;
                    case 2:
                        tri.vertices[2] = tri.vertices[0];
                        tri.vertices[0] = tri.vertices[1];
                        tri.vertices[1] = v;
                        tri.neighbors[2] = tri.neighbors[0];
                        tri.neighbors[0] = tri.neighbors[1];
                        tri.neighbors[1] = n;
                        break;
                    default:
                        break;
                    }
                    it.StepOver(insert);
                }

                polys[p].push_back(f);
                insert = 0; // everything but first should use 0
            }

            TriInd probe = cdt.triangles[f].neighbors[it.around];
            if (probe == noNeighbor || (map[probe] & mask) != p)
            {
                // check if we've covered current vertex 
                // with some face before stepping over
                assert(step_on); 

                // we're on last tri inside poly (marked or unmarked)
                // step on other leg:
                static const int other_leg[3] = { 1,2,0 };
                it.StepOver(other_leg[it.around]);
                step_on = false;
                continue;
            }

            f = it.Clockwise();
            if (f == first)
                break;
        }
    }

    return polys;
}
}

#endif

PFNGLBINDBUFFERPROC glBindBuffer = 0;
PFNGLDELETEBUFFERSPROC glDeleteBuffers = 0;
PFNGLGENBUFFERSPROC glGenBuffers = 0;
PFNGLBUFFERDATAPROC glBufferData = 0;
PFNGLBUFFERSUBDATAPROC glBufferSubData = 0;
PFNGLMAPBUFFERPROC glMapBuffer = 0;
PFNGLUNMAPBUFFERPROC glUnmapBuffer = 0;
PFNGLPRIMITIVERESTARTINDEXPROC glPrimitiveRestartIndex = 0;

PFNGLCREATESHADERPROC  glCreateShader = 0;
PFNGLDELETESHADERPROC  glDeleteShader = 0;
PFNGLCREATEPROGRAMPROC glCreateProgram = 0;
PFNGLDELETEPROGRAMPROC glDeleteProgram = 0;
PFNGLSHADERSOURCEPROC  glShaderSource = 0;
PFNGLCOMPILESHADERPROC glCompileShader = 0;
PFNGLATTACHSHADERPROC  glAttachShader = 0;
PFNGLDETACHSHADERPROC  glDetachShader = 0;
PFNGLLINKPROGRAMPROC   glLinkProgram = 0;
PFNGLUSEPROGRAMPROC    glUseProgram = 0;

PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog = 0;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog = 0;

PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation = 0;
PFNGLUNIFORM4FPROC glUniform4f = 0;
//PFNGLUNIFORMMATRIX4DVPROC glUniformMatrix4dv = 0;
PFNGLUNIFORM4DVPROC glUniform4dv = 0;
PFNGLVERTEXATTRIBLPOINTERPROC glVertexAttribLPointer = 0;
PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray = 0;
PFNGLDISABLEVERTEXATTRIBARRAYPROC glDisableVertexAttribArray = 0;

bool BindGL(bool prim_restart, bool shaders)
{
	#define BINDGL(proc) if ((*(void**)&proc = SDL_GL_GetProcAddress(#proc)) == 0) return false;
    
	BINDGL(glBindBuffer);
	BINDGL(glDeleteBuffers);
	BINDGL(glGenBuffers);
	BINDGL(glBufferData);
	BINDGL(glBufferSubData);
	BINDGL(glMapBuffer);
	BINDGL(glUnmapBuffer);

    if (shaders)
    {
        BINDGL(glCreateShader);
        BINDGL(glDeleteShader);
        BINDGL(glCreateProgram);
        BINDGL(glDeleteProgram);
        BINDGL(glShaderSource);
        BINDGL(glCompileShader);
        BINDGL(glAttachShader);
        BINDGL(glDetachShader);
        BINDGL(glLinkProgram);
        BINDGL(glUseProgram);
        BINDGL(glGetShaderInfoLog);
        BINDGL(glGetProgramInfoLog);
        //BINDGL(glUniformMatrix4dv);
        BINDGL(glUniform4dv);
        BINDGL(glUniform4f);
        BINDGL(glVertexAttribLPointer);
        BINDGL(glEnableVertexAttribArray);
        BINDGL(glDisableVertexAttribArray);
        BINDGL(glGetUniformLocation);
    }

    if (prim_restart)
	    BINDGL(glPrimitiveRestartIndex);

	#undef BINDGL
	return true;
}


int errlog(void* stream, const char* fmt, ...)
{
	va_list arg;
	va_start(arg,fmt);
	int ret = vfprintf((FILE*)stream, fmt, arg);
	va_end(arg);
    fflush((FILE*)stream);
	return ret;
}

static uint64_t uSec()
{
	#ifdef _WIN32
	LARGE_INTEGER c;
	static LARGE_INTEGER f;
	static BOOL bf = QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&c);
	uint64_t n = c.QuadPart;
	uint64_t d = f.QuadPart;
	uint64_t m = 1000000;
	// calc microseconds = n*m/d carefully!
	// naive mul/div would work only for upto 5h on 1GHz freq
	// we exploit fact that m*d fits in uint64 (upto 18THz freq)
	// so n%d*m fits as well,
	return n / d * m + n % d * m / d;
	#else
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (uint64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
	#endif
}

// wrap GL buffer, 
// so mapping works even if it doesnt
struct Buf
{
    GLuint buf;
    GLenum target;
    GLsizei size;
    void* map;
    bool mapped;

    Buf() : buf(0),target(0),size(0),map(0),mapped(false) {}

    GLuint Gen(GLenum t, GLsizei s)
    {
        target = t;
        size = s;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glGenBuffers(1, &buf);
        glBindBuffer(target, buf);
        glBufferData(target, size, 0, GL_STATIC_DRAW);

        glBindBuffer(target, push);
        return buf;
    }

    void Del()
    {
        if (buf)
        {
            Unmap();
            glDeleteBuffers(1,&buf);
            buf = 0;
        }
    }

    void* Map()
    {
        if (map || !buf)
            return 0;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        map = glMapBuffer(target, GL_WRITE_ONLY);
        mapped = true;
        if (!map)
        {
            map = malloc(size);
            mapped = false;
        }

        glBindBuffer(target, push);
        return map;
    }

    void Unmap()
    {
        if (!map || !buf)
            return;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        if (mapped)
            glUnmapBuffer(target);
        else
        {
            glBufferSubData(target, 0, size, map);
            free(map);
        }

        map = 0;
        mapped = false;

        glBindBuffer(target, push);
    }

    void Bind()
    {
        glBindBuffer(target,buf);
    }
};

typedef DELABELLA_LEGACY MyCoord;

struct MyPoint
{
    MyCoord x;
    MyCoord y;

    // operators used by result comparator

    bool operator == (const MyPoint& p) const
    {
        return x == p.x && y == p.y;
    }

    bool operator < (const MyPoint& p) const
    {
        if (x < p.x)
            return true;
        if (x == p.x)
            return y < p.y;
        return false;
    }
};

struct MyEdge
{
    MyEdge(int a, int b) : a(a), b(b) {}
    int a, b;
};

struct GfxStuffer
{
    GfxStuffer() : prg(0) {}

    GLenum type;
    Buf vbo, ibo_delabella, ibo_constrain;
    Buf vbo_voronoi, ibo_voronoi;
    #ifdef WITH_CDT
    Buf ibo_cdt;
    #endif
    MyCoord box[4];

    // for GLdouble only
    GLuint prg; 
    GLint tfm;
    GLint clr;

    void VertexPointer(GLint s, GLenum t, GLsizei d, const GLvoid* ptr)
    {
        if (type == GL_DOUBLE)
            glVertexAttribLPointer(0,s,t,d,ptr);
        else
            glVertexPointer(s,t,d,ptr);
    }

    void EnableVertexArray(bool on)
    {
        if (type == GL_DOUBLE)
            if (on)
                glEnableVertexAttribArray(0);
            else
                glDisableVertexAttribArray(0);
        else
            if (on)
                glEnableClientState(GL_VERTEX_ARRAY);
            else
                glDisableClientState(GL_VERTEX_ARRAY);
    }

    void LoadProj(int vpw, int vph, double cx, double cy, double scale)
    {
        glViewport(0,0,vpw,vph);
        if (type == GL_DOUBLE)
        {
            glUseProgram(prg);

            double left = cx - vpw/scale;
            double right = cx + vpw/scale;
            double bottom = cy - vph/scale;
            double top = cy + vph/scale;

            GLdouble mat[4];

            mat[0] = scale / vpw;
            mat[1] = scale / vph;
            mat[2] = -cx * mat[0];
            mat[3] = -cy * mat[1];

            //glUniformMatrix4dv(tfm, 1, GL_FALSE, mat);
            glUniform4dv(tfm, 1, mat);
        }
        else
        {
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glOrtho(cx - vpw/scale, cx + vpw/scale, cy - vph/scale, cy + vph/scale, -1, +1);

            const static GLfloat z2w[16] =
            {
                1,0,0,0,
                0,1,0,0,
                0,0,0,1,
                0,0,0,0
            };

            glMultMatrixf(z2w);
        }
    }

    void SetColor(float r, float g, float b, float a)
    {
        if (type == GL_DOUBLE)
            glUniform4f(clr, r,g,b,a);
        else
            glColor4f(r,g,b,a);
    }

    void Destroy()
    {
        if (type == GL_DOUBLE)
        {
            glUseProgram(0);
            if (prg)
                glDeleteProgram(prg);
        }

        vbo.Del();
        ibo_delabella.Del();
        ibo_constrain.Del();

        #ifdef WITH_CDT
        ibo_cdt.Del();
        #endif

        #ifdef VORONOI
        vbo_voronoi.Del();
        ibo_voronoi.Del();
        #endif
    }    

    void Upload(
        GLenum gl_e,
        const IDelaBella* idb,
        int points,
        const MyPoint* cloud,
        int constrain_edges,
        const MyEdge* force
        #ifdef VORONOI
        , int voronoi_vertices,
        const MyPoint* voronoi_vtx_buf,
        int voronoi_indices,
        const int* voronoi_idx_buf
        #endif
        #ifdef WITH_CDT
        , const CDT::Triangulation<MyCoord>& cdt,
        const CDT::DuplicatesInfo& dups
        #endif
    )
    {
        Destroy();
        
        type = gl_e;

        if (gl_e == GL_DOUBLE)
        {
            #define CODE(...) #__VA_ARGS__

            static const char* vs_src[] = {CODE(#version 410\n
                uniform dvec4 tfm;
                layout (location = 0) in dvec3 v;
                void main()
                {
                    gl_Position = vec4(
                        tfm.x * v.x + tfm.z * v.z, 
                        tfm.y * v.y + tfm.w * v.z,
                        0.0lf, v.z);
                }
            )};

            static const char* fs_src[] = {CODE(#version 410\n
                uniform vec4 clr;
                layout (location = 0) out vec4 c;
                void main()
                {
                    c = clr;
                }
            )};

            #undef CODE

            char nfolog[1025];
            int nfolen;

            prg = glCreateProgram();
            
            GLuint vs = glCreateShader(GL_VERTEX_SHADER);
            glShaderSource(vs, 1, vs_src, 0);
            glCompileShader(vs);
            glAttachShader(prg,vs);

            glGetShaderInfoLog(vs,1024,&nfolen,nfolog);
            nfolog[nfolen]=0;
            printf("VS:\n%s\n\n",nfolog);

            GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
            glShaderSource(fs, 1, fs_src, 0);
            glCompileShader(fs);
            glAttachShader(prg,fs);

            glGetShaderInfoLog(vs,1024,&nfolen,nfolog);
            nfolog[nfolen]=0;
            printf("FS:\n%s\n\n",nfolog);

            glLinkProgram(prg);

            glDeleteShader(vs);
            glDeleteShader(fs);

            tfm = glGetUniformLocation(prg, "tfm");
            clr = glGetUniformLocation(prg, "clr");
        }

        int gl_s = gl_e == GL_DOUBLE ? sizeof(GLdouble) : sizeof(GLfloat);
        int tris_delabella = idb->GetNumOutputIndices() / 3;
        int contour = idb->GetNumBoundaryVerts();

        if (constrain_edges)
        {
            ibo_constrain.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[2]) * constrain_edges);
            GLuint* map = (GLuint*)ibo_constrain.Map();
            for (int i = 0; i < constrain_edges; i++)
            {
                map[2 * i + 0] = (GLuint)force[i].a;
                map[2 * i + 1] = (GLuint)force[i].b;
            }
            ibo_constrain.Unmap();
        }

        vbo.Gen(GL_ARRAY_BUFFER, gl_s * 3 * points);
        void* vbo_ptr = vbo.Map();

        // let's give a hand to gpu by centering vertices around 0,0
        box[0] = box[2] = cloud[0].x;
        box[1] = box[3] = cloud[0].y;
        for (int i = 0; i<points; i++)
        {
            box[0] = std::min(box[0], cloud[i].x);
            box[1] = std::min(box[1], cloud[i].y);
            box[2] = std::max(box[2], cloud[i].x);
            box[3] = std::max(box[3], cloud[i].y);
        }

        MyCoord vbo_x = (box[0]+box[2]) / 2;
        MyCoord vbo_y = (box[1]+box[3]) / 2;

        #if 0 // vbo centering?
        box[0] -= vbo_x;
        box[1] -= vbo_y;
        box[2] -= vbo_x;
        box[3] -= vbo_y;
        #else
        vbo_x = 0;
        vbo_y = 0;
        #endif

        if (gl_e == GL_DOUBLE)
        {
            GLdouble* p = (GLdouble*)vbo_ptr;
            for (int i = 0; i<points; i++)
            {
                p[3*i+0] = (GLdouble)(cloud[i].x - vbo_x);
                p[3*i+1] = (GLdouble)(cloud[i].y - vbo_y);
                p[3*i+2] = (GLdouble)(1.0); // i%5; // color
            }
        }
        else
        {
            GLfloat* p = (GLfloat*)vbo_ptr;
            for (int i = 0; i<points; i++)
            {
                p[3*i+0] = (GLfloat)(cloud[i].x - vbo_x);
                p[3*i+1] = (GLfloat)(cloud[i].y - vbo_y);
                p[3*i+2] = (GLfloat)(1.0); // i%5; // color
            }
        }
        vbo.Unmap();

        ibo_delabella.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delabella + sizeof(GLuint) * contour);	
        GLuint* ibo_ptr = (GLuint*)ibo_delabella.Map();
        const DelaBella_Triangle* dela = idb->GetFirstDelaunaySimplex();
        for (int i = 0; i<tris_delabella; i++)
        {
            int v0 = dela->v[0]->i;
            int v1 = dela->v[1]->i;
            int v2 = dela->v[2]->i;

            ibo_ptr[3*i+0] = (GLuint)v0;
            ibo_ptr[3*i+1] = (GLuint)v1;
            ibo_ptr[3*i+2] = (GLuint)v2;

            dela = dela->next;
        }

        const DelaBella_Vertex* vert = idb->GetFirstBoundaryVertex();
        for (int i = 0; i<contour; i++)    
        {
            ibo_ptr[i + 3*tris_delabella] = (GLuint)vert->i;
            vert = vert->next;
        }

        ibo_delabella.Unmap();

        #ifdef VORONOI
        vbo_voronoi.Gen(GL_ARRAY_BUFFER, gl_s * 3 * voronoi_vertices);
        ibo_voronoi.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)* voronoi_indices);
        void* vbo_voronoi_ptr = vbo_voronoi.Map();
        GLuint* ibo_voronoi_ptr = (GLuint*)ibo_voronoi.Map();

        if (gl_e == GL_DOUBLE)
        {
            GLdouble* p = (GLdouble*)vbo_voronoi_ptr;
            for (int i = 0; i < voronoi_vertices; i++)
            {
                if (i < voronoi_vertices - contour)
                {
                    p[3 * i + 0] = (GLdouble)(voronoi_vtx_buf[i].x - vbo_x);
                    p[3 * i + 1] = (GLdouble)(voronoi_vtx_buf[i].y - vbo_y);
                    p[3 * i + 2] = (GLdouble)1;
                }
                else
                {
                    p[3 * i + 0] = (GLdouble)voronoi_vtx_buf[i].x;
                    p[3 * i + 1] = (GLdouble)voronoi_vtx_buf[i].y;
                    p[3 * i + 2] = (GLdouble)0;
                }
            }
        }
        else
        {
            GLfloat* p = (GLfloat*)vbo_voronoi_ptr;
            for (int i = 0; i < voronoi_vertices; i++)
            {
                if (i < voronoi_vertices - contour)
                {
                    p[3 * i + 0] = (GLfloat)(voronoi_vtx_buf[i].x - vbo_x);
                    p[3 * i + 1] = (GLfloat)(voronoi_vtx_buf[i].y - vbo_y);
                    p[3 * i + 2] = (GLfloat)1;
                }
                else
                {
                    p[3 * i + 0] = (GLfloat)voronoi_vtx_buf[i].x;
                    p[3 * i + 1] = (GLfloat)voronoi_vtx_buf[i].y;
                    p[3 * i + 2] = (GLfloat)0;
                }
            }
        }

        for (int i = 0; i < voronoi_indices; i++)
            ibo_voronoi_ptr[i] = (GLuint)(voronoi_idx_buf[i]);

        vbo_voronoi.Unmap();
        ibo_voronoi.Unmap();
        #endif

        #ifdef WITH_CDT
        int tris_cdt = (int)cdt.triangles.size();
        ibo_cdt.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_cdt);
        {
            ibo_ptr = (GLuint*)ibo_cdt.Map();

            if (dups.duplicates.size())
            {
                // let's make inverse mapping first
                int* invmap = 0;
                int invmap_size = (int)dups.mapping.size() - (int)dups.duplicates.size();
                invmap = (int*)malloc(sizeof(int) * invmap_size);
                for (int i = 0; i < (int)dups.mapping.size(); i++)
                    invmap[dups.mapping[i]] = i;

                for (int i = 0; i < tris_cdt; i++)
                {
                    int a = cdt.triangles[i].vertices[0];
                    int b = cdt.triangles[i].vertices[1];
                    int c = cdt.triangles[i].vertices[2];
                    ibo_ptr[3 * i + 0] = (GLuint)invmap[c];
                    ibo_ptr[3 * i + 1] = (GLuint)invmap[b];
                    ibo_ptr[3 * i + 2] = (GLuint)invmap[a];
                }

                free(invmap);
            }
            else
            {
                // 1:1 mapping (no dups)
                for (int i = 0; i < tris_cdt; i++)
                {
                    int a = cdt.triangles[i].vertices[0];
                    int b = cdt.triangles[i].vertices[1];
                    int c = cdt.triangles[i].vertices[2];
                    ibo_ptr[3 * i + 0] = (GLuint)c;
                    ibo_ptr[3 * i + 1] = (GLuint)b;
                    ibo_ptr[3 * i + 2] = (GLuint)a;
                }
            }
            ibo_cdt.Unmap();
        }
        #endif

    }
};

int main(int argc, char* argv[])
{
	#ifdef _WIN32
	SetProcessDPIAware();
	#endif

	if (argc<2)
	{
        printf("usage: %s <input.xy> [output.abc]\n", argc > 0 ? argv[0] : "<executable>");
		printf("required argument (.xy file name / number of random points) missing, terminating!\n");
		return -1;
	}

	std::vector<MyPoint> cloud;
    std::vector<MyEdge> force;

	FILE* f = fopen(argv[1],"r");
    int n = atoi(argv[1]);
    if (!f)
    {
        if (n <= 2)
        {
            printf("can't open %s file, terminating!\n", argv[1]);
            return -1;
        }
        printf("generating random %d points\n", n);
        std::random_device rd{};
        std::mt19937_64 gen{ 0x12345678 /*rd()*/};

        //std::uniform_real_distribution<double> d(-2.503515625, +2.503515625);
        //std::normal_distribution<double> d{0.0,2.0};
        std::gamma_distribution<double> d(0.1,2.0);

        MyCoord max_coord = sizeof(MyCoord) < 8 ? /*float*/0x1.p31 : /*double*/0x1.p255;

        for (int i = 0; i < n; i++)
        {
            MyPoint p = { (d(gen) + 50.0), (d(gen) + 50.0) };
            //MyPoint p = { d(gen), d(gen) };
            
            //p.x *= 0x1.p250;
            //p.y *= 0x1.p250;

            // please, leave some headroom for arithmetics!
            assert(std::abs(p.x) <= max_coord && std::abs(p.y) <= max_coord);

            cloud.push_back(p);
        }
        
        
        /*
        {
            const double x = 0x5af2efc1.p-30;
            const double y = 0x348268e0.p-30;
            const double r = 0x6904d1c1.p-30;
            const MyPoint p[4] = 
            { 
                {0,0}, 
                {0,(MyCoord)(2*y)}, 
                {(MyCoord)x,(MyCoord)(r+y)}, 
                {(MyCoord)x,(MyCoord)(r+3*y)} 
            };

            const MyCoord dx = (MyCoord)(2*x);
            const MyCoord dy = (MyCoord)(2*(r+y));

            int rows = (int)ceil(sqrt(n / 7.0));
            int cols = (n + 2 * rows) / (4*rows);

            n = rows * cols * 4;

            printf("%d\n", rows* cols * 4);

            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < cols; col++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        MyPoint q = p[i];
                        q.x += col * dx;
                        q.y += row * dy;
                        cloud.push_back(q);
                    }
                }
            }
        }
        */
        
        if (1)
        {
            int m = n / 10;

            // init sub[] with all n point indices
            int* sub = (int*)malloc(sizeof(int) * n);
            for (int i = 0; i < n; i++)
                sub[i] = i;

            // pick m random ones from n
            // place them as first m items of sub[]
            for (int i = 0; i < m; i++)
            {
                int j = i + gen() % ((size_t)n - i);
                int r = sub[j];
                sub[j] = sub[i];
                sub[i] = r;
            }

            std::vector<MyPoint> xxx;
            for (int i = 0; i < m; i++)
                xxx.push_back(cloud[sub[i]]);

            IDelaBella* helper = IDelaBella::Create();
            helper->Triangulate(m, &xxx.data()->x, &xxx.data()->y, sizeof(MyPoint));

            // to avoid forth-and-back edges repeatitions,
            // traverse all faces but use edges with 
            // ascending y or in case of flat y use only if ascending x

            const IDelaBella2<MyCoord>::Simplex* dela = helper->GetFirstDelaunaySimplex();
            while (dela)
            {
                for (int a = 0, b = 1, c = 2; a < 3; b = c, c = a, a++)
                {
                    //if (dela->f[c]->index < 0 || dela->v[a]->i < dela->v[b]->i)
                    if (dela->f[c]->index >= 0 && dela->v[a]->i < dela->v[b]->i)
                        force.push_back(MyEdge(sub[dela->v[a]->i], sub[dela->v[b]->i]));
                }

                dela = dela->next;
            }

            helper->Destroy();
            free(sub);
        }


        /*
        int a = gen() % n;
        for (int i = 0; i < 1000; i++)
        {
            int b = (a + gen() % (n-1)) % n;
            force.push_back(MyEdge(a, b));
        }
        */
     
        /*
        for (int i = 0; i < n; i++)
        {
            MyPoint p = { d(gen), 0 };
            //p.y = p.x;
            cloud.push_back(p);
        }
        MyPoint p = { 0, 1};
        cloud.push_back(p);
        */
        
	}
    else
    {
        int r,n,c;
        r = fscanf(f, "%d %d", &n, &c);

        for (int i=0; i<n; i++)
        {
            double dbl_x, dbl_y;
            // allow variety of separators and extra fields till end of the line
            int n = fscanf(f,"%lf%*[,; \v\t]%lf%*[^\n]", &dbl_x, &dbl_y);

            MyCoord x = (MyCoord)dbl_x;
            MyCoord y = (MyCoord)dbl_y;

            MyPoint p = {x,y};
            cloud.push_back(p);
        }

        for (int i = 0; i < c; i++)
        {
            int a, b;
            r = fscanf(f, "%d %d", &a, &b);
            MyEdge e = {a, b};
            force.push_back(e);
        }

        fclose(f);
    }

    int points = (int)cloud.size();

    #ifdef WITH_CDT

        std::vector<CDT::V2d<MyCoord>> nodups;
        for (size_t i = 0; i < cloud.size(); i++)
        {
            CDT::V2d<double> v;
            v.x = cloud[i].x;
            v.y = cloud[i].y;
            nodups.push_back(v);
        }

        std::vector<CDT::Edge> edges;
        for (size_t c = 0; c < force.size(); c++)
        {
            CDT::Edge e
            { 
                (CDT::VertInd)force[c].a,
                (CDT::VertInd)force[c].b
            };
            edges.push_back(e);
        }

        printf("cdt removing dups... ");
        uint64_t t0 = uSec();

        CDT::DuplicatesInfo dups = CDT::RemoveDuplicatesAndRemapEdges(nodups,edges);

        uint64_t t1 = uSec();
        printf("%d ms\n", (int)((t1 - t0) / 1000));

        printf("cdt triangulation... ");
        CDT::Triangulation<MyCoord> cdt;
        cdt.insertVertices(nodups);

        uint64_t t2 = uSec();
        printf("%d ms\n", (int)((t2 - t1) / 1000));

        if (force.size()>0)
        {
            printf("cdt forcing edges... ");
            cdt.insertEdges(edges);
            uint64_t t3 = uSec();
            printf("%d ms\n", (int)((t3 - t2) / 1000));
        }

        uint64_t t4 = uSec();
        printf("cdt erasing super... ");
        cdt.eraseSuperTriangle();
        uint64_t t5 = uSec();
        printf("%d ms\n", (int)((t5 - t4) / 1000));
        printf("CDT TOTAL: %d\n", (int)((t5 - t0) / 1000));

        printf("CDT triangles = %d\n", (int)cdt.triangles.size());

        uint64_t t_p0 = uSec();
        std::vector<CDT::Poly> cdt_polys = CDT::Polygonize(cdt);
        uint64_t t_p1 = uSec();

        printf("CDT POLYS = %d (in %d ms)\n", (int)cdt_polys.size(), (int)((t_p1 - t_p0) / 1000));

    #endif


    #ifdef WITH_DELAUNATOR
    {
        std::vector<double> coords;
        for (int i=0; i<points; i++)
        {
            coords.push_back(cloud[i].x);
            coords.push_back(cloud[i].y);
        }
        uint64_t t0 = uSec();
        printf("running delaunator...\n");
        delaunator::Delaunator* d = 0;
        try
        {
            d = new delaunator::Delaunator(coords);
        }
        catch (...)
        {
            printf("delaunator threw an exception!\n");
            d = 0;
        }
        int tris_delaunator = d ? (int)d->triangles.size() / 3 : 0;
        uint64_t t1 = uSec();
        printf("elapsed %d ms\n", (int)((t1-t0)/1000));
        printf("delaunator triangles: %d\n", tris_delaunator);
    }
    #endif

	IDelaBella* idb = IDelaBella::Create();
	idb->SetErrLog(errlog, stdout);
	
    printf("running delabella...\n");
    uint64_t t6 = uSec();
	int verts = idb->Triangulate(points, &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
	int tris_delabella = verts > 0 ? verts / 3 : 0;
    int contour = idb->GetNumBoundaryVerts();
    int non_contour = idb->GetNumInternalVerts();
	int vert_num = contour + non_contour;

    #ifdef VORONOI
    printf("Polygonizing for VD\n");
    //idb->Polygonize(); // optional

    printf("Generating VD vertices\n");
    int voronoi_vertices = idb->GenVoronoiDiagramVerts(0, 0, 0);
    MyPoint* voronoi_vtx_buf = (MyPoint*)malloc(voronoi_vertices * sizeof(MyPoint));
    idb->GenVoronoiDiagramVerts(&voronoi_vtx_buf->x, &voronoi_vtx_buf->y, sizeof(MyPoint));

    printf("Generating VD indices\n");
    #ifdef VORONOI_POLYS
    // testing... will remove
    int voronoi_closed_indices;
    int voronoi_indices = idb->GenVoronoiDiagramPolys(0, 0, 0);
    int* voronoi_idx_buf = (int*)malloc(voronoi_indices * sizeof(int));
    idb->GenVoronoiDiagramPolys(voronoi_idx_buf, sizeof(int), &voronoi_closed_indices);
    #else
    int voronoi_closed_indices = 0;
    int voronoi_indices = idb->GenVoronoiDiagramEdges(0, 0);
    int* voronoi_idx_buf = (int*)malloc(voronoi_indices * sizeof(int));
    idb->GenVoronoiDiagramEdges(voronoi_idx_buf, sizeof(int));
    #endif

    printf("VD vertices = %d, indices = %d\n", voronoi_vertices, voronoi_indices);
    #endif

    if (force.size()>0)
        idb->ConstrainEdges((int)force.size(), &force.data()->a, &force.data()->b, (int)sizeof(MyEdge), false);

    const DelaBella_Triangle** dela_polys = (const DelaBella_Triangle**)malloc(sizeof(const DelaBella_Triangle*) * tris_delabella);
    int polys_delabella = idb->Polygonize(dela_polys);

    printf("delabella triangles: %d\n", tris_delabella);
    printf("delabella contour: %d\n", contour);
    printf("delabella polygons: %d\n", polys_delabella);

    //return 0;

	// if positive, all ok 
	if (verts<=0)
	{
        printf("nothing interesting to show, exiting!\n");
		// no points given or all points are colinear
		// make emergency call ...
		idb->Destroy();
		return -2;
	}

    // COMPARE CDT with Dela
    // if (0)
    {
        #ifdef WITH_CDT

        int tris_cdt = (int)cdt.triangles.size();

        if (tris_delabella != tris_cdt || polys_delabella == cdt_polys.size())
            printf("WARNING! Results are not comparable - different number of tris or polys\n");
        else
        {
            // note:
            // currently, any difference between constrained edges passing 
            // internal edge of a polygon won't be detected !!!

            int poly_indices = 2 * polys_delabella + tris_delabella;

            struct MyPoly
            {
                int size;
                int offs;
            };

            struct PolyPred
            {
                PolyPred(const std::vector<MyPoint>& v) : v(v) {}
                bool operator () (const MyPoly& p, const MyPoly& q) const
                {
                    if (p.size < q.size)
                        return true;
                    if (p.size == q.size)
                    {
                        for (size_t i = 0; i < (size_t)p.size; i++)
                        {
                            if (v[i + p.offs] == v[i + q.offs])
                                continue;
                            return v[i + p.offs] < v[i + q.offs];
                        }
                    }
                    return false;
                }
                const std::vector<MyPoint>& v;
            };

            printf("preping cdt for cmp ...\n");
            std::vector<MyPoint> cdt_v(poly_indices);
            std::vector<MyPoly> cdt_p(cdt_polys.size());
            for (int p = 0, n = 0; p < (int)cdt_polys.size(); p++)
            {
                int s = 0;

                for (int i = 0; i < 3; i++)
                {
                    int j = n + s;
                    int k = cdt.triangles[cdt_polys[p][0]].vertices[i];
                    cdt_v[j].x = cdt.vertices[k].x;
                    cdt_v[j].y = cdt.vertices[k].y;
                    s++;
                }

                for (int i = 1; i < (int)cdt_polys[p].size(); i++)
                {
                    int j = n + s;
                    int k = cdt.triangles[cdt_polys[p][i]].vertices[0];
                    cdt_v[j].x = cdt.vertices[k].x;
                    cdt_v[j].y = cdt.vertices[k].y;
                    s++;
                }

                // convert from ccw to cw
                std::reverse(cdt_v.begin() + n, cdt_v.begin() + n + s);

                // find smallest point, make it leftmost
                std::vector<MyPoint>::iterator smallest = std::min_element(cdt_v.begin() + n, cdt_v.begin() + n + s);
                std::rotate(cdt_v.begin() + n, smallest, cdt_v.begin() + n + s);


                cdt_p[p].size = s;
                cdt_p[p].offs = n;

                n += s;
            }
            printf("sorting cdt ...\n");
            std::sort(cdt_p.begin(), cdt_p.end(), PolyPred(cdt_v));

            printf("preping idb for cmp ...\n");
            std::vector<MyPoint> idb_v(poly_indices);
            std::vector<MyPoly> idb_p(polys_delabella);
            for (int p = 0, n = 0; p < polys_delabella; p++)
            {
                int s = 0;

                const DelaBella_Triangle* f = dela_polys[p];

                for (int i = 0; i < 3; i++)
                {
                    int j = n + s;
                    const DelaBella_Vertex* k = f->v[i];
                    idb_v[j].x = k->x;
                    idb_v[j].y = k->y;
                    s++;
                }

                f = f->next;

                while (f && f->index == p)
                {
                    int j = n + s;
                    const DelaBella_Vertex* k = f->v[0];
                    idb_v[j].x = k->x;
                    idb_v[j].y = k->y;
                    s++;
                    f = f->next;
                }

                // find smallest point, rotate it to the left
                std::vector<MyPoint>::iterator smallest = std::min_element(idb_v.begin() + n, idb_v.begin() + n + s);
                std::rotate(idb_v.begin() + n, smallest, idb_v.begin() + n + s);

                idb_p[p].size = s;
                idb_p[p].offs = n;

                n += s;
            }
            printf("sorting idb ...\n");
            std::sort(idb_p.begin(), idb_p.end(), PolyPred(idb_v));

            printf("COMPARING... ");
            bool compare_ok = true;
            for (int p = 0; p < polys_delabella; p++)
            {
                MyPoly p_idb = idb_p[p];
                MyPoly p_cdt = cdt_p[p];

                MyPoint* data_idb = idb_v.data() + p_idb.offs;
                MyPoint* data_cdt = cdt_v.data() + p_cdt.offs;

                if (p_idb.size != p_cdt.size ||
                    memcmp(data_idb, data_cdt, sizeof(MyPoint) * p_idb.size))
                {
                    compare_ok = false;
                    break;
                }
            }
            printf(compare_ok ? "OK\n" : "DIFFERENT!\n");
        }
        #endif
    }

    if (argc>=3)
    {
        f = fopen(argv[2],"w");
        if (f)
        {
            for (int i=0; i<tris_delabella; i++)
            {
                const DelaBella_Triangle* dela = idb->GetFirstDelaunaySimplex();
                fprintf(f,"%d %d %d\n", 
                    dela->v[0]->i,
                    dela->v[1]->i,
                    dela->v[2]->i);
                dela = dela->next;
            }
            fclose(f);
        }
    }

    SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1");
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    //SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	// create viewer wnd
    int width = 800, height = 600;
    const char* title = "delablella-sdl2";
    SDL_Window * window = SDL_CreateWindow( title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_ALLOW_HIGHDPI | SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (!window)
    {
        printf("SDL_CreateWindow failed, terminating!\n");
        idb->Destroy();
        return -1;
    }

    SDL_GLContext context = SDL_GL_CreateContext( window );
    if (!context)
    {
        printf("SDL_GL_CreateContext failed, terminating!\n");
        idb->Destroy();
        return -1;
    }

    bool prim_restart = false;
    {
        const char* ext = (const char*)glGetString(GL_EXTENSIONS);
        while (ext)
        {
            ext = strstr(ext,"GL_NV_primitive_restart");
            if (ext)
            {
                if (ext[23]==0 || ext[23]==' ')
                {
                    prim_restart = true;
                    break;
                }
                ext += 23;
            }
        }
    }

    // check if we can use shaders with dmat/dvec/double support
    // it is standard in GLSL 4.0 and above
    int glsl_ver = 0;
    const char* glsl_str = (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);
    if (glsl_str && glGetError() == GL_NO_ERROR)
    {
        int v[2];
        int from,to;
        if (2 == sscanf(glsl_str, "%d.%n%d%n", v+0,&from,v+1,&to))
        {
            int num = to-from;
            if (num == 1)
                v[1] *= 10;
            else
            while (num>2)
            {
                v[1] /= 10;
                num--;
            }
            glsl_ver = v[0]*100+v[1];
        }
    }

    //glsl_ver = 0;

	if (!BindGL(prim_restart, glsl_ver >= 410))
	{
		printf("Can't bind to necessary GL functions, terminating!\n");
		idb->Destroy();
		return -1;
	}

	printf("preparing graphics...\n");

    GfxStuffer gfx;

    gfx.Upload( glsl_ver >= 410 ? GL_DOUBLE : GL_FLOAT,
                idb,
                (int)cloud.size(),
                cloud.data(),
                (int)force.size(),
                force.data()
                #ifdef VORONOI
                , voronoi_vertices,
                voronoi_vtx_buf,
                voronoi_indices,
                voronoi_idx_buf
                #endif
                #ifdef WITH_CDT
                , cdt,
                dups
                #endif
                );

    #ifdef VORONOI
    free(voronoi_idx_buf);
    free(voronoi_vtx_buf);
    #endif

    int vpw, vph;
    SDL_GL_GetDrawableSize(window, &vpw, &vph);

    double cx = 0.5 * (gfx.box[0]+gfx.box[2]);
    double cy = 0.5 * (gfx.box[1]+gfx.box[3]);
    double scale = 2.0 * fmin((double)vpw/(gfx.box[2]-gfx.box[0]),(double)vph/(gfx.box[3]-gfx.box[1]));
    int zoom = -3+(int)round(log(scale) / log(1.01));

    int drag_x, drag_y, drag_zoom;
    double drag_cx, drag_cy;
    int drag = 0;

    if (prim_restart)
    {
        glPrimitiveRestartIndex((GLuint)~0);
        glEnable(GL_PRIMITIVE_RESTART);
    }

	printf("going interactive.\n");

    bool show_f = true; // fill
    bool show_b = true; // boundary
    bool show_v = true; // voronoi
    bool show_c = true; // constraints
    bool show_x = true; // cross-compare with cdt
    bool show_d = true; // delaunay

    printf("\n");
    printf(
        "Change layers visibility while graphics window is in focus:\n"
        "[F]ill, [B]oundary, "
        #ifdef VORONOI
        "[V]oronoi, "
        #endif
        #ifdef WITH_CDT
        "[C]onstraints, [D]elaunay, [X]compare\n"
        #else
        "[C]onstraints, [D]elaunay\n"
        #endif
        );

    printf("\n");
    printf(
        "Mouse controls:\n"
        "[LMB]pan, [RMB/wheel]zoom\n");

    printf("\n");

    for( ;; )
    {
        bool x = false;
        SDL_Event event;
        while( SDL_PollEvent( &event ) )
        {
            switch( event.type )
            {
                case SDL_MOUSEWHEEL:
                {
                    if (drag==0)
                    {
                        scale = pow(1.01, zoom);
                        int vpw, vph;
                        SDL_GL_GetDrawableSize(window, &vpw, &vph);

                        SDL_GetMouseState(&drag_x,&drag_y);

                        drag_cx = cx + 2.0 * (drag_x - vpw*0.5) / scale;
                        drag_cy = cy + 2.0 * (vph*0.5 - drag_y) / scale;

                        zoom += event.wheel.y * 10;
                        scale = pow(1.01, zoom);

                        cx = (drag_cx - 2.0 * (drag_x - vpw*0.5) / scale);
                        cy = (drag_cy - 2.0 * (vph*0.5 - drag_y) / scale);
                    }

                    break;
                }

                case SDL_MOUSEMOTION:
                {
                    if (drag == 1)
                    {
                        int dx = event.motion.x - drag_x;
                        int dy = event.motion.y - drag_y;

                        double scale = pow(1.01, zoom);
                        cx = drag_cx - 2*dx/scale;
                        cy = drag_cy + 2*dy/scale;
                    }

                    if (drag == 2)
                    {
                        int dx = event.motion.x - drag_x;
                        int dy = event.motion.y - drag_y;

                        zoom = drag_zoom - dy;

                        scale = pow(1.01, zoom);
                        int vpw, vph;
                        SDL_GL_GetDrawableSize(window, &vpw, &vph);

                        cx = (drag_cx - 2.0 * (drag_x - vpw*0.5) / scale);
                        cy = (drag_cy - 2.0 * (vph*0.5 - drag_y) / scale);                        
                    }
                    break;
                }

                case SDL_MOUSEBUTTONDOWN:
                {
                    if (!drag)
                    {
                        if (event.button.button == SDL_BUTTON_LEFT)
                        {
                            drag = 1;
                            drag_x = event.button.x;
                            drag_y = event.button.y;
                            drag_cx = cx;
                            drag_cy = cy;
                            drag_zoom = zoom;
                            SDL_CaptureMouse(SDL_TRUE);
                        }
                        if (event.button.button == SDL_BUTTON_RIGHT)
                        {
                            scale = pow(1.01, zoom);
                            int vpw, vph;
                            SDL_GL_GetDrawableSize(window, &vpw, &vph);
                            double dx = cx + 2.0 * (event.button.x - vpw*0.5) / scale;
                            double dy = cy + 2.0 * (vph*0.5 - event.button.y) / scale;

                            drag = 2;
                            drag_x = event.button.x;
                            drag_y = event.button.y;
                            drag_cx = dx;
                            drag_cy = dy;
                            drag_zoom = zoom;
                            SDL_CaptureMouse(SDL_TRUE);
                        }
                    }
                    break;
                }

                case SDL_MOUSEBUTTONUP:
                {
                    if ((event.button.button == SDL_BUTTON_LEFT && drag == 1) ||
                        (event.button.button == SDL_BUTTON_RIGHT && drag == 2))
                    {
                        drag = 0;
                        SDL_CaptureMouse(SDL_FALSE);
                    }
                    break;
                }

                case SDL_WINDOWEVENT:
                {
                    switch (event.window.event) {

                        case SDL_WINDOWEVENT_CLOSE:   // exit game
                            x = true;
                            break;

                        default:
                            break;
                    }
                    break;
                }

                case SDL_KEYDOWN:
                    if (event.key.keysym.sym == SDLK_f)
                        show_f = !show_f;
                    if (event.key.keysym.sym == SDLK_b)
                        show_b = !show_b;
                    if (event.key.keysym.sym == SDLK_v)
                        show_v = !show_v;
                    if (event.key.keysym.sym == SDLK_c)
                        show_c = !show_c;
                    if (event.key.keysym.sym == SDLK_x)
                        show_x = !show_x;
                    if (event.key.keysym.sym == SDLK_d)
                        show_d = !show_d;
                    break;
            }

            if (x)
                break;
        }

        if (x)
            break;

        int vpw, vph;
        SDL_GL_GetDrawableSize(window, &vpw, &vph);

        double scale = pow(1.01, zoom);

        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);

        gfx.LoadProj(vpw,vph, cx,cy, scale);

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);

        gfx.vbo.Bind();
		gfx.VertexPointer(3, gfx.type, 0, 0);
		gfx.EnableVertexArray(true);

        gfx.ibo_delabella.Bind();

        // grey fill
        if (show_f)
        {
            gfx.SetColor(0.2f, 0.2f, 0.2f, 1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delabella * 3, GL_UNSIGNED_INT, 0);
        }

        // paint constraints
        int constrain_indices = 2*(int)force.size();
        if (constrain_indices && show_c)
        {
            gfx.ibo_constrain.Bind();
            glLineWidth(3.0f);
            gfx.SetColor(.9f, .9f, .9f, 1.0f);
            glDrawElements(GL_LINES, constrain_indices, GL_UNSIGNED_INT, (GLuint*)0);
            glLineWidth(1.0f);

            // oops
            gfx.ibo_delabella.Bind();
        }

        if (show_d)
        {
            gfx.SetColor(1.0f, 0.0f, 0.0f, 1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_TRIANGLES, tris_delabella * 3, GL_UNSIGNED_INT, 0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (show_b)
        {
            gfx.SetColor(0.0f, 0.0f, 1.0f, 1.0f);
            glLineWidth(3.0f);
            glDrawElements(GL_LINE_LOOP, contour, GL_UNSIGNED_INT, (GLuint*)0 + tris_delabella * 3);
            glLineWidth(1.0f);
        }

        // compare with CDT
        #ifdef WITH_CDT
        if (show_x)
        {
            int tris_cdt = (int)cdt.triangles.size();

            glEnable(GL_BLEND);
            glBlendFunc(GL_ONE, GL_ONE);
            gfx.ibo_cdt.Bind();
            gfx.SetColor(0.0f,0.0f,1.0f,1.0f);
            glLineWidth(1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_cdt * 3, GL_UNSIGNED_INT, 0);
            glDisable(GL_BLEND);
        }
        #endif

        // voronoi!
        #ifdef VORONOI
        if (show_v)
        {
            gfx.vbo_voronoi.Bind();
            gfx.VertexPointer(3, gfx.type, 0, 0);
            gfx.EnableVertexArray(true);

            // voro-verts in back
            gfx.SetColor(1.0f, 1.0f, 0.0f, 1.0f);
            glPointSize(3.0f);
            glDrawArrays(GL_POINTS, 0, voronoi_vertices - contour);
            glPointSize(1.0f);
            gfx.ibo_voronoi.Bind();

            gfx.SetColor(0.0f, 0.75f, 0.0f, 1.0f);

            #ifdef VORONOI_POLYS
            // draw structured polys, note: open polys are silently closed (at infinity)
            // glDrawElements(GL_LINE_LOOP, voronoi_indices, GL_UNSIGNED_INT, (GLuint*)0);
            // if you wanna be a ganan: after first M closed polys switch from line_loops to line_strips
            // and draw remaining N open polygons
            glDrawElements(GL_LINE_LOOP, voronoi_closed_indices, GL_UNSIGNED_INT, (GLuint*)0);
            glDrawElements(GL_LINE_STRIP, voronoi_indices - voronoi_closed_indices, GL_UNSIGNED_INT, (GLuint*)0+voronoi_closed_indices);
            #else
            // draw edge soup
            glDrawElements(GL_LINES, voronoi_indices, GL_UNSIGNED_INT, (GLuint*)0);
            #endif
        }
        #endif


        gfx.vbo.Bind();
        gfx.ibo_delabella.Bind();
        gfx.VertexPointer(3, gfx.type, 0, 0);
        gfx.EnableVertexArray(true);

        // put verts over everything else
        gfx.SetColor(1.0f, 1.0f, 0.0f, 1.0f);
        glPointSize(3.0f);
        glDrawArrays(GL_POINTS, 0, points);
        glPointSize(1.0f);

        SDL_GL_SwapWindow(window);
        SDL_Delay(15);
    }

    gfx.Destroy();

    SDL_GL_DeleteContext( context );
    SDL_DestroyWindow( window );
    SDL_Quit();

	printf("exiting.\n");

	return 0;
}

