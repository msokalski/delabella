/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#define _CRT_SECURE_NO_WARNINGS

#define BENCH

//#define CULLING

//#define VORONOI
//#define VORONOI_POLYS
// otherwise EDGES

// override build define
#undef WITH_DELAUNATOR 
//#define WITH_DELAUNATOR

// override build define
#undef WITH_CDT
#define WITH_CDT

// override build define
#undef WITH_FADE 
//#define WITH_FADE


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


#ifdef WITH_FADE
#ifdef _WIN32
#pragma comment(lib,"libgmp-10.lib")
#ifdef _DEBUG
#pragma comment(lib,"fade2D_x64_v142_Debug.lib")
#else
#pragma comment(lib,"fade2D_x64_v142_Release.lib")
#endif
#endif
#include "fade/include_fade2d/Fade_2D.h"
using namespace GEOM_FADE2D;
#endif

// competitors
#ifdef WITH_DELAUNATOR
#include "delaunator-cpp/include/delaunator-header-only.hpp"
#endif

#ifdef WITH_CDT
#include "CDT/CDT/include/CDT.h"

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

        const auto& v1 = cdt.vertices[tri.vertices[0]];
        const auto& v2 = cdt.vertices[tri.vertices[1]];
        const auto& v3 = cdt.vertices[tri.vertices[2]];

        // compare i'th tri with its adjacent tris
        int e_from = 0, e_to = 1;
        for (const auto adj : tri.neighbors)
        {
            // but only if there's adjacent tri and it is processed already ...
            if (adj == noNeighbor || adj > iT)
            {
                e_from = e_to;
                e_to = e_to == 2 ? 0 : e_to + 1;
                continue;
            }

            // ... and edge between these 2 faces is not marked as fixed!
            auto fix_it = cdt.fixedEdges.find(Edge(tri.vertices[e_from], tri.vertices[e_to]));

            e_from = e_to;
            e_to = e_to == 2 ? 0 : e_to + 1;

            if (fix_it != cdt.fixedEdges.end())
                continue;

            // locate reflex vert in adj triangle
            const auto& vr = cdt.vertices[opposedVertex(cdt.triangles[adj], iT)];

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

#else
#include "predicates.h" // we use them for panning expansions
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
PFNGLUNIFORM1IPROC glUniform1i = 0;
PFNGLUNIFORM4FPROC glUniform4f = 0;
//PFNGLUNIFORMMATRIX4DVPROC glUniformMatrix4dv = 0;
PFNGLUNIFORM4DVPROC glUniform4dv = 0;
PFNGLUNIFORM4FVPROC glUniform4fv = 0;
PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer = 0;
PFNGLVERTEXATTRIBLPOINTERPROC glVertexAttribLPointer = 0;
PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray = 0;
PFNGLDISABLEVERTEXATTRIBARRAYPROC glDisableVertexAttribArray = 0;

PFNGLGENVERTEXARRAYSPROC glGenVertexArrays = 0;
PFNGLDELETEVERTEXARRAYSPROC glDeleteVertexArrays = 0;
PFNGLBINDVERTEXARRAYPROC glBindVertexArray = 0;

PFNGLTEXBUFFERPROC glTexBuffer = 0;

bool BindGL()
{
	#define BINDGL(proc) if ((*(void**)&proc = SDL_GL_GetProcAddress(#proc)) == 0) return false;
    
	BINDGL(glBindBuffer);
	BINDGL(glDeleteBuffers);
	BINDGL(glGenBuffers);
	BINDGL(glBufferData);
	BINDGL(glBufferSubData);
	BINDGL(glMapBuffer);
	BINDGL(glUnmapBuffer);

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
    //BINDGL(glUniformMatrix4dv); // packed proj to vec4 
    BINDGL(glUniform4dv);
    BINDGL(glUniform4fv);
    BINDGL(glUniform4f);
    BINDGL(glUniform1i);
    BINDGL(glVertexAttribPointer);
    BINDGL(glVertexAttribLPointer);
    BINDGL(glEnableVertexAttribArray);
    BINDGL(glDisableVertexAttribArray);
    BINDGL(glGetUniformLocation);

    BINDGL(glGenVertexArrays);
    BINDGL(glDeleteVertexArrays);
    BINDGL(glBindVertexArray);

    BINDGL(glPrimitiveRestartIndex);

    BINDGL(glTexBuffer);

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
    GLsizeiptr size;
    void* map;
    bool mapped;

    Buf() : buf(0),target(0),size(0),map(0),mapped(false) {}

    GLuint Gen(GLenum t, GLsizeiptr s)
    {
        target = t;
        size = s;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_TEXTURE_BUFFER)
            glGetIntegerv(GL_TEXTURE_BUFFER_BINDING, &push);

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
        else
        if (target == GL_TEXTURE_BUFFER)
            glGetIntegerv(GL_TEXTURE_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        map = glMapBuffer(target, GL_WRITE_ONLY);
        mapped = true;
        if (!map)
        {
            map = malloc(size);
            mapped = false;
        }

        assert(map);
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
        else
        if (target == GL_TEXTURE_BUFFER)
            glGetIntegerv(GL_TEXTURE_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        if (mapped)
        {
            GLboolean buf_unmap_ok = glUnmapBuffer(target);
            assert(buf_unmap_ok);
        }
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


typedef double MyCoord;
//typedef intptr_t MyIndex;
//#define IDXF "%zd"
typedef int32_t MyIndex;
#define IDXF "%d"

typedef IDelaBella2<MyCoord, MyIndex> IDelaBella;
typedef IDelaBella::Vertex DelaBella_Vertex;
typedef IDelaBella::Simplex DelaBella_Triangle;

struct MyPoint
{
    MyPoint() {}
    MyPoint(MyCoord x, MyCoord y) : x(x), y(y) {}
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
    MyEdge(MyIndex a, MyIndex b) : a(a), b(b) {}
    MyIndex a, b;
};

struct GfxStuffer
{
    GfxStuffer() { memset(this, 0, sizeof(GfxStuffer));  }

    GLenum type;
    Buf vbo, ibo_delabella, ibo_constraint;
    Buf vbo_voronoi, ibo_voronoi;

    GLuint texbuf;
    Buf tbo;

    #ifdef WITH_CDT
    Buf ibo_cdt;
    #endif
    MyCoord box[4];

    #ifdef CULLING
    MyCoord* max_tri_len;
    MyCoord* max_vor_len;
    MyCoord* max_con_len;


    MyIndex ConsByScale(MyIndex num, double scale)
    {
        const double thr = 6;
        if (num < 2)
            return num;
        MyCoord* len = max_con_len;
        if (len[num - 1] * scale > thr)
            return num;
        if (len[0] * scale <= thr)
            return 0;

        MyIndex lo = 1, hi = num - 2;

        while (lo < hi)
        {
            MyIndex med = (lo + hi) / 2;
            if (len[med] * scale > thr)
                lo = med + 1;
            else
                hi = med - 1;
        }

        return lo + 1;
    }

    MyIndex VoroByScale(MyIndex num, double scale)
    {
        const double thr = 6;
        if (num < 2)
            return num;
        MyCoord* len = max_vor_len;
        if (len[num - 1] * scale > thr)
            return num;
        if (len[0] * scale <= thr)
            return 0;

        MyIndex lo = 1, hi = num - 2;

        while (lo < hi)
        {
            MyIndex med = (lo + hi) / 2;
            if (len[med] * scale > thr)
                lo = med + 1;
            else
                hi = med - 1;
        }

        return lo + 1;
    }

    MyIndex TrisByScale(MyIndex num, double scale)
    {
        const double thr = 6;
        if (num < 2)
            return num;
        MyCoord* len = max_tri_len;
        if (len[num-1] * scale > thr)
            return num;
        if (len[0] * scale <= thr)
            return 0;

        MyIndex lo = 1, hi = num - 2;

        while (lo < hi)
        {
            MyIndex med = (lo + hi) / 2;
            if (len[med] * scale > thr)
                lo = med + 1;
            else
                hi = med - 1;
        }

        return lo + 1;
    }
    #endif

    struct Vao
    {
        Vao() : vao(0) {}

        GLuint Gen()
        {
            glGenVertexArrays(1, &vao);
            return vao;
        }

        void Bind()
        {
            glBindVertexArray(vao);
        }

        void Del()
        {
            if (vao)
                glDeleteVertexArrays(1, &vao);
            vao = 0;
        }

        GLuint vao;
    };

    Vao vao_main, vao_constraint, vao_voronoi, vao_cdt;

    GLuint prg; 
    GLint tfm;
    GLint low;
    GLint clr;
    GLint tex;

    void LoadProj(int vpw, int vph, double cx, double cy, double scale, double lx, double ly)
    {
        glViewport(0,0,vpw,vph);
        if (type == GL_DOUBLE)
        {
            glUseProgram(prg);

            GLdouble mat[4];

            mat[0] = scale / vpw;
            mat[1] = scale / vph;
            mat[2] = -cx;
            mat[3] = -cy;

            glUniform4dv(tfm, 1, mat);

            GLdouble lxy[4] = { -lx, -ly, 0.0, 0.0 }; // zw-spare
            glUniform4dv(low, 1, lxy);
        }
        else
        {
            glUseProgram(prg);

            GLfloat mat[4];

            mat[0] = (GLfloat)(scale / vpw);
            mat[1] = (GLfloat)(scale / vph);
            mat[2] = (GLfloat)(-cx);
            mat[3] = (GLfloat)(-cy);

            glUniform4fv(tfm, 1, mat);

            GLfloat lxy[4] = { -(GLfloat)lx, -(GLfloat)ly, 0.0f, 0.0f }; // zw-spare
            glUniform4fv(low, 1, lxy);
        }
    }

    void SetColor(float r, float g, float b, float a)
    {
        glUniform4f(clr, r,g,b,a);
    }

    void Destroy()
    {
        glBindTexture(GL_TEXTURE_BUFFER, 0);
        glDeleteTextures(1,&texbuf);
        tbo.Del();

        #ifdef CULLING
        if (max_tri_len)
            free(max_tri_len);

        #ifdef VORONOI
        if (max_vor_len)
            free(max_vor_len);
        #endif
        #endif

        vao_main.Del();
        vao_constraint.Del();
        vao_voronoi.Del();
        vao_cdt.Del();

        glUseProgram(0);
        if (prg)
            glDeleteProgram(prg);

        vbo.Del();
        ibo_delabella.Del();
        ibo_constraint.Del();

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
        MyIndex points,
        const MyPoint* cloud,
        MyIndex constrain_edges,
        const MyEdge* force
        #ifdef VORONOI
        , MyIndex voronoi_vertices,
        const MyPoint* voronoi_vtx_buf,
        MyIndex voronoi_indices,
        const MyIndex* voronoi_idx_buf
        #endif
        #ifdef WITH_CDT
        , const CDT::Triangulation<MyCoord>& cdt,
        const CDT::DuplicatesInfo& dups
        #endif
    )
    {
        Destroy();

        assert(points >= 0);
        assert(constrain_edges >= 0);
        #ifdef VORONOI
        assert(voronoi_vertices >= 0);
        assert(voronoi_indices >= 0);
        #endif

        vao_main.Gen();

        type = gl_e;

        static const char* vs_src[1] = { 0 };
        static const char* fs_src[1] = { 0 };

        #define CODE(...) #__VA_ARGS__

        fs_src[0] = CODE(#version 410\n
            uniform vec4 clr;
            uniform usamplerBuffer tex;
            layout (location = 0) out vec4 c;
            void main()
            {
                uint flags = texelFetch(tex, gl_PrimitiveID).r;
                if ((flags & 0x40) != 0)
                    c = clr + vec4(0.25);
                else
                    c = clr;
            }
        );

        if (type == GL_DOUBLE)
        {
            vs_src[0] = CODE(#version 410\n
                uniform dvec4 tfm;
                uniform dvec4 low;
                layout (location = 0) in dvec3 v;
                void main()
                {
                    gl_Position = vec4(dvec4(tfm.xy*(v.xy + tfm.zw * v.z + low.xy * v.z), 0.0lf, v.z));
                }
            );
        }
        else
        {
            vs_src[0] = CODE(#version 330\n
                uniform vec4 tfm;
                uniform vec4 low;
                layout (location = 0) in vec3 v;
                void main()
                {
                    gl_Position = vec4(tfm.xy*(v.xy + tfm.zw * v.z + low.xy * v.z), 0.0, v.z);
                }
            );
        }

        #undef CODE

        /*
        char nfolog[1025];
        int nfolen;
        */

        prg = glCreateProgram();
            
        GLuint vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, vs_src, 0);
        glCompileShader(vs);
        glAttachShader(prg,vs);

        /*
        glGetShaderInfoLog(vs,1024,&nfolen,nfolog);
        nfolog[nfolen]=0;
        printf("VS:\n%s\n\n",nfolog);
        */

        GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, fs_src, 0);
        glCompileShader(fs);
        glAttachShader(prg,fs);

        /*
        glGetShaderInfoLog(vs,1024,&nfolen,nfolog);
        nfolog[nfolen]=0;
        printf("FS:\n%s\n\n",nfolog);
        */

        glLinkProgram(prg);

        glDeleteShader(vs);
        glDeleteShader(fs);

        tfm = glGetUniformLocation(prg, "tfm");
        low = glGetUniformLocation(prg, "low");
        clr = glGetUniformLocation(prg, "clr");
        tex = glGetUniformLocation(prg, "tex");
        
        size_t gl_s = type == GL_DOUBLE ? sizeof(GLdouble) : sizeof(GLfloat);
        MyIndex tris_delabella = idb->GetNumOutputIndices() / 3;
        MyIndex contour = idb->GetNumBoundaryVerts();

        if (constrain_edges)
        {
            vao_constraint.Gen();
            ibo_constraint.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[2]) * (size_t)constrain_edges);
            GLuint* map = (GLuint*)ibo_constraint.Map();

            #ifdef CULLING
            struct ConSort
            {
                MyIndex e;
                MyCoord weight;
                bool operator < (const ConSort& b) const
                {
                    return weight > b.weight;
                }
            };

            ConSort* consort = (ConSort*)malloc(sizeof(ConSort) * (size_t)constrain_edges);
            assert(consort);

            for (MyIndex i = 0; i < constrain_edges; i++)
            {
                consort[i].e = i;
                MyIndex i0 = force[i].a;
                MyIndex i1 = force[i].b;
                MyCoord v01[2] = { cloud[i1].x - cloud[i0].x, cloud[i1].y - cloud[i0].y };
                MyCoord sqr = v01[0] * v01[0] + v01[1] * v01[1];
                consort[i].weight = sqrt(sqr);
            }

            std::sort(consort, consort + constrain_edges);

            max_con_len = (MyCoord*)malloc(sizeof(MyCoord) * (size_t)constrain_edges);
            assert(max_con_len);

            for (MyIndex i = 0; i < constrain_edges; i++)
            {
                MyIndex e = consort[i].e;
                map[2 * i + 0] = (GLuint)force[e].a;
                map[2 * i + 1] = (GLuint)force[e].b;
               
                max_con_len[i] = consort[i].weight;
            }

            free(consort);

            #else

            for (int i = 0; i < constrain_edges; i++)
            {
                map[2 * i + 0] = (GLuint)force[i].a;
                map[2 * i + 1] = (GLuint)force[i].b;
            }
            #endif
            ibo_constraint.Unmap();
        }

        vbo.Gen(GL_ARRAY_BUFFER, gl_s * 3 * (size_t)points);
        void* vbo_ptr = vbo.Map();

        // let's give a hand to gpu by centering vertices around 0,0
        box[0] = box[2] = cloud[0].x;
        box[1] = box[3] = cloud[0].y;
        for (MyIndex i = 0; i<points; i++)
        {
            box[0] = std::min(box[0], cloud[i].x);
            box[1] = std::min(box[1], cloud[i].y);
            box[2] = std::max(box[2], cloud[i].x);
            box[3] = std::max(box[3], cloud[i].y);
        }

        MyCoord vbo_x = (box[0]+box[2]) / 2;
        MyCoord vbo_y = (box[1]+box[3]) / 2;

        #if 0 // vbo centering? it's lossy!
        box[0] -= vbo_x;
        box[1] -= vbo_y;
        box[2] -= vbo_x;
        box[3] -= vbo_y;
        #else
        vbo_x = 0;
        vbo_y = 0;
        #endif

        if (type == GL_DOUBLE)
        {
            GLdouble* p = (GLdouble*)vbo_ptr;
            for (MyIndex i = 0; i<points; i++)
            {
                p[3*i+0] = (GLdouble)(cloud[i].x - vbo_x);
                p[3*i+1] = (GLdouble)(cloud[i].y - vbo_y);
                p[3*i+2] = (GLdouble)(1.0); // i%5; // color
            }
        }
        else
        {
            GLfloat* p = (GLfloat*)vbo_ptr;
            for (MyIndex i = 0; i<points; i++)
            {
                p[3*i+0] = (GLfloat)(cloud[i].x - vbo_x);
                p[3*i+1] = (GLfloat)(cloud[i].y - vbo_y);
                p[3*i+2] = (GLfloat)(1.0); // i%5; // color
            }
        }
        vbo.Unmap();

        tbo.Gen(GL_TEXTURE_BUFFER, tris_delabella);
        uint8_t* tbo_ptr = (uint8_t*)tbo.Map();

        GLuint* ibo_ptr = 0;
        #ifdef CULLING
        {
            struct TriSort
            {
                const DelaBella_Triangle* tri;
                MyCoord weight;
                bool operator < (const TriSort& b) const
                {
                    return weight > b.weight;
                }
            };
            TriSort* trisort = (TriSort*)malloc(sizeof(TriSort) * (size_t)tris_delabella);
            assert(trisort);
            const DelaBella_Triangle* dela = idb->GetFirstDelaunaySimplex();
            for (MyIndex i = 0; i < tris_delabella; i++)
            {
                trisort[i].tri = dela;
                MyCoord v01[2] = { dela->v[1]->x - dela->v[0]->x, dela->v[1]->y - dela->v[0]->y };
                MyCoord v12[2] = { dela->v[2]->x - dela->v[1]->x, dela->v[2]->y - dela->v[1]->y };
                MyCoord v20[2] = { dela->v[0]->x - dela->v[2]->x, dela->v[0]->y - dela->v[2]->y };

                MyCoord sqr = v01[0] * v01[0] + v01[1] * v01[1];
                sqr = std::max(sqr, v12[0] * v12[0] + v12[1] * v12[1]);
                sqr = std::max(sqr, v20[0] * v20[0] + v20[1] * v20[1]);

                trisort[i].weight = sqrt(sqr);
                dela = dela->next;
            }

            std::sort(trisort, trisort + tris_delabella);

            max_tri_len = (MyCoord*)malloc(sizeof(MyCoord) * (size_t)tris_delabella);
            assert(max_tri_len);

            ibo_delabella.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3])* (size_t)tris_delabella + sizeof(GLuint) * (size_t)contour);
            ibo_ptr = (GLuint*)ibo_delabella.Map();
            for (MyIndex i = 0; i < tris_delabella; i++)
            {
                dela = trisort[i].tri;
                MyIndex v0 = dela->v[0]->i;
                MyIndex v1 = dela->v[1]->i;
                MyIndex v2 = dela->v[2]->i;

                ibo_ptr[3 * i + 0] = (GLuint)v0;
                ibo_ptr[3 * i + 1] = (GLuint)v1;
                ibo_ptr[3 * i + 2] = (GLuint)v2;

                tbo_ptr[i] = dela->flags;

                max_tri_len[i] = trisort[i].weight;
            }

            free(trisort);
        }
        #else
        {
            ibo_delabella.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delabella + sizeof(GLuint) * contour);
            ibo_ptr = (GLuint*)ibo_delabella.Map();
            const DelaBella_Triangle* dela = idb->GetFirstDelaunaySimplex();
            for (int i = 0; i < tris_delabella; i++)
            {
                MyIndex v0 = dela->v[0]->i;
                MyIndex v1 = dela->v[1]->i;
                MyIndex v2 = dela->v[2]->i;

                ibo_ptr[3 * i + 0] = (GLuint)v0;
                ibo_ptr[3 * i + 1] = (GLuint)v1;
                ibo_ptr[3 * i + 2] = (GLuint)v2;

                tbo_ptr[i] = dela->flags;

                dela = dela->next;
            }
        }
        #endif

        tbo.Unmap();
        tbo_ptr = 0;

        glGenTextures(1, &texbuf);
        glBindTexture(GL_TEXTURE_BUFFER, texbuf);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_R8UI, tbo.buf);

        typedef GLuint tri_in_ibo[3];

        const DelaBella_Vertex* vert = idb->GetFirstBoundaryVertex();
        for (MyIndex i = 0; i<contour; i++)
        {
            ibo_ptr[i + 3*tris_delabella] = (GLuint)vert->i;
            vert = vert->next;
        }

        ibo_delabella.Unmap();

        #ifdef VORONOI
        vao_voronoi.Gen();

        vbo_voronoi.Gen(GL_ARRAY_BUFFER, gl_s * 3 * (size_t)voronoi_vertices);
        ibo_voronoi.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * (size_t)voronoi_indices);
        void* vbo_voronoi_ptr = vbo_voronoi.Map();
        GLuint* ibo_voronoi_ptr = (GLuint*)ibo_voronoi.Map();

        if (type == GL_DOUBLE)
        {
            GLdouble* p = (GLdouble*)vbo_voronoi_ptr;
            for (MyIndex i = 0; i < voronoi_vertices; i++)
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
            for (MyIndex i = 0; i < voronoi_vertices; i++)
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

        #ifdef CULLING
        {
            #ifdef VORONOI_POLYS
            // todo
            for (int i = 0; i < voronoi_indices; i++)
                ibo_voronoi_ptr[i] = (GLuint)(voronoi_idx_buf[i]);
            max_vor_len = 0;
            #else
            MyIndex edges = voronoi_indices / 2;
            struct VorSort
            {
                MyIndex e;
                MyCoord weight;
                bool operator < (const VorSort& b) const
                {
                    return weight > b.weight;
                }
            };
            VorSort* vorsort = (VorSort*)malloc(sizeof(VorSort) * (size_t)edges);
            assert(vorsort);
            for (MyIndex i = 0; i < edges; i++)
            {
                vorsort[i].e = i;
                MyIndex i0 = voronoi_idx_buf[2 * i];
                MyIndex i1 = voronoi_idx_buf[2 * i + 1];

                if (i0 >= tris_delabella || i1 >= tris_delabella)
                {
                    // inf edge
                    vorsort[i].weight = INFINITY; // -1
                }
                else
                {
                    MyCoord v01[2] = { voronoi_vtx_buf[i1].x - voronoi_vtx_buf[i0].x, voronoi_vtx_buf[i1].y - voronoi_vtx_buf[i0].y };
                    MyCoord sqr = v01[0] * v01[0] + v01[1] * v01[1];
                    vorsort[i].weight = sqrt(sqr);
                }
            }

            std::sort(vorsort, vorsort + edges);

            max_vor_len = (MyCoord*)malloc(sizeof(MyCoord)* (size_t)edges);
            assert(max_vor_len);

            for (MyIndex i = 0; i < edges; i++)
            {
                MyIndex e = vorsort[i].e;
                ibo_voronoi_ptr[2*i] = (GLuint)(voronoi_idx_buf[2*e+0]);
                ibo_voronoi_ptr[2*i+1] = (GLuint)(voronoi_idx_buf[2*e+1]);

                max_vor_len[i] = vorsort[i].weight;
            }

            free(vorsort);

            #endif
        }
        #else
        {
            for (int i = 0; i < voronoi_indices; i++)
                ibo_voronoi_ptr[i] = (GLuint)(voronoi_idx_buf[i]);
        }
        #endif

        vbo_voronoi.Unmap();
        ibo_voronoi.Unmap();
        #endif

        #ifdef WITH_CDT
        vao_cdt.Gen();

        MyIndex tris_cdt = (MyIndex)cdt.triangles.size();
        ibo_cdt.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_cdt);
        {
            ibo_ptr = (GLuint*)ibo_cdt.Map();

            if (dups.duplicates.size())
            {
                // let's make inverse mapping first
                MyIndex* invmap = 0;
                MyIndex invmap_size = (MyIndex)dups.mapping.size() - (MyIndex)dups.duplicates.size();
                invmap = (MyIndex*)malloc(sizeof(MyIndex) * invmap_size);
                assert(invmap);

                for (MyIndex i = 0; i < (MyIndex)dups.mapping.size(); i++)
                    invmap[dups.mapping[i]] = i;

                for (MyIndex i = 0; i < tris_cdt; i++)
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
                for (MyIndex i = 0; i < tris_cdt; i++)
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

        vao_main.Bind();
        vbo.Bind();
        ibo_delabella.Bind();
        if (type == GL_DOUBLE)
            glVertexAttribLPointer(0, 3, GL_DOUBLE, 0, 0);
        else
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        vao_constraint.Bind();
        vbo.Bind();
        ibo_constraint.Bind();
        if (type == GL_DOUBLE)
            glVertexAttribLPointer(0, 3, GL_DOUBLE, 0, 0);
        else
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        #ifdef VORONOI
        vao_voronoi.Bind();
        vbo_voronoi.Bind();
        ibo_voronoi.Bind();
        if (type == GL_DOUBLE)
            glVertexAttribLPointer(0, 3, GL_DOUBLE, 0, 0);
        else
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);
        #endif

        #ifdef WITH_CDT
        vao_cdt.Bind();
        vbo.Bind();
        ibo_cdt.Bind();
        if (type == GL_DOUBLE)
            glVertexAttribLPointer(0, 3, GL_DOUBLE, 0, 0);
        else
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);
        #endif

        glBindVertexArray(0);
    }
};

#ifdef BENCH
extern uint64_t sorting_bench;
struct Bench
{
    void operator += (const Bench& b)
    {
        removing_dups += b.removing_dups;
        triangulation += b.triangulation;
        constrain_edges += b.constrain_edges;
        erase_super += b.erase_super;
        flood_fill += b.flood_fill;
        polygons += b.polygons;
    }

    void operator /= (int b)
    {
        removing_dups /= b;
        triangulation /= b;
        constrain_edges /= b;
        erase_super /= b;
        flood_fill /= b;
        polygons /= b;
    }

    uint64_t removing_dups;
    uint64_t triangulation;
    uint64_t constrain_edges;
    uint64_t erase_super;
    uint64_t flood_fill;
    uint64_t polygons;
};
#endif

#ifdef BENCH
int bench_main(int argc, char* argv[])
{
    Bench* idb_bench = (Bench*)argv[2];
    Bench* cdt_bench = (Bench*)argv[3];
    Bench* fad_bench = (Bench*)argv[4];
    Bench* fmt_bench = (Bench*)argv[5];
    Bench* del_bench = (Bench*)argv[6];
    char* dist = argv[7];
    char* bias = argv[8];

#else
int main(int argc, char* argv[])
{
    const char* dist = "cir";
    const char* bias = "";
#endif

	#ifdef _WIN32
	SetProcessDPIAware();
    struct handler
    {
        static BOOL WINAPI routine(DWORD CtrlType)
        {
            printf("\nterminating!\n");
            exit(-1);
            return TRUE;
        }
    };
    SetConsoleCtrlHandler(handler::routine, TRUE);
	#else
    // ...
    #endif

	if (argc<2)
	{
        printf("usage: %s <input.xy> [output.abc]\n", argc > 0 ? argv[0] : "<executable>");
		printf("required argument (.xy file name / number of random points) missing, terminating!\n");
		return -1;
	}

	std::vector<MyPoint> cloud;
    std::vector<MyEdge> force;

    #ifdef BENCH
    FILE* f = 0;
    #else
	FILE* f = fopen(argv[1],"r");
    #endif

    if (!f)
    {
        MyIndex n = atoi(argv[1]);

        if (n <= 2)
        {
            printf("can't open %s file, terminating!\n", argv[1]);
            return -1;
        }
        printf("generating random " IDXF " points\n", n);
        std::random_device rd{};
        uint64_t seed = rd();
        std::mt19937_64 gen{ seed };
        printf("SEED = 0x%016llX\n", (long long unsigned int)seed);

        std::uniform_real_distribution<MyCoord> d_uni((MyCoord)-1.0, (MyCoord)+1.0);
        std::normal_distribution<MyCoord> d_std{(MyCoord)0.0,(MyCoord)2.0};
        std::gamma_distribution<MyCoord> d_gam((MyCoord)0.1,(MyCoord)2.0);

        MyCoord max_coord = sizeof(MyCoord) < 8 ? /*float*/0x1.p31 : /*double*/0x1.p255;

        MyCoord bias_xy[] = { 0,0 };
        if (bias[0])
        {
            bias_xy[0] = 50;
            bias_xy[1] = 50;
        }

        if (strcmp(dist, "uni") == 0)
        {
            for (MyIndex i = 0; i < n; i++)
            {
                MyPoint p = { d_uni(gen) + bias_xy[0], d_uni(gen) + bias_xy[1] };
                // please, leave some headroom for arithmetics!
                assert(std::abs(p.x) <= max_coord && std::abs(p.y) <= max_coord);
                cloud.push_back(p);
            }
        }
        else
        if (strcmp(dist, "std") == 0)
        {
            for (MyIndex i = 0; i < n; i++)
            {
                MyPoint p = { d_std(gen) + bias_xy[0], d_std(gen) + bias_xy[1] };
                // please, leave some headroom for arithmetics!
                assert(std::abs(p.x) <= max_coord && std::abs(p.y) <= max_coord);
                cloud.push_back(p);
            }
        }
        else
        if (strcmp(dist, "gam") == 0)
        {
            for (MyIndex i = 0; i < n; i++)
            {
                MyPoint p = { d_gam(gen) + bias_xy[0], d_gam(gen) + bias_xy[1] };

                // please, leave some headroom for arithmetics!
                assert(std::abs(p.x) <= max_coord && std::abs(p.y) <= max_coord);
                cloud.push_back(p);
            }
        }
        else
        if (strcmp(dist, "sym") == 0)
        {
            n /= 8;
            for (MyIndex i = 0; i < n; i++)
            {
                MyPoint p = { d_gam(gen), d_gam(gen) };

                // please, leave some headroom for arithmetics!
                assert(std::abs(p.x) <= max_coord && std::abs(p.y) <= max_coord);

                cloud.push_back(MyPoint(p.x+bias_xy[0], p.y+bias_xy[1]));
                p.x = -p.x;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));
                p.y = -p.y;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));
                p.x = -p.x;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));

                MyCoord s = p.x;
                p.x = -p.y;
                p.y = s;

                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));
                p.x = -p.x;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));
                p.y = -p.y;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));
                p.x = -p.x;
                cloud.push_back(MyPoint(p.x + bias_xy[0], p.y + bias_xy[1]));

            }
            n *= 8;
        }
        else
        if (strcmp(dist, "cir") == 0)
        {
            for (MyIndex i = 0; i < n; i++)
            {
                MyCoord x = d_std(gen);
                MyCoord y = d_std(gen);

                MyCoord l = sqrt(x * x + y * y);

                x = x / l;
                y = y / l;

                MyPoint p{ x + bias_xy[0],y + bias_xy[1] };
                cloud.push_back(p);
            }
        }
        else
        // HEXGRID
        {
            const double x = 0x5af2efc1.p-30;
            const double y = 0x348268e0.p-30;
            const double r = 0x6904d1c1.p-30;
            const MyPoint p[4] =
            {
                {0,0},
                {0,(MyCoord)(2 * y)},
                {(MyCoord)x,(MyCoord)(r + y)},
                {(MyCoord)x,(MyCoord)(r + 3 * y)}
            };

            const MyCoord dx = (MyCoord)(2 * x);
            const MyCoord dy = (MyCoord)(2 * (r + y));

            int rows = (int)ceil(sqrt(n / 7.0));
            int cols = (n + 2 * rows) / (4 * rows);

            n = rows * cols * 4;

            printf("%d\n", rows * cols * 4);

            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < cols; col++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        MyPoint q = p[i];
                        q.x += col * dx + bias_xy[0];
                        q.y += row * dy + bias_xy[1];
                        cloud.push_back(q);
                    }
                }
            }
        }
             
        if (1) // generate 1/10 constrain edges
        {
            MyIndex m = n / 10;

            // init sub[] with all n point indices
            MyIndex* sub = (MyIndex*)malloc(sizeof(MyIndex) * (size_t)n);
            for (MyIndex i = 0; i < n; i++)
                sub[i] = i;

            // pick m random ones from n
            // place them as first m items of sub[]
            for (MyIndex i = 0; i < m; i++)
            {
                MyIndex j = i + gen() % ((size_t)n - i);
                MyIndex r = sub[j];
                sub[j] = sub[i];
                sub[i] = r;
            }

            std::vector<MyPoint> xxx;
            for (MyIndex i = 0; i < m; i++)
                xxx.push_back(cloud[(size_t)sub[i]]);

            IDelaBella* helper = IDelaBella::Create();
            helper->Triangulate(m, &xxx.data()->x, &xxx.data()->y, sizeof(MyPoint));

            // to avoid forth-and-back edges repeatitions,
            // traverse all faces but use edges with 
            // ascending y or in case of flat y use only if ascending x

            const DelaBella_Triangle* dela = helper->GetFirstDelaunaySimplex();
            while (dela)
            {
                for (int a = 0, b = 1, c = 2; a < 3; b = c, c = a, a++)
                {
                    //if ((dela->f[c]->flags & 0x80) || dela->v[a]->i < dela->v[b]->i)
                    if (!(dela->f[c]->flags & 0x80) && dela->v[a]->i < dela->v[b]->i)
                        force.push_back(MyEdge(sub[dela->v[a]->i], sub[dela->v[b]->i]));
                }

                dela = dela->next;
            }

            helper->Destroy();
            free(sub);
        }
        
	}
    else
    {
        int r=0,n=0,c=0;
        r = fscanf(f, "%d %d", &n, &c);

        if (n == -1 && c == -1)
        {
            // auto constraining flavour
            int start = 0;
            int current = 0;
            while (1)
            {
                double dbl_x, dbl_y;
                // allow variety of separators and extra fields till end of the line
                int n = fscanf(f, "%lf%*[,; \v\t]%lf%*[^\n]", &dbl_x, &dbl_y);

                if (n <= 0)
                {
                    char check[2];
                    char* res = fgets(check,2,f);

                    if (start == current)
                        break; // eof?

                    // end of contour
                    // generate edges loop
                    for (int i = start, p = current-1; i < current; p=i, i++)
                    {
                        MyEdge e = { p, i };
                        force.push_back(e);
                    }

                    start = current;
                }
                else
                if (n==2)
                {
                    MyCoord x = (MyCoord)dbl_x;
                    MyCoord y = (MyCoord)dbl_y;

                    MyPoint p = { x,y };
                    cloud.push_back(p);
                    current++;
                }
                else
                    break;
            }
        }
        else
        {
            for (int i = 0; n < 0 || i < n; i++)
            {
                double dbl_x, dbl_y;
                // allow variety of separators and extra fields till end of the line
                int n = fscanf(f, "%lf%*[,; \v\t]%lf%*[^\n]", &dbl_x, &dbl_y);

                MyCoord x = (MyCoord)dbl_x;
                MyCoord y = (MyCoord)dbl_y;

                MyPoint p = { x,y };
                cloud.push_back(p);
            }

            for (int i = 0; i < c; i++)
            {
                MyIndex a, b;
                r = fscanf(f, "" IDXF " " IDXF "", &a, &b);
                MyEdge e = { a, b };
                force.push_back(e);
            }
        }

        fclose(f);
    }

    MyIndex points = (MyIndex)cloud.size();

    #ifdef WITH_CDT

        std::vector<CDT::V2d<MyCoord>> nodups;
        for (size_t i = 0; i < cloud.size(); i++)
        {
            CDT::V2d<MyCoord> v;
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
        
        #ifdef BENCH
        cdt_bench->removing_dups = t1 - t0;
        #endif

        printf("cdt triangulation... ");
        CDT::Triangulation<MyCoord> cdt;
        cdt.insertVertices(nodups);

        uint64_t t2 = uSec();
        printf("%d ms\n", (int)((t2 - t1) / 1000));
        
        #ifdef BENCH
        cdt_bench->triangulation = t2 - t1;
        #endif

        if (force.size()>0)
        {
            printf("cdt forcing edges... ");
            cdt.insertEdges(edges);
            uint64_t t3 = uSec();
            printf("%d ms\n", (int)((t3 - t2) / 1000));
            
            #ifdef BENCH
            cdt_bench->constrain_edges = t3 - t2;
            #endif

            printf("cdt copying... ");
            uint64_t t4 = uSec();
            CDT::Triangulation<MyCoord> cdt2 = cdt;
            uint64_t t5 = uSec();
            printf("%d ms\n", (int)((t5 - t4) / 1000));

            printf("cdt erasing outer and holes... ");
            uint64_t t6 = uSec();
            cdt2.eraseOuterTrianglesAndHoles();
            uint64_t t7 = uSec();

            printf("%d ms\n", (int)((t7 - t6) / 1000));
            
            #ifdef BENCH
            cdt_bench->flood_fill = t7 - t6;
            #endif

            printf("CDT has %d faces after eraseOuterTrianglesAndHoles()\n", (int)cdt2.triangles.size());

        }

        uint64_t t4 = uSec();
        printf("cdt erasing super... ");
        cdt.eraseSuperTriangle();
        uint64_t t5 = uSec();
        printf("%d ms\n", (int)((t5 - t4) / 1000));
        
        #ifdef BENCH
        cdt_bench->erase_super = t5 - t4;
        #endif

        printf("CDT triangles = %d\n", (int)cdt.triangles.size());

        uint64_t t_p0 = uSec();
        std::vector<CDT::Poly> cdt_polys = CDT::Polygonize(cdt);
        uint64_t t_p1 = uSec();

        printf("CDT POLYS = " IDXF " (in %d ms)\n", (MyIndex)cdt_polys.size(), (int)((t_p1 - t_p0) / 1000));
        
        #ifdef BENCH
        cdt_bench->polygons = t_p1 - t_p0;
        #endif

    #endif

    #ifdef WITH_FADE
    if (points < 500000 || strcmp(dist,"gam") || strcmp(bias,"+"))
    {
        printf("running fade ...");
        uint64_t t0 = uSec(), t1, t2;
        Point2** handles = (Point2**)malloc(sizeof(Point2*) * cloud.size());
        MyIndex tris_fade;
        try {
            Fade_2D dt;
            dt.insert((int)cloud.size(), &cloud.data()->x, handles);
            tris_fade = (MyIndex)dt.numberOfTriangles();

            t1 = uSec();

            // we need to map using handles
            std::vector<Segment2> ve;
            ve.reserve(force.size());
            for (int i = 0; i < force.size(); i++)
            {
                Segment2 e{ *handles[force[i].a], *handles[force[i].b] };
                ve.push_back(e);
            }
            auto cgr = dt.createConstraint(ve, ConstraintInsertionStrategy::CIS_CONSTRAINED_DELAUNAY, true);

            t2 = uSec();

            auto zon = dt.createZone(cgr, ZL_INSIDE, false);

            //#ifndef BENCH
            //zon->show("zone.ps", false, false);
            //#endif

            printf("ZONE has %d triangles", (int)zon->getNumberOfTriangles());

            //delete cgr;
            //delete zon;
        }
        catch (...)
        {
            tris_fade = 0;
        }
        free(handles);
        uint64_t t3 = uSec();
        printf("(%d ms)\nFADE %d tris\n", (int)((t1 - t0) / 1000), (int)tris_fade);

        #ifdef BENCH
        fad_bench->triangulation = t1 - t0;
        if (tris_fade != (MyIndex)cdt.triangles.size())
            fad_bench->removing_dups = 1000000000;
        else
            fad_bench->removing_dups = 0;
        fad_bench->constrain_edges = t2-t1;
        fad_bench->flood_fill = t3 - t2;
        #endif
    }
    {
        printf("running fade_mt ...");
        uint64_t t0 = uSec(), t1;
        Point2** handles = (Point2**)malloc(sizeof(Point2*) * cloud.size());
        MyIndex tris_fade;
        try {
            Fade_2D dt;
            dt.setNumCPU(0);
            dt.insert((int)cloud.size(), &cloud.data()->x, handles);
            tris_fade = (MyIndex)dt.numberOfTriangles();

            t1 = uSec();

            // we need to map using handles
            if (points < 500000 || strcmp(dist, "gam"))
            {
                std::vector<Segment2> ve;
                ve.reserve(force.size());
                for (int i = 0; i < force.size(); i++)
                {
                    Segment2 e{ *handles[force[i].a], *handles[force[i].b] };
                    ve.push_back(e);
                }
                auto cgr = dt.createConstraint(ve, ConstraintInsertionStrategy::CIS_CONSTRAINED_DELAUNAY, false);

                t2 = uSec();

                auto zon = dt.createZone(0, ZL_GLOBAL, false);
            }
            else
                t2 = t1;

            //delete cgr;
            //delete zon;
        }
        catch (...)
        {
            tris_fade = 0;
        }
        free(handles);
        uint64_t t3 = uSec();
        printf("(%d ms)\nFADE_MT %d tris\n", (int)((t1 - t0) / 1000), (int)tris_fade);

        #ifdef BENCH
        fmt_bench->triangulation = t1 - t0;
        if (tris_fade != (MyIndex)cdt.triangles.size())
            fmt_bench->removing_dups = 1000000000;
        else
            fmt_bench->removing_dups = 0;
        fmt_bench->constrain_edges = t2 - t1;
        fmt_bench->flood_fill = t3 - t2;
        #endif
    }

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
        MyIndex tris_delaunator = d ? (MyIndex)d->triangles.size() / 3 : 0;
        uint64_t t1 = uSec();
        printf("elapsed %d ms\n", (int)((t1-t0)/1000));

        #ifdef BENCH
        del_bench->triangulation = t1 - t0;
        if (tris_delaunator != (MyIndex)cdt.triangles.size())
            del_bench->removing_dups = 1000000000;
        else
            del_bench->removing_dups = 0;
        #endif
        
        printf("delaunator triangles: " IDXF "\n", tris_delaunator);

        delete d;
    }
    #endif

	IDelaBella* idb = IDelaBella::Create();

    #ifndef BENCH
	idb->SetErrLog(errlog, stdout);
    #endif

    printf("running delabella...\n");
    uint64_t t6 = uSec();
    MyIndex verts = idb->Triangulate(points, &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
    //idb->CheckTopology();
    #ifdef BENCH
    idb_bench->removing_dups = sorting_bench;
    idb_bench->triangulation = uSec()-t6 - sorting_bench;
    #endif

    MyIndex tris_delabella = verts > 0 ? verts / 3 : 0;
    MyIndex contour = idb->GetNumBoundaryVerts();
    MyIndex non_contour = idb->GetNumInternalVerts();
    MyIndex vert_num = contour + non_contour;

    #ifdef VORONOI
    //printf("Polygonizing for VD\n");
    //idb->Polygonize(); // optional

    printf("Generating VD vertices\n");
    MyIndex voronoi_vertices = idb->GenVoronoiDiagramVerts(0, 0, 0);
    MyPoint* voronoi_vtx_buf = (MyPoint*)malloc((size_t)voronoi_vertices * sizeof(MyPoint));
    assert(voronoi_vtx_buf);
    idb->GenVoronoiDiagramVerts(&voronoi_vtx_buf->x, &voronoi_vtx_buf->y, sizeof(MyPoint));

    printf("Generating VD indices\n");
    #ifdef VORONOI_POLYS
    // testing... will remove
    MyIndex voronoi_closed_indices;
    MyIndex voronoi_indices = idb->GenVoronoiDiagramPolys(0, 0, 0);
    MyIndex* voronoi_idx_buf = (MyIndex*)malloc(voronoi_indices * sizeof(MyIndex));
    assert(voronoi_idx_buf);
    idb->GenVoronoiDiagramPolys(voronoi_idx_buf, sizeof(MyIndex), &voronoi_closed_indices);
    #else
    MyIndex voronoi_closed_indices = 0;
    MyIndex voronoi_indices = idb->GenVoronoiDiagramEdges(0, 0);
    MyIndex* voronoi_idx_buf = (MyIndex*)malloc((size_t)voronoi_indices * sizeof(MyIndex));
    assert(voronoi_idx_buf);
    idb->GenVoronoiDiagramEdges(voronoi_idx_buf, sizeof(MyIndex));
    #endif

    printf("VD vertices = " IDXF ", indices = " IDXF "\n", voronoi_vertices, voronoi_indices);
    #endif

    if (force.size() > 0)
    {
        #ifdef BENCH
        idb_bench->constrain_edges = uSec();
        #endif

        idb->ConstrainEdges((MyIndex)force.size(), &force.data()->a, &force.data()->b, (int)sizeof(MyEdge));
        //idb->CheckTopology();

        #ifdef BENCH
        idb_bench->constrain_edges = uSec() - idb_bench->constrain_edges;
        #endif

        uint64_t ff0 = uSec();
        MyIndex num_interior = idb->FloodFill(false, 0);
        //idb->CheckTopology();
        uint64_t ff1 = uSec();

        printf("interior %d faces in %d ms\n", num_interior, (int)((ff1 - ff0) / 1000));

        #ifdef BENCH
        idb_bench->flood_fill = ff1-ff0;
        #endif
    }

    #ifdef BENCH
    idb_bench->erase_super = 0;
    #endif

    const DelaBella_Triangle** dela_polys = (const DelaBella_Triangle**)malloc(sizeof(const DelaBella_Triangle*) * (size_t)tris_delabella);

    #ifdef BENCH
    idb_bench->polygons = uSec();
    #endif    
    
    MyIndex polys_delabella = idb->Polygonize(dela_polys);
    //idb->CheckTopology();
    
    #ifdef BENCH
    idb_bench->polygons = uSec() - idb_bench->polygons;
    #endif

    printf("delabella triangles: " IDXF "\n", tris_delabella);
    printf("delabella contour: " IDXF "\n", contour);
    printf("delabella polygons: " IDXF "\n", polys_delabella);

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

        MyIndex tris_cdt = (MyIndex)cdt.triangles.size();

        if (tris_delabella != tris_cdt || polys_delabella != cdt_polys.size())
            printf("WARNING! Results are not comparable - different number of tris or polys\n");
        else
        {
            // note:
            // currently, any difference between constrained edges passing 
            // internal edge of a polygon won't be detected !!!

            MyIndex poly_indices = 2 * polys_delabella + tris_delabella;

            struct MyPoly
            {
                MyIndex size;
                MyIndex offs;
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
                        for (MyIndex i = 0; i < p.size; i++)
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
            for (MyIndex p = 0, n = 0; p < (MyIndex)cdt_polys.size(); p++)
            {
                MyIndex s = 0;

                for (int i = 0; i < 3; i++)
                {
                    MyIndex j = n + s;
                    MyIndex k = cdt.triangles[cdt_polys[p][0]].vertices[i];
                    cdt_v[j].x = cdt.vertices[k].x;
                    cdt_v[j].y = cdt.vertices[k].y;
                    s++;
                }

                for (MyIndex i = 1; i < (MyIndex)cdt_polys[p].size(); i++)
                {
                    MyIndex j = n + s;
                    MyIndex k = cdt.triangles[cdt_polys[p][i]].vertices[0];
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
            for (MyIndex p = 0, n = 0; p < polys_delabella; p++)
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
            for (MyIndex p = 0; p < polys_delabella; p++)
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

    free(dela_polys);

    if (argc>=3)
    {
        f = fopen(argv[2],"w");
        if (f)
        {
            for (MyIndex i=0; i<tris_delabella; i++)
            {
                const DelaBella_Triangle* dela = idb->GetFirstDelaunaySimplex();
                fprintf(f,"" IDXF " " IDXF " " IDXF "\n",
                    dela->v[0]->i,
                    dela->v[1]->i,
                    dela->v[2]->i);
                dela = dela->next;
            }
            fclose(f);
        }
    }

    if (argc >= 4)
    {
        // write back input
        f = fopen(argv[3], "w");
        if (f)
        {
            fprintf(f, "%d %d\n", (int)cloud.size(), (int)force.size());
            for (MyIndex i = 0; i < cloud.size(); i++)
                fprintf(f, "%f %f\n", (double)cloud[i].x, (double)cloud[i].y);
            for (MyIndex i = 0; i < force.size(); i++)
                fprintf(f, "%d %d\n", (int)force[i].a, (int)force[i].b);
            fclose(f);
        }
    }

#ifdef BENCH
    #ifdef VORONOI
    free(voronoi_idx_buf);
    free(voronoi_vtx_buf);
    #endif
    idb->Destroy();
    return 0;
#endif

    SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1");
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

// we want it at least 3.3 but would be nice to have 4.1 or above
//    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
//    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	// create viewer wnd
    int width = 800, height = 800;
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

    int glsl_ver = 0;
    const char* glsl_str = (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);
    if (glsl_str && glGetError() == GL_NO_ERROR)
    {
        int v[2];
        int from, to;
        if (2 == sscanf(glsl_str, "%d.%n%d%n", v + 0, &from, v + 1, &to))
        {
            int num = to - from;
            if (num == 1)
                v[1] *= 10;
            else
            while (num > 2)
            {
                v[1] /= 10;
                num--;
            }
            glsl_ver = v[0] * 100 + v[1];
        }
    }

    if (glsl_ver < 330)
    {
        printf("GLSL %d is too weak - terminating!", glsl_ver);
        idb->Destroy();
        return -1;
    }

	if (!BindGL())
	{
		printf("Can't bind to necessary GL functions, terminating!\n");
		idb->Destroy();
		return -1;
	}

    GfxStuffer gfx;

    //glsl_ver = 330;
    printf("preparing graphics for GLSL %d...\n", glsl_ver);

    gfx.Upload( glsl_ver >= 410 ? GL_DOUBLE : GL_FLOAT,
                idb,
                (MyIndex)cloud.size(),
                cloud.data(),
                (MyIndex)force.size(),
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
    double lx = 0.0;
    double ly = 0.0;
    double scale = 2.0 * fmin((double)vpw/(gfx.box[2]-gfx.box[0]),(double)vph/(gfx.box[3]-gfx.box[1]));
    int zoom = -3+(int)round(log(scale) / log(1.01));

    int drag_x, drag_y, drag_zoom;
    double drag_cx, drag_cy;
    double drag_lx, drag_ly;
    int drag = 0;

    glPrimitiveRestartIndex(~(GLuint)0);
    glEnable(GL_PRIMITIVE_RESTART);

    glEnable(GL_BLEND);

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

    float lohi[2];

    glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lohi);
    float thin = 1.0f, thick = 3.0f;
    if (thick > lohi[1])
        thick = lohi[1];

    glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, lohi);
    float dot = 1.0f, blob = 3.0f;
    if (blob > lohi[1])
        blob = lohi[1];

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

                        if (gfx.type == GL_DOUBLE)
                        {
                            predicates::detail::Expansion<GLdouble, 1> adx, ady;
                            adx.push_back((GLdouble)(-2.0 * dx / scale));
                            ady.push_back((GLdouble)(2.0 * dy / scale));

                            auto edx =
                                predicates::detail::ExpansionBase<GLdouble>::Plus((GLdouble)drag_cx, (GLdouble)drag_lx) + adx;

                            auto edy =
                                predicates::detail::ExpansionBase<GLdouble>::Plus((GLdouble)drag_cy, (GLdouble)drag_ly) + ady;

                            cx = edx.m_size > 0 ? (double)edx[edx.m_size - 1] : 0.0;
                            lx = edx.m_size > 1 ? (double)edx[edx.m_size - 2] : 0.0;

                            cy = edy.m_size > 0 ? (double)edy[edy.m_size - 1] : 0.0;
                            ly = edy.m_size > 1 ? (double)edy[edy.m_size - 2] : 0.0;
                        }
                        else
                        {
                            predicates::detail::Expansion<GLfloat, 1> adx, ady;
                            adx.push_back((GLfloat)( - 2.0 * dx / scale));
                            ady.push_back((GLfloat)(2.0 * dy / scale));

                            auto edx =
                                predicates::detail::ExpansionBase<GLfloat>::Plus((GLfloat)drag_cx, (GLfloat)drag_lx) + adx;

                            auto edy =
                                predicates::detail::ExpansionBase<GLfloat>::Plus((GLfloat)drag_cy, (GLfloat)drag_ly) + ady;

                            cx = edx.m_size > 0 ? (double)edx[edx.m_size - 1] : 0.0;
                            lx = edx.m_size > 1 ? (double)edx[edx.m_size - 2] : 0.0;

                            cy = edy.m_size > 0 ? (double)edy[edy.m_size - 1] : 0.0;
                            ly = edy.m_size > 1 ? (double)edy[edy.m_size - 2] : 0.0;
                        }
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
                            drag_lx = lx;
                            drag_ly = ly;
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

        gfx.LoadProj(vpw,vph, cx,cy, scale, lx,ly);
        glUniform1i(gfx.tex, 1); // nothing is bound there

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);

        gfx.vao_main.Bind();

        #ifdef CULLING
            MyIndex cull = gfx.TrisByScale(tris_delabella, scale);
            MyIndex cons_cull = gfx.ConsByScale((MyIndex)force.size(), scale);
            #ifdef VORONOI
            MyIndex voro_cull = 2 * gfx.VoroByScale(voronoi_indices / 2, scale);
            #endif
        #else
            MyIndex cull = tris_delabella;
            MyIndex cons_cull = (MyIndex)force.size();
            #ifdef VORONOI
            MyIndex voro_cull = voronoi_indices;
            #endif
        #endif

        // grey fill
        // TODO: switch to trifan using contour indices
        if (show_f)
        {
            glUniform1i(gfx.tex, 0);
            gfx.SetColor(0.2f, 0.2f, 0.2f, 1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDrawElements(GL_TRIANGLES, (GLsizei)cull/*tris_delabella*/ * 3, GL_UNSIGNED_INT, 0);
            glUniform1i(gfx.tex, 1); // nothing is bound there
        }

        // paint constraints
        MyIndex constrain_indices = 2*cons_cull;
        if (constrain_indices && show_c)
        {
            gfx.vao_constraint.Bind();

            glLineWidth(thick);
            gfx.SetColor(.9f, .9f, .9f, 1.0f);
            glDrawElements(GL_LINES, (GLsizei)constrain_indices, GL_UNSIGNED_INT, (GLuint*)0);
            glLineWidth(thin);

            // oops
            gfx.vao_main.Bind();
            //gfx.ibo_delabella.Bind();
        }

        if (show_d)
        {
            gfx.SetColor(1.0f, 0.0f, 0.0f, 1.0f);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_TRIANGLES, (GLsizei)cull/*tris_delabella*/ * 3, GL_UNSIGNED_INT, 0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if (show_b)
        {
            gfx.SetColor(0.0f, 0.0f, 1.0f, 1.0f);
            glLineWidth(thick);
            glDrawElements(GL_LINE_LOOP, (GLsizei)contour, GL_UNSIGNED_INT, (GLuint*)0 + (intptr_t)tris_delabella * 3);
            glLineWidth(thin);
        }

        // compare with CDT
        #ifdef WITH_CDT
        if (show_x)
        {
            gfx.vao_cdt.Bind();

            MyIndex tris_cdt = (MyIndex)cdt.triangles.size();

            //glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE);

            gfx.SetColor(0.0f,0.0f,1.0f,1.0f);
            glLineWidth(thin);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_TRIANGLES, /*0,points-1,*/ (GLsizei)tris_cdt * 3, GL_UNSIGNED_INT, 0);

            //glDisable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        }
        #endif

        // voronoi!
        #ifdef VORONOI
        if (show_v)
        {
            gfx.vao_voronoi.Bind();

            // voro-verts in back
            gfx.SetColor(1.0f, 1.0f, 0.0f, 1.0f);
            glPointSize(blob);
            glDrawArrays(GL_POINTS, 0, (GLsizei)(voronoi_vertices - contour));
            glPointSize(dot);
            gfx.ibo_voronoi.Bind();

            gfx.SetColor(0.0f, 0.75f, 0.0f, 1.0f);

            #ifdef VORONOI_POLYS
            // draw structured polys, note: open polys are silently closed (at infinity)
            // glDrawElements(GL_LINE_LOOP, voronoi_indices, GL_UNSIGNED_INT, (GLuint*)0);
            // if you wanna be a ganan: after first M closed polys switch from line_loops to line_strips
            // and draw remaining N open polygons
            glDrawElements(GL_LINE_LOOP, (GLsizei)voronoi_closed_indices, GL_UNSIGNED_INT, (GLuint*)0);
            glDrawElements(GL_LINE_STRIP, (GLsizei)(voronoi_indices - voronoi_closed_indices), GL_UNSIGNED_INT, (GLuint*)0 + (intptr_t)voronoi_closed_indices);
            #else
            // draw edge soup
            glDrawElements(GL_LINES, (GLsizei)voro_cull/*voronoi_indices*/, GL_UNSIGNED_INT, (GLuint*)0);
            #endif
        }
        #endif

        // put verts over everything else
        gfx.vao_main.Bind();
        gfx.SetColor(1.0f, 1.0f, 0.0f, 1.0f);
        glPointSize(blob);
        glDrawArrays(GL_POINTS, 0, (GLsizei)points);
        glPointSize(dot);

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

#ifdef BENCH

int main(int argc, char* argv[])
{
    setlocale(LC_ALL, "");

    const char* test_dist[] = { "uni","std","gam","sym","cir","hex",0};
    const char* test_bias[] = { "","+", 0 };

    int test_size[] =
    {
        1000,2500,5000,
        10000,25000,50000, 
        100000,250000,500000,
        1000000, /*2500000,5000000,
        10000000,25000000,50000000,*/
        0
    };

    const int players = 5;
    Bench bench[players];

    char bin[2] = "";
    char num[16];

    char test_path[100];

    for (int d = 0; test_dist[d]; d++)
    {
        for (int b = 0; test_bias[b]; b++)
        {
            sprintf(test_path, "bench_%s%s.txt", test_dist[d], test_bias[b]);
            FILE* bench_file = fopen(test_path, "w");

            fprintf(bench_file, "        DLB[us]     CDT[us]     FAD[us]   FADMT[us]     DEL[us]\n");

            char* args[] = 
            { 
                bin,
                num, 
                (char*)(bench + 0), 
                (char*)(bench + 1), 
                (char*)(bench + 2), 
                (char*)(bench + 3),
                (char*)(bench + 4),
                (char*)test_dist[d], 
                (char*)test_bias[b]
            };

            for (int i = 0; test_size[i]; i++)
            {
                Bench accum[players];
                memset(accum, 0, sizeof(accum));

                sprintf(num,"%d",test_size[i]);

                uint64_t t0 = uSec();
                int acc = 0;
                do
                {
                    memset(bench, 0, sizeof(bench));
                    bench_main(2, args);
                    for (int i=0; i< players; i++)
                        accum[i] += bench[i];
                    acc++;
                } while (uSec() - t0 < 5000000);

                for (int i = 0; i < players; i++)
                    accum[i] /= acc;

                struct F
                {
                    const char* operator () (int n, uint64_t v)
                    {
                        int len = sprintf(buf, "%llu", (long long unsigned int)v);
                        if (len <= 0)
                            return 0;
                        int sep = (len - 1) / 3;

                        int len2 = len + sep;
                        int pad = n > len2 ? n - len2 : 0;

                        int ofs = pad + len2;

                        if (ofs >= sizeof(buf))
                            return 0;

                        buf[ofs] = 0;
                        ofs--;

                        for (int i = len - 1, j = 0; i >= 0; i--, j++)
                        {
                            if (j == 3)
                            {
                                buf[ofs] = ',';
                                ofs--;
                                j = 0;
                            }

                            buf[ofs] = buf[i];
                            ofs--;
                        }

                        for (int i = 0; i < pad; i++)
                            buf[i] = ' ';

                        return buf;
                    }

                    char buf[32];
                } f[players] = {0};

                // write bench results...
                fprintf(bench_file, "N=%s\n", f[0](0, test_size[i]));
                fprintf(bench_file, "RD: %s %s %s %s %s\n", 
                    f[0](11, accum[0].removing_dups), 
                    f[1](11, accum[1].removing_dups), 
                    accum[2].removing_dups ? "    INVALID" : f[2](11, accum[2].removing_dups), 
                    accum[3].removing_dups ? "    INVALID" : f[3](11, accum[3].removing_dups), 
                    accum[4].removing_dups ? "    INVALID" : f[4](11, accum[4].removing_dups));
                fprintf(bench_file, "TR: %s %s %s %s %s\n", 
                    f[0](11, accum[0].triangulation), 
                    f[1](11, accum[1].triangulation), 
                    f[2](11, accum[2].triangulation), 
                    f[3](11, accum[3].triangulation), 
                    f[4](11, accum[4].triangulation));
                fprintf(bench_file, "CE: %s %s %s %s %s\n", 
                    f[0](11, accum[0].constrain_edges), 
                    f[1](11, accum[1].constrain_edges), 
                    f[2](11, accum[2].constrain_edges), 
                    f[3](11, accum[3].constrain_edges), 
                    f[4](11, accum[4].constrain_edges));
                fprintf(bench_file, "ES: %s %s %s %s %s\n", 
                    f[0](11, accum[0].erase_super), 
                    f[1](11, accum[1].erase_super),
                    f[2](11, accum[2].erase_super),
                    f[3](11, accum[3].erase_super),
                    f[4](11, accum[4].erase_super));
                fprintf(bench_file, "FF: %s %s %s %s %s\n", 
                    f[0](11, accum[0].flood_fill), 
                    f[1](11, accum[1].flood_fill), 
                    f[2](11, accum[2].flood_fill), 
                    f[3](11, accum[3].flood_fill), 
                    f[4](11, accum[4].flood_fill));
                fprintf(bench_file, "PL: %s %s %s %s %s\n", 
                    f[0](11, accum[0].polygons), 
                    f[1](11, accum[1].polygons),
                    f[2](11, accum[2].polygons),
                    f[3](11, accum[3].polygons),
                    f[4](11, accum[4].polygons));
                fprintf(bench_file, "\n");

                // live view
                fflush(bench_file);
            }

            fclose(bench_file);
        }
    }

    return 0;
}
#endif