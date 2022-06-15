/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#define _CRT_SECURE_NO_WARNINGS

#define DELAUNATOR

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

#include <GL/gl.h>

#ifndef _WIN32
#include <GL/glext.h>
#else
#pragma comment(lib,"OpenGL32.lib")
#pragma comment(lib,"SDL2.lib")
#undef main // on windows SDL does this weird thing
#endif

// competitor
#ifdef DELAUNATOR
#include "delaunator/delaunator-header-only.hpp"
#endif

PFNGLBINDBUFFERPROC glBindBuffer = 0;
PFNGLDELETEBUFFERSPROC glDeleteBuffers = 0;
PFNGLGENBUFFERSPROC glGenBuffers = 0;
PFNGLBUFFERDATAPROC glBufferData = 0;
PFNGLBUFFERSUBDATAPROC glBufferSubData = 0;
PFNGLMAPBUFFERPROC glMapBuffer = 0;
PFNGLUNMAPBUFFERPROC glUnmapBuffer = 0;
PFNGLPRIMITIVERESTARTINDEXPROC glPrimitiveRestartIndex = 0;
bool BindGL(bool prim_restart)
{
	#define BINDGL(proc) if ((*(void**)&proc = SDL_GL_GetProcAddress(#proc)) == 0) return false;
	BINDGL(glBindBuffer);
	BINDGL(glDeleteBuffers);
	BINDGL(glGenBuffers);
	BINDGL(glBufferData);
	BINDGL(glBufferSubData);
	BINDGL(glMapBuffer);
	BINDGL(glUnmapBuffer);
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
        Unmap();
        glDeleteBuffers(1,&buf);
    }

    void* Map()
    {
        if (map)
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
        if (!map)
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

	struct MyPoint
	{
		/*long*/ double x;
		/*long*/ double y;
	};
	
	std::vector<MyPoint> cloud;

	FILE* f = fopen(argv[1],"r");
    int n = atoi(argv[1]);
	if (!f)
	{
        if (n<=2)
        {
		    printf("can't open %s file, terminating!\n", argv[1]);
            return -1;
        }
        printf("generating random %d points\n", n);
       	std::random_device rd{};
    	std::mt19937_64 gen{rd()};

        //std::uniform_real_distribution</*long*/ double> d(-1.0,1.0);
        //std::normal_distribution</*long*/ double> d{0.0,2.0};
        std::gamma_distribution</*long*/ double> d(0.1,2.0);

        for (int i=0; i<n; i++)
        {
            MyPoint p = { d(gen) - 5.0, d(gen) - 5.0 };
			//MyPoint p = { d(gen), d(gen) };
            cloud.push_back(p);
        }
	}
    else
    {
        while (1)
        {
            /*
            long double x,y;
            if (fscanf(f,"%Lf %Lf", &x, &y) != 2)
                break;
            */
            double x,y;
            if (fscanf(f,"%lf %lf", &x, &y) != 2)
                break;
            MyPoint p = {x,y};
            cloud.push_back(p);
        }
        fclose(f);
    }

    int points = (int)cloud.size();

    #ifdef DELAUNATOR
    std::vector<double> coords;
    for (int i=0; i<points; i++)
    {
        coords.push_back(cloud[i].x);
        coords.push_back(cloud[i].y);
    }
    uint64_t t0 = uSec();
    printf("running delaunator...\n");
    delaunator::Delaunator d(coords);
    int tris_delaunator = d.triangles.size() / 3;
    uint64_t t1 = uSec();
    printf("elapsed %d ms\n", (int)((t1-t0)/1000));
    printf("delaunator triangles: %d\n", tris_delaunator);
    /*
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) 
    {
        printf("%d %d %d\n", (int)d.triangles[i], (int)d.triangles[i+1], (int)d.triangles[i+2]);
    }
    */
    #endif

	#ifdef CRUDE_XA
	xa_pool_alloc(1000);
	#endif

	IDelaBella* idb = IDelaBella::Create();
	idb->SetErrLog(errlog, stdout);
	
    printf("running delabella...\n");
    uint64_t t2 = uSec();
	int verts = idb->Triangulate(points, &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
	int tris_delabella = verts / 3;
    int contour = idb->GetNumBoundaryVerts();
    int non_contour = idb->GetNumInternalVerts();
	int vert_num = contour + non_contour;

    uint64_t t3 = uSec();
    printf("elapsed %d ms\n", (int)((t3-t2)/1000));
    printf("delabella triangles: %d\n", tris_delabella);
    printf("delabella contour: %d\n", contour);
    printf("delabella is %s\n", tris_delabella == 2* vert_num - 4 - (contour-2) ? "CORRECT" : "WRONG");


	// if positive, all ok 
	if (verts<=0)
	{
		// no points given or all points are colinear
		// make emergency call ...
		idb->Destroy();
		return -2;
	}

    if (argc>=3)
    {
        f = fopen(argv[2],"w");
        if (f)
        {
            for (int i=0; i<tris_delabella; i++)
            {
                const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
                fprintf(f,"%d %d %d\n", 
                    dela->v[0]->i,
                    dela->v[1]->i,
                    dela->v[2]->i);
                dela = dela->next;
            }
            fclose(f);
        }
    }

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
    #ifdef CRUDE_XA
    const char* title = "delablella-xa";
    #else
    const char* title = "delablella";
    #endif
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

	if (!BindGL(prim_restart))
	{
		printf("Can't bind to necessary GL functions, terminating!\n");
		idb->Destroy();
		return -1;
	}

	printf("preparing graphics...\n");

	typedef GLfloat gl_t;
	GLenum gl_e = GL_FLOAT;
	#define glLoadMatrix(m) glLoadMatrixf(m)

	// useless until glVertexAttribLPointer and accompanying shader using dmat/dvec/double arithmetic is being used
	// typedef GLdouble gl_t;
	// GLenum gl_e = GL_DOUBLE;
	// #define glLoadMatrix(m) glLoadMatrixd(m)


    // create vbo and ibo
    Buf vbo, ibo_delabella;

	vbo.Gen(GL_ARRAY_BUFFER, sizeof(gl_t[3]) * points);
    gl_t* vbo_ptr = (gl_t*)vbo.Map();
    float box[4]={(float)cloud[0].x, (float)cloud[0].y, (float)cloud[0].x, (float)cloud[0].y};
	for (int i = 0; i<points; i++)
    {
        vbo_ptr[3*i+0] = (gl_t)cloud[i].x;
        vbo_ptr[3*i+1] = (gl_t)cloud[i].y;
        vbo_ptr[3*i+2] = (gl_t)(1.0); // i%5; // color

        box[0] = fmin(box[0], (float)cloud[i].x);
        box[1] = fmin(box[1], (float)cloud[i].y);
        box[2] = fmax(box[2], (float)cloud[i].x);
        box[3] = fmax(box[3], (float)cloud[i].y);
    }
    vbo.Unmap();

	printf("box: x0:%f y0:%f x1:%f y1:%f\n", box[0], box[1], box[2], box[3]);

    // pure indices, without: center points, restarts, loop closing
    // points may be a bit too much (cuza duplicates)
    int voronoi_indices = 2 * (vert_num + tris_delabella - 1) + contour;
    int voronoi_vertices = tris_delabella + contour;

    if (prim_restart)
        voronoi_indices += vert_num; // add primitive restarts
    else
        voronoi_indices = 4 * (vert_num + tris_delabella - 1); // almost 2x bigger ehh

    int ibo_voronoi_idx = 0;
    Buf vbo_voronoi, ibo_voronoi;
	vbo_voronoi.Gen(GL_ARRAY_BUFFER, sizeof(gl_t[3]) * voronoi_vertices);
	ibo_voronoi.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * voronoi_indices);
    gl_t* vbo_voronoi_ptr = (gl_t*)vbo_voronoi.Map();
    GLuint* ibo_voronoi_ptr = (GLuint*)ibo_voronoi.Map();

    ibo_delabella.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delabella + sizeof(GLuint) * contour);	
    GLuint* ibo_ptr = (GLuint*)ibo_delabella.Map();
    const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
	for (int i = 0; i<tris_delabella; i++)
	{
        int v0 = dela->v[0]->i;
        int v1 = dela->v[1]->i;
        int v2 = dela->v[2]->i;

        ibo_ptr[3*i+0] = (GLuint)v0;
        ibo_ptr[3*i+1] = (GLuint)v1;
        ibo_ptr[3*i+2] = (GLuint)v2;

        // put it into vbo_voronoi at 'i'

		// almost exact in hybrid
		/*
		vbo_voronoi_ptr[3 * i + 0] = (gl_t)dela->n.x;
		vbo_voronoi_ptr[3 * i + 1] = (gl_t)dela->n.y;
		vbo_voronoi_ptr[3 * i + 2] = (gl_t)(-2.0*dela->n.z);
		*/

		// less jumping on extreme zooming
        vbo_voronoi_ptr[3*i+0] = (gl_t)(-0.5 * (double)dela->n.x / (double)dela->n.z);
        vbo_voronoi_ptr[3*i+1] = (gl_t)(-0.5 * (double)dela->n.y / (double)dela->n.z);
        vbo_voronoi_ptr[3*i+2] = (gl_t)(1.0);

		dela = dela->next;
	}

	const DelaBella_Vertex* prev = idb->GetFirstBoundaryVertex();
    const DelaBella_Vertex* vert = prev->next;
    int contour_min = points-1;
    int contour_max = 0;
    for (int i = 0; i<contour; i++)    
    {
        ibo_ptr[i + 3*tris_delabella] = (GLuint)vert->i;
        contour_min = vert->i < contour_min ? vert->i : contour_min;
        contour_max = vert->i > contour_max ? vert->i : contour_max;

        // put infinite edge normal to vbo_voronoi at tris_delabella + 'i'

		// coords are overwritten with edge normals in edge postproc!
        vbo_voronoi_ptr[3*(tris_delabella+i)+0] = (gl_t)(prev->x);
        vbo_voronoi_ptr[3*(tris_delabella+i)+1] = (gl_t)(prev->y);
        vbo_voronoi_ptr[3*(tris_delabella+i)+2] = (gl_t)(0.0);

        // create special-fan / line_strip in ibo_voronoi around this boundary vertex
        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)i + tris_delabella; // begin

        // iterate all dela faces around prev
        // add their voro-vert index == dela face index
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = prev->StartIterator(&it); 

        // it starts at random face, so lookup the prev->vert edge
        while (1)
        {
            if (t->n.z<0)
            {
                if (t->v[0] == prev && t->v[1] == vert ||
                    t->v[1] == prev && t->v[2] == vert ||
                    t->v[2] == prev && t->v[0] == vert)
                    break;
            }
            t = it.Next();
        }

        // now iterate around, till we're inside the boundary
        while (t->n.z<0)
        {
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // end
            if (!prim_restart)
                ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // begin of next line segment
            t = it.Next();
        }

        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)( i==0 ? contour-1 : i-1 ) + tris_delabella; // loop-wrapping!
        
        if (prim_restart)
            ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart

        prev = vert;
        vert = vert->next;
    }

    int voronoi_strip_indices = ibo_voronoi_idx;

    // finally, for all internal vertices
    vert = idb->GetFirstInternalVertex();
    for (int i = 0; i<non_contour; i++)
    {
        // create regular-fan / line_loop in ibo_voronoi around this internal vertex
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = vert->StartIterator(&it);
        const DelaBella_Triangle* e = t;
        do
        {
            assert(t->n.z<0);
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // begin
            t = it.Next();
            if (!prim_restart)
                ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // end
        } while (t!=e);
        
        if (prim_restart)
            ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart

        vert = vert->next;
    }

    int voronoi_loop_indices = ibo_voronoi_idx - voronoi_strip_indices;

    // TODO: FIX CONTOUR BUILDING FOR TRIVIAL 3-POINTS CASE !!!
    assert(ibo_voronoi_idx == voronoi_indices);

    ibo_delabella.Unmap();
    vbo_voronoi.Unmap();
    ibo_voronoi.Unmap();

    #ifdef DELAUNATOR
    Buf ibo_delaunator;
    ibo_delaunator.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delaunator);
    ibo_ptr = (GLuint*)ibo_delaunator.Map();
    for (int i = 0; i<tris_delaunator; i++)    
    {
        ibo_ptr[3*i+0] = (GLuint)d.triangles[3*i+0];
        ibo_ptr[3*i+1] = (GLuint)d.triangles[3*i+1];
        ibo_ptr[3*i+2] = (GLuint)d.triangles[3*i+2];
    }
    ibo_delaunator.Unmap();
    #endif

    // now, everything is copied to gl, free delabella
	#ifdef CRUDE_XA
	// close pool before destroy so we won't move 
	// everything unneccessarily to the pool
	xa_pool_free(); 
	#endif

    idb->Destroy();

    int vpw, vph;
    SDL_GL_GetDrawableSize(window, &vpw, &vph);

    double cx = 0.5 * (box[0]+box[2]);
    double cy = 0.5 * (box[1]+box[3]);
    double scale = 2.0 * fmin((double)vpw/(box[2]-box[0]),(double)vph/(box[3]-box[1]));
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

                case SDL_KEYUP:
                    if( event.key.keysym.sym == SDLK_ESCAPE )
                        return 0;
                    break;
            }

            if (x)
                break;
        }

        if (x)
            break;

        int vpw, vph;
        SDL_GL_GetDrawableSize(window, &vpw, &vph);
        glViewport(0,0,vpw,vph);

        double scale = pow(1.01, zoom);

        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(cx - vpw/scale, cx + vpw/scale, cy - vph/scale, cy + vph/scale, -1, +1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glScalef(1,1,0); // flatten z when not shaded

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);

        // hello ancient world
        vbo.Bind();
        //glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)
		glVertexPointer(3, gl_e, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

        ibo_delabella.Bind();

        glColor4f(0.2f,0.2f,0.2f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delabella*3, GL_UNSIGNED_INT, 0);

        // put verts over fill
        glColor4f(1.0f,1.0f,0.0f,1.0f);
        glPointSize(3.0f);
        glDrawArrays(GL_POINTS, 0, points);
        glPointSize(1.0f);

        glColor4f(0.5f,0.5f,0.5f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delabella * 3, GL_UNSIGNED_INT, 0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glColor4f(1.0f,0.0f,0.0f,1.0f);
        //glLineWidth(3.0f);
        glDrawElements(GL_LINE_LOOP, /*contour_min, contour_max,*/ contour, GL_UNSIGNED_INT, (GLuint*)0 + tris_delabella*3);
        //glLineWidth(1.0f);

        // delaunator
		#if 0
        #ifdef DELAUNATOR
        ibo_delaunator.Bind();
        glColor4f(1.0f,1.0f,1.0f,1.0f);
        glLineWidth(1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delaunator * 3, GL_UNSIGNED_INT, 0);
        #endif
		#endif

        // voronoi!
        vbo_voronoi.Bind();
        //glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)
		glVertexPointer(3, gl_e, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

        const static gl_t z2w[16]=
        {
            1,0,0,0,
            0,1,0,0,
            0,0,0,1,
            0,0,0,0
        };

        glLoadMatrix(z2w);

		// voro-verts in back
		glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
		glPointSize(3.0f);
		glDrawArrays(GL_POINTS, 0, tris_delabella);
		glPointSize(1.0f);
        ibo_voronoi.Bind();

        glColor4f(0.0f,0.75f,0.0f,1.0f);

        if (prim_restart)
        {
            // first, draw open cells
            glDrawElements(GL_LINE_STRIP, voronoi_strip_indices, GL_UNSIGNED_INT, (GLuint*)0);

            // then closed cells
            glDrawElements(GL_LINE_LOOP, voronoi_loop_indices, GL_UNSIGNED_INT, (GLuint*)0 + voronoi_strip_indices);
        }
        else
        {
            // first, draw open cells
            glDrawElements(GL_LINES, voronoi_strip_indices, GL_UNSIGNED_INT, (GLuint*)0);

            // then closed cells
            glDrawElements(GL_LINES, voronoi_loop_indices, GL_UNSIGNED_INT, (GLuint*)0 + voronoi_strip_indices);
        }

        SDL_GL_SwapWindow(window);
        SDL_Delay(15);
    }

    vbo.Del();
    ibo_delabella.Del();

    #ifdef DELAUNATOR
    ibo_delaunator.Del();
    #endif

    vbo_voronoi.Del();
    ibo_voronoi.Del();

    SDL_GL_DeleteContext( context );
    SDL_DestroyWindow( window );
    SDL_Quit();

    #ifdef XA_VAL_LEAKS
    printf("LEAKED %d allocs\n", xa_leaks(0));
    #endif

	printf("exiting!\n");

	return 0;
}