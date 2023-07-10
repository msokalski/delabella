#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cmath> // INFINITY

#include "delabella.h"

struct dba_point
{
    double x, y;
};

struct dba_edge
{
    int a, b;
};

int errlog(void* stream, const char* fmt, ...)
{
    va_list arg;
    va_start(arg,fmt);
    int ret = vfprintf((FILE*)stream, fmt, arg);
    va_end(arg);
    //fflush((FILE*)stream);
    return ret;
}

int main( int argc, char* argv[] )
{
    FILE *fp = NULL;

    if (argc<2)
    {
        printf( "No input given.\n" );
        exit( 1 );
    }

    fp = fopen( argv[1], "r" );

    if ( !fp )
    {
        printf( "File not found.\n" );
        exit( 1 );
    }

    // Read input file.
    int npt;
    fscanf( fp, "%d", &npt );

    dba_point* cloud = new dba_point[npt];

    for ( int j = 0; j < npt; j++ )
    {
        int indx;
        fscanf( fp, "%d %lf %lf", &indx, &(cloud[ j ].x), &(cloud[ j ].y) );
    }

    int nedg;
    fscanf( fp, "%d", &nedg );

    dba_edge* bounds = new dba_edge[nedg];

    for ( int i = 0 ; i < nedg ; i++ )
    {
        int indx;
        fscanf( fp, "%d %d %d", &indx, &(bounds[i].a), &(bounds[i].b) );
    }
    fclose( fp );

    // Echo input file to verify read success.
    printf( "%d\n", npt );
    for ( int j = 0; j < npt; j++ )
    {
        printf( "%d %.18e %.18e\n", j, cloud[ j ].x, cloud[ j ].y );
    }
    printf( "%d\n", nedg );
    for ( int i = 0 ; i < ( int )nedg ; i++ )
    {
        printf( "%d %d %d\n", i, bounds[i].a, bounds[i].b );
    }
    printf( "\n" );

    FILE* svg = 0;
    double svg_scale = 1;
    double svg_trans_x = 0;
    double svg_trans_y = 0;

    if (argc>3)    
        svg = fopen(argv[3],"w");

    if (svg)
    {
        double x_min = INFINITY, x_max = -INFINITY;
        double y_min = INFINITY, y_max = -INFINITY;
        for ( int j = 0; j < npt; j++ )
        {
            if (cloud[j].x < x_min)
                x_min = cloud[j].x;
            if (cloud[j].x > x_max)
                x_max = cloud[j].x;
            if (cloud[j].y < y_min)
                y_min = cloud[j].y;
            if (cloud[j].y > y_max)
                y_max = cloud[j].y;
        }

        double x_size = x_max - x_min;
        double y_size = y_max - y_min;
        double in_size = x_size > y_size ? x_size : y_size;
        double out_size = 640;

        svg_scale = out_size / (in_size * 1.2);
        svg_trans_x = out_size * 0.5 - (x_max + x_min) * svg_scale * 0.5;
        svg_trans_y = out_size * 0.5 - (y_max + y_min) * svg_scale * 0.5;

        char svg_size[200];
        sprintf(svg_size,"width=\"%f\" height=\"%f\" viewBox=\"%f %f %f %f\"", out_size, out_size, 0.0, 0.0, out_size, out_size);
        const char* svg_attr = "fill=\"none\" stroke=\"currentcolor\" stroke-linecap=\"round\" stroke-linejoin=\"round\" stroke-width=\"2\"";

        fprintf(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" %s %s>\n", svg_size, svg_attr);

        fprintf(svg, "  <rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" rx=\"15\" fill=\"#fff\"/>\n", 0.0, 0.0, out_size, out_size);

        //fprintf(svg, "</svg>\n");
        //fclose(svg);
    }

    // Process points
    IDelaBella2 < double > * idb = IDelaBella2 < double > ::Create();
    idb->SetErrLog( errlog, stdout );

    printf( "Triangulate\n" );

    int tris = 0;
    int verts = idb->Triangulate( npt, &cloud->x, &cloud->y, sizeof( dba_point ) );

    if ( verts > 0 )
    {
        printf( "ConstrainEdges\n" );

        idb->ConstrainEdges( nedg, &bounds->a, &bounds->b, sizeof( dba_edge ) );

        printf( "Removing constrained ears\n" );

        /*
        // concyclic checking ...
        int polys = idb->GetNumPolygons();
        assert(polys == idb->Polygonize());
        */

        tris = idb->FloodFill(false, 0);      

        printf( "Done\n" );
  
        const IDelaBella2<double>::Simplex* dela = idb->GetFirstDelaunaySimplex();

        struct Face
        {
            int v[3];
            static int cmp(const void* a, const void* b)
            {
                const Face* l = (const Face*)a;
                const Face* r = (const Face*)b;
                for (int i=0; i<3; i++)
                {
                    int q = l->v[i] - r->v[i];
                    if (q)
                        return q;
                }
                return 0;
            }
        };

        assert(tris <= 1000);
        Face t[1000];

        for ( int i = 0; i < tris; i++ )
        {
            Face face;

            for ( int j = 0; j < 3; j++ )
            {
                if ( dela->v[ j ]->i >= npt || dela->v[ j ]->i < 0 )
                {
                    printf( "Invalid index.\n" );
                }

                // flip CW->CCW
                // face.v[2-j] = dela->v[j]->i;

                // no flip
                face.v[j] = dela->v[j]->i;
            }

            // make lowest index first, keep direction

            int min_i = 0;
            if (face.v[1] < face.v[min_i])
                min_i = 1;
            if (face.v[2] < face.v[min_i])
                min_i = 2;

            t[i].v[0] = face.v[(min_i+0)%3];
            t[i].v[1] = face.v[(min_i+1)%3];
            t[i].v[2] = face.v[(min_i+2)%3];

            dela = dela->next;
        }

        // sort faces by smaller indices
        qsort(t, tris, sizeof(Face), Face::cmp);

        FILE* o = 0;
        if (argc>2)
        {
            o = fopen(argv[2],"w");
            if (!o)
                printf("Can't open out file\n");
        }

        printf("%d\n", tris);
        if (o)
            fprintf(o,"%d\n", tris);


        int pal[6][3] = 
        {
            {255,0,0}, {255,255,0}, {0,255,0},
            {0,255,255}, {0,0,255}, {255,0,255},
        };

        for (int i=0; i<tris; i++)
        {
            printf("%d %d %d %d\n", i, t[i].v[0],t[i].v[1],t[i].v[2]);
            if (o)
                fprintf(o,"%d %d %d %d\n", i, t[i].v[0],t[i].v[1],t[i].v[2]);

            if (svg)
            {
                fprintf(svg, "  <path stroke-width=\"0.25\" fill=\"RGB(%d,%d,%d)\" fill-opacity=\"0.25\" d=\"", pal[i%6][0], pal[i%6][1], pal[i%6][2]);
                fprintf(svg, "M %f %f L %f %f  L %f %f  Z", 
                    cloud[t[i].v[0]].x*svg_scale + svg_trans_x, cloud[t[i].v[0]].y*svg_scale + svg_trans_y,
                    cloud[t[i].v[1]].x*svg_scale + svg_trans_x, cloud[t[i].v[1]].y*svg_scale + svg_trans_y,
                    cloud[t[i].v[2]].x*svg_scale + svg_trans_x, cloud[t[i].v[2]].y*svg_scale + svg_trans_y);
                fprintf(svg, "\" />\n");
            }
        }

        if (o)
            fclose(o);
    }
    else
    {
        printf( "DLB Error! %d\n", verts );
    }

    if (svg)
    {
        // overlay input
        for (int j=0; j<npt; j++)
        {
            fprintf(svg, "  <circle cx=\"%f\" cy=\"%f\" r=\"3\" />\n", cloud[j].x*svg_scale + svg_trans_x, cloud[j].y*svg_scale + svg_trans_y);
        }

        if (nedg)
        {
            fprintf(svg, "  <path stroke-width=\"2\" d=\"");
            for (int i=0; i<nedg; i++)
            {
                fprintf(svg, "M %f %f L %f %f ", 
                    cloud[bounds[i].a].x*svg_scale + svg_trans_x, cloud[bounds[i].a].y*svg_scale + svg_trans_y,
                    cloud[bounds[i].b].x*svg_scale + svg_trans_x, cloud[bounds[i].b].y*svg_scale + svg_trans_y);
            }

            fprintf(svg, "\" />\n");
        }

        for (int j=0; j<npt; j++)
        {
            fprintf(svg, "  <text stroke=\"none\" fill=\"#800\" font-size=\"16\" font-family=\"Arial, Helvetica, sans-serif\" x=\"%f\" y=\"%f\">%d</text>\n", cloud[j].x*svg_scale + svg_trans_x - 5, cloud[j].y*svg_scale + svg_trans_y - 4, j);
        }


        fprintf(svg, "  <text stroke=\"none\" fill=\"#000\" font-size=\"24\" font-family=\"Arial, Helvetica, sans-serif\" x=\"%f\" y=\"%f\">%d %s</text>\n", 24.0, 24.0, tris, "triangles");

        fprintf(svg, "</svg>\n");
        fclose(svg);        
    }

    delete[] cloud;
    delete[] bounds;

    idb->Destroy();

    return 0;
}
