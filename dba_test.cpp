#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>

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


    // Process points
    IDelaBella2 < double > * idb = IDelaBella2 < double > ::Create();
    idb->SetErrLog( errlog, stdout );

    printf( "Triangulate\n" );

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

        int tris = idb->FloodFill(false, 0);      

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

        for (int i=0; i<tris; i++)
        {
            printf("%d %d %d %d\n", i, t[i].v[0],t[i].v[1],t[i].v[2]);
            if (o)
                fprintf(o,"%d %d %d %d\n", i, t[i].v[0],t[i].v[1],t[i].v[2]);
        }

        if (o)
            fclose(o);
    }
    else
    {
        printf( "DLB Error! %d\n", verts );
    }

    delete[] cloud;
    delete[] bounds;

    idb->Destroy();

    return 0;
}
