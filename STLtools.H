#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include<AMReX_REAL.H>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

using namespace amrex;
namespace STLtools
{
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real *tri_pts; 
    extern Gpu::ManagedVector<Real>* tri_pts_vec; 
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real *tri_normals; 
    extern Gpu::ManagedVector<Real>* tri_normals_vec; 
    extern AMREX_GPU_DEVICE_MANAGED int num_tri; 

    extern AMREX_GPU_DEVICE_MANAGED int ndata_per_tri;
    extern AMREX_GPU_DEVICE_MANAGED int ndata_per_normal;

    extern int nlines_per_facet;
    void read_stl_file(std::string fname);

    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real Distance2(Real P1[3],Real P2[3])
    {
        return( (P1[0]-P2[0])*(P1[0]-P2[0]) +
                (P1[1]-P2[1])*(P1[1]-P2[1]) +
                (P1[2]-P2[2])*(P1[2]-P2[2]) );
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real DotProd(Real v1[3],Real v2[3])
    {
        return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real side_op(Real L1[6],Real L2[6])
    {
        return(    L1[0]*L2[4] 
                +  L1[1]*L2[5] 
                +  L1[2]*L2[3] 
                +  L1[3]*L2[2] 
                +  L1[4]*L2[0] 
                +  L1[5]*L2[1] );
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE void getvec(Real P1[3],Real P2[3],Real v[3])
    {
        v[0]=P2[0]-P1[0];
        v[1]=P2[1]-P1[1];
        v[2]=P2[2]-P1[2];
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE void CrossProd(Real v1[3],Real v2[3],Real v[3])
    {
        v[0]=v1[1]*v2[2]-v1[2]*v2[1];
        v[1]=v1[2]*v2[0]-v1[0]*v2[2];
        v[2]=v1[0]*v2[1]-v1[1]*v2[0];
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE void get_plucker_coords(Real v1[3],Real v2[3],Real L[6])
    {
        L[0] = v1[0]*v2[1] - v1[1]*v2[0];
        L[1] = v1[0]*v2[2] - v1[2]*v2[0];
        L[2] = v1[0]       - v2[0];
        L[3] = v1[1]*v2[2] - v1[2]*v2[1];
        L[4] = v1[2]       - v2[2];
        L[5] = v2[1]       - v1[1];
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE void side_op3(Real v1[3],Real v2[3],
            Real t1[3],Real t2[3],Real t3[3],
            Real &S1, Real &S2, Real &S3)
    {

        Real L[6],e1[6],e2[6],e3[6];

        get_plucker_coords(v1,v2,L);
        get_plucker_coords(t1,t2,e1);
        get_plucker_coords(t2,t3,e2);
        get_plucker_coords(t3,t1,e3);

        S1=side_op(L,e1);
        S2=side_op(L,e2);
        S3=side_op(L,e3);
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real triangle_area(Real P1[3],Real P2[3],Real P3[3])
    {
        Real v1[3],v2[3],area[3];

        getvec(P1,P2,v1);
        getvec(P1,P3,v2);
        CrossProd(v1,v2,area);
        return(0.5 * sqrt(DotProd(area,area)) );
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE int lineseg_tri_intersect(Real v1[3],Real v2[3],
            Real t1[3],Real t2[3],Real t3[3])
    {
        //see plucker coordinates based method
        //https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/Line_Triangle.html

        Real S1,S2,S3;
        Real tri_area,area1,area2;
        Real L2[6],L3[6],L4[6],ls_s1,ls_s2;

        side_op3(v1,v2,t1,t2,t3,S1,S2,S3);

        //we are assuming there are no intersections initially
        int no_intersections=1;

        Real eps = std::numeric_limits<Real>::epsilon();

        //coplanar (S1,S2,S3 = 0)
        if(fabs(S1) < eps and fabs(S2) < eps and fabs(S3) < eps) 
        {
            //Print()<<"line segment and triangle are in the same plane\t"<<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
            tri_area=triangle_area(t1,t2,t3);

            /*if(tri_area == 0)
              {
              amrex::Abort("problem with triangle\n");
              }*/
            area1=(triangle_area(t1,t2,v1)+triangle_area(t2,t3,v1)+triangle_area(t3,t1,v1));
            area2=(triangle_area(t1,t2,v2)+triangle_area(t2,t3,v2)+triangle_area(t3,t1,v2));

            if( fabs(area1-tri_area)>eps or fabs(area2-tri_area)>eps )
            {
                no_intersections = 0;
            }
        }
        //proper and edge intersection
        else if( (S1 < 0.0 and S2 < 0.0 and S3 < 0.0) or 
                (S1 > 0.0 and S2 > 0.0 and S3 > 0.0) or
                (fabs(S1) < eps and S2*S3 > 0.0) or     //S1=0
                (fabs(S2) < eps and S3*S1 > 0.0) or     //S2=0
                (fabs(S3) < eps and S1*S2 > 0.0) )      //S3=0
        {

            get_plucker_coords(v1,t1,L2);
            get_plucker_coords(t1,v2,L3);
            get_plucker_coords(t2,t3,L4);

            /*if(fabs(S1*S2*S3) < eps)
              {
              Print()<<"edge intersection S1,S2,S3:"
              <<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
              }*/

            ls_s1 = side_op(L4,L3);
            ls_s2 = side_op(L4,L2);

            if(ls_s1*ls_s2 > 0.0)
            {
                no_intersections = 0;
            }
        }
        //vertex intersection
        else if( (fabs(S1) < eps and fabs(S2) < eps) or  //S1,S2=0
                (fabs(S2) < eps and fabs(S3) < eps) )   //S2,S3=0
        {

            //Print()<<"vertex intersection type 1\n";   
            //don't chose vertex 2 or 3
            get_plucker_coords(v2,t1,L2);
            get_plucker_coords(t1,v1,L3);
            get_plucker_coords(t2,t3,L4);

            ls_s1=side_op(L4,L3);
            ls_s2=side_op(L4,L2);

            if(ls_s1*ls_s2 > 0)
            {
                no_intersections = 0;
            }
        }
        else if(fabs(S3) < eps and fabs(S1) < eps) //S3,S1=0
        {

            //Print()<<"vertex intersection type 2\n";
            //don't chose vertex 1
            get_plucker_coords(v2,t2,L2);
            get_plucker_coords(t2,v1,L3);
            get_plucker_coords(t3,t1,L4);

            ls_s1=side_op(L4,L3);
            ls_s2=side_op(L4,L2);

            if(ls_s1*ls_s2 > 0)
            {
                no_intersections=0;
            }
        }

        return(no_intersections);

    }

}
