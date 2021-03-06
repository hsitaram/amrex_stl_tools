#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_utils.H>
#include <STLtools.H>

using namespace amrex;
AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_tri(Real P1[3],Real P2[3],
        Real t1[3],Real t2[3],Real t3[3])
{
    bool success=false;

    for(int tr=0;tr<STLtools::num_tri;tr++)
    {
        t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
        t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
        t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

        t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
        t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
        t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

        t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
        t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
        t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

        if(STLtools::lineseg_tri_intersect(P1,P2,t1,t2,t3) == 0)
        {
            success=true;
            break;
        }
    }

    return(success);
}


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        int nghost = 1;
        int max_grid_size=64;
        MultiFab lsphi,apx,apy,apz;
        std::string stl_fname;

        Vector<Real> plo;
        Vector<Real> phi;
        Vector<int> ncells;
        Vector<Real> pointoutside;
        Real dx[3];

        ParmParse pp;
        pp.getarr("prob_lo",plo);
        pp.getarr("prob_hi",phi);
        pp.getarr("ncells",ncells);
        pp.get("stl_file",stl_fname);
        pp.getarr("outside_point",pointoutside);
        pp.query("max_grid_size",max_grid_size);

        RealBox real_box({AMREX_D_DECL(plo[0], plo[1], plo[2])},
                {AMREX_D_DECL(phi[0], phi[1], phi[2])});

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));

        dx[0]=(phi[0]-plo[0])/ncells[0];
        dx[1]=(phi[1]-plo[1])/ncells[1];
        dx[2]=(phi[2]-plo[2])/ncells[2];

        Box domain(domain_lo, domain_hi);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);
        DistributionMapping dm(ba);
        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());
        BoxArray facex_ba = amrex::convert(ba, {1,0,0});
        BoxArray facey_ba = amrex::convert(ba, {0,1,0});
        BoxArray facez_ba = amrex::convert(ba, {0,0,1});

        lsphi.define(nodal_ba, dm, 1, nghost);
        apx.define(facex_ba,dm,1,nghost);
        apy.define(facey_ba,dm,1,nghost);
        apz.define(facez_ba,dm,1,nghost);

        STLtools::read_stl_file(stl_fname);

        Array<Real,AMREX_SPACEDIM> plo_arr{plo[0],plo[1],plo[2]};
        Array<Real,AMREX_SPACEDIM> po_arr{pointoutside[0],pointoutside[1],pointoutside[2]};
        Real inside=-1.0;
        Real outside=1.0;

        for (MFIter mfi(lsphi); mfi.isValid(); ++mfi) // Loop over grids
        {
            const Box& bx = mfi.validbox();

            auto lsphi_arr=lsphi[mfi].array();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Real coords[3];
                Real t1[3],t2[3],t3[3];
                Real outp[]={po_arr[0],po_arr[1],po_arr[2]};

                coords[0]=plo_arr[0]+i*dx[0];
                coords[1]=plo_arr[1]+j*dx[1];
                coords[2]=plo_arr[2]+k*dx[2];

                int num_intersects=0;

                for(int tr=0;tr<STLtools::num_tri;tr++)
                {
                    t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
                    t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
                    t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

                    t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
                    t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
                    t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

                    t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
                    t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
                    t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

                    num_intersects += (1-STLtools::lineseg_tri_intersect(outp,coords,t1,t2,t3));
                }
                if(num_intersects%2 == 0)
                {
                    lsphi_arr(i,j,k)=outside;
                }
                else
                {
                    lsphi_arr(i,j,k)=inside;
                }

            }); 
        }

        lsphi.FillBoundary(geom.periodicity());

        #include"faceFrac_x.H"

        //write plot file
        std::string pltfile;
        pltfile = "pltphi";
        WriteSingleLevelPlotfile(pltfile, lsphi, {"phi"}, geom, 0.0, 0);

        pltfile = "pltapx";
        WriteSingleLevelPlotfile(pltfile, apx, {"apx"}, geom, 0.0, 0);

    }

    amrex::Finalize();
}
