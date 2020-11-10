#include<STLtools.H>

namespace STLtools
{
    AMREX_GPU_DEVICE_MANAGED Real *tri_pts=NULL;
    AMREX_GPU_DEVICE_MANAGED Real *tri_normals=NULL;
    Gpu::ManagedVector<Real>* tri_pts_vec=NULL; 
    Gpu::ManagedVector<Real>* tri_normals_vec=NULL; 
    AMREX_GPU_DEVICE_MANAGED int num_tri=0;

    AMREX_GPU_DEVICE_MANAGED int ndata_per_tri=9;
    AMREX_GPU_DEVICE_MANAGED int ndata_per_normal=3;
    int nlines_per_facet=7;

    void read_stl_file(std::string fname)
    {
        std::string tmpline,tmp1,tmp2;
        int nlines=0;

        std::ifstream infile(fname.c_str());
        Print()<<"STL file name:"<<fname<<"\n";

        std::getline(infile,tmpline); //solid <solidname>
        while(!infile.eof())
        {
            std::getline(infile,tmpline);
            if(tmpline.find("endsolid")!=std::string::npos)
            {
                break;
            }
            nlines++;
        }

        if(nlines%nlines_per_facet!=0)
        {
            amrex::Abort("may be there are blank lines in the STL file\n");
        }

        num_tri=nlines/nlines_per_facet;
        Print()<<"number of triangles:"<<num_tri<<"\n";

        tri_pts_vec=new Gpu::ManagedVector<Real>;
        tri_normals_vec=new Gpu::ManagedVector<Real>;

        tri_pts_vec->resize(num_tri*ndata_per_tri);
        tri_normals_vec->resize(num_tri*ndata_per_normal);

        infile.seekg(0);
        std::getline(infile,tmpline); //solid <solidname>

        for(int i=0;i<num_tri;i++)
        {
            std::getline(infile,tmpline);  //facet normal
            std::istringstream fcnormal(tmpline);
            fcnormal>>tmp1>>tmp2
                >>(*tri_normals_vec)[i*ndata_per_normal+0]
                >>(*tri_normals_vec)[i*ndata_per_normal+1]
                >>(*tri_normals_vec)[i*ndata_per_normal+2];

            std::getline(infile,tmpline); // outer loop

            std::getline(infile,tmpline); //vertex 1
            std::istringstream vertex1(tmpline);
            vertex1>>tmp1
                >>(*tri_pts_vec)[i*ndata_per_tri+0]
                >>(*tri_pts_vec)[i*ndata_per_tri+1]
                >>(*tri_pts_vec)[i*ndata_per_tri+2];

            std::getline(infile,tmpline); //vertex 2
            std::istringstream vertex2(tmpline);
            vertex2>>tmp1 
                >>(*tri_pts_vec)[i*ndata_per_tri+3]
                >>(*tri_pts_vec)[i*ndata_per_tri+4]
                >>(*tri_pts_vec)[i*ndata_per_tri+5];

            std::getline(infile,tmpline); //vertex 3
            std::istringstream vertex3(tmpline);
            vertex3>>tmp1 //vertex
                >>(*tri_pts_vec)[i*ndata_per_tri+6]
                >>(*tri_pts_vec)[i*ndata_per_tri+7]
                >>(*tri_pts_vec)[i*ndata_per_tri+8];

            std::getline(infile,tmpline); //end loop
            std::getline(infile,tmpline); //end facet
        }

        tri_pts     = tri_pts_vec->dataPtr();
        tri_normals = tri_normals_vec->dataPtr();

        /*for(int i=0;i<num_tri;i++)
        {
            Print()<<"Normals:"
                <<tri_normals[i*ndata_per_normal+0]<<"\t"
                <<tri_normals[i*ndata_per_normal+1]<<"\t"
                <<tri_normals[i*ndata_per_normal+2]<<"\n";

            for(int j=0;j<3;j++)
            {
                Print()<<"point "<<j<<" :"
                    <<tri_pts[i*ndata_per_tri+3*j+0]<<"\t"
                    <<tri_pts[i*ndata_per_tri+3*j+1]<<"\t"
                    <<tri_pts[i*ndata_per_tri+3*j+2]<<"\n";
            }
        }*/
    }
}
