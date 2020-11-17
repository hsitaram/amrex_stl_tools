for (MFIter mfi(apx); mfi.isValid(); ++mfi) // Loop over x faces
{
    const Box& bx = mfi.validbox();
    auto lsphi_arr = lsphi[mfi].array();
    auto apx_arr   = apx[mfi].array();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Real t1[3],t2[3],t3[3],ip[3],P1[3],P2[3],P3[3],P4[3];
        Real vals[4],testval;
        Real l1,l2;
        bool success,found_tri;

        apx_arr(i,j,k)=0.0;
        Real sum_ls = lsphi_arr(i,j,k)+lsphi_arr(i,j+1,k)+lsphi_arr(i,j+1,k+1)+lsphi_arr(i,j,k+1);

        if(sum_ls != 4.0*inside and sum_ls != 4.0*outside)
        { 
           P1[0]=plo_arr[0]+i*dx[0];
           P1[1]=plo_arr[1]+j*dx[1];
           P1[2]=plo_arr[2]+k*dx[2];

           P2[0]=plo_arr[0]+i*dx[0];
           P2[1]=plo_arr[1]+(j+1)*dx[1];
           P2[2]=plo_arr[2]+k*dx[2];
           
           P3[0]=plo_arr[0]+i*dx[0];
           P3[1]=plo_arr[1]+(j+1)*dx[1];
           P3[2]=plo_arr[2]+(k+1)*dx[2];
           
           P4[0]=plo_arr[0]+i*dx[0];
           P4[1]=plo_arr[1]+j*dx[1];
           P4[2]=plo_arr[2]+(k+1)*dx[2];

            if(sum_ls == 0)
            {
                if(lsphi_arr(i,j,k)==lsphi_arr(i,j+1,k+1))
                {
                    Print()<<"error multi cuts:"<<lsphi_arr(i,j,k)<<"\t"<<
                        lsphi_arr(i,j+1,k)<<"\t"<<lsphi_arr(i,j+1,k+1)<<"\t"<<
                        lsphi_arr(i,j,k+1)<<"\n";
                }

                //P1 and P2 are inside
                if(lsphi_arr(i,j,k)==inside and lsphi_arr(i,j+1,k)==inside)
                {

                    found_tri=get_tri(P1,P4,t1,t2,t3);
                    success=STLtools::find_intersection_point(P1,P4,t1,t2,t3,ip);
                    l1=fabs(ip[2]-P1[2]);
                    
                    found_tri=get_tri(P2,P3,t1,t2,t3);
                    success=STLtools::find_intersection_point(P2,P3,t1,t2,t3,ip);
                    l2=fabs(ip[2]-P2[2]);

                    apx_arr(i,j,k)=1.0 - 0.5*dx[1]*(l1+l2)/(dx[1]*dx[2]);
                }
                //P2 and P3 are inside
                if(lsphi_arr(i,j+1,k)==inside and lsphi_arr(i,j+1,k+1)==inside)
                {
                    found_tri=get_tri(P2,P1,t1,t2,t3);
                    success=STLtools::find_intersection_point(P2,P1,t1,t2,t3,ip);
                    l1=fabs(ip[1]-P2[1]);

                    found_tri=get_tri(P3,P4,t1,t2,t3);
                    success=STLtools::find_intersection_point(P3,P4,t1,t2,t3,ip);
                    l2=fabs(ip[1]-P3[1]);

                    apx_arr(i,j,k)=1.0 - 0.5*dx[2]*(l1+l2)/(dx[1]*dx[2]);
                }
                //P3 and P4 are inside
                if(lsphi_arr(i,j+1,k+1)==inside and lsphi_arr(i,j,k+1)==inside)
                {

                    found_tri=get_tri(P3,P2,t1,t2,t3);
                    success=STLtools::find_intersection_point(P3,P2,t1,t2,t3,ip);
                    l1=fabs(ip[2]-P3[2]);

                    found_tri=get_tri(P4,P1,t1,t2,t3);
                    success=STLtools::find_intersection_point(P4,P1,t1,t2,t3,ip);
                    l2=fabs(ip[2]-P4[2]);

                    apx_arr(i,j,k)=1.0 - 0.5*dx[1]*(l1+l2)/(dx[1]*dx[2]);
                }
                //P4 and P1 are inside
                if(lsphi_arr(i,j,k+1)==inside and lsphi_arr(i,j,k)==inside)
                {

                    found_tri=get_tri(P4,P3,t1,t2,t3);
                    success=STLtools::find_intersection_point(P4,P3,t1,t2,t3,ip);
                    l1=fabs(ip[1]-P4[1]);
                    
                    found_tri=get_tri(P1,P2,t1,t2,t3);
                    success=STLtools::find_intersection_point(P1,P2,t1,t2,t3,ip);
                    l2=fabs(ip[1]-P1[1]);

                    apx_arr(i,j,k)=1.0 - 0.5*dx[2]*(l1+l2)/(dx[1]*dx[2]);
                }

            }
            if(fabs(sum_ls)==2.0)
            {
                //find the odd-man out
                vals[0]=lsphi_arr(i,j,k);
                vals[1]=lsphi_arr(i,j+1,k);
                vals[2]=lsphi_arr(i,j+1,k+1);
                vals[3]=lsphi_arr(i,j,k+1);

                if(sum_ls < 0.0) //mostly inside
                {
                    testval=outside;
                }
                else
                {
                    testval=inside;
                }
                int c;
                for(c=0;c<4;c++)
                {
                    if((vals[c]+testval) != 0.0)
                    {
                        break;
                    }
                }
                if(c==0) //P1
                {
                    found_tri=get_tri(P2,P1,t1,t2,t3);
                    success=STLtools::find_intersection_point(P2,P1,t1,t2,t3,ip);
                    l1=fabs(ip[1]-P1[1]);
                    
                    found_tri=get_tri(P4,P1,t1,t2,t3);
                    success=STLtools::find_intersection_point(P4,P1,t1,t2,t3,ip);
                    l2=fabs(ip[1]-P1[1]);
                    
                    apx_arr(i,j,k)=0.5*l1*l2/(dx[1]*dx[2]);
                    if(sum_ls > 0.0)
                    {
                       apx_arr(i,j,k)=1.0-apx_arr(i,j,k);
                    }
                }
                if(c==1) //P2
                {
                    found_tri=get_tri(P1,P2,t1,t2,t3);
                    success=STLtools::find_intersection_point(P1,P2,t1,t2,t3,ip);
                    l1=fabs(ip[1]-P2[1]);
                    
                    found_tri=get_tri(P3,P2,t1,t2,t3);
                    success=STLtools::find_intersection_point(P3,P2,t1,t2,t3,ip);
                    l2=fabs(ip[2]-P2[2]);

                    apx_arr(i,j,k)=0.5*l1*l2/(dx[1]*dx[2]);
                    if(sum_ls > 0.0)
                    {
                       apx_arr(i,j,k)=1.0-apx_arr(i,j,k);
                    }
                }
                if(c==2) //P3
                {
                    found_tri=get_tri(P2,P3,t1,t2,t3);
                    success=STLtools::find_intersection_point(P2,P3,t1,t2,t3,ip);
                    l1=fabs(ip[2]-P3[2]);
                    
                    found_tri=get_tri(P4,P3,t1,t2,t3);
                    success=STLtools::find_intersection_point(P4,P3,t1,t2,t3,ip);
                    l2=fabs(ip[1]-P3[1]);
                    
                    apx_arr(i,j,k)=0.5*l1*l2/(dx[1]*dx[2]);
                    if(sum_ls > 0.0)
                    {
                       apx_arr(i,j,k)=1.0-apx_arr(i,j,k);
                    }
                }
                if(c==3) //P4
                {
                    found_tri=get_tri(P1,P4,t1,t2,t3);
                    success=STLtools::find_intersection_point(P1,P4,t1,t2,t3,ip);
                    l1=fabs(ip[2]-P4[2]);
                    
                    found_tri=get_tri(P3,P4,t1,t2,t3);
                    success=STLtools::find_intersection_point(P3,P4,t1,t2,t3,ip);
                    l2=fabs(ip[1]-P4[1]);
                    
                    apx_arr(i,j,k)=0.5*l1*l2/(dx[1]*dx[2]);
                    if(sum_ls > 0.0)
                    {
                       apx_arr(i,j,k)=1.0-apx_arr(i,j,k);
                    }
                }

            }

        }
    });
}
