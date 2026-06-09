  //// DEBUG    
       // debug      
//------------------------------------------------
      Real xx = prob_lo[0] + (i + Real(0.5)) * dx[0]; 
      Real yy = prob_lo[1] + (j + Real(0.5)) * dx[1]; 
      Real zz = prob_lo[2] + (k + Real(0.5)) * dx[2]; 
      Real rad = std::sqrt(xx*xx+yy*yy);

      //bool cellprob = (zz > 0.132) && (zz < 0.138) && (xx < -0.070) && (xx >  -0.076) && (j==64);
      bool cellprob = (i == 3 && j == 64 && k == 107);
      
      //bool cellprob= std::abs(sumnorm)< 1.e-6;
      //cellprob = cellprob || std::abs(sumError) > 1.e-6 || max(std::abs(Err_A[0]), max(std::abs(Err_A[1]), std::abs(Err_A[2]))) > 1.e-6;

      
      if (cellprob){  //cell to debug   

        
        printf(" -------------------------------------------- \n");                
        printf(" ERROR LARGE  in CELL i=%d j=%d k=%d \n",i,j,k);
        if(std::abs(sumnorm)< 1.e-6) printf(" SUMNORMis very small = %e \n", sumnorm);
        printf("  X=%f Y=%f Z=%f  R=%f \n",xx,yy,zz,rad);
        printf( " marker (0) = %d marker(1) = %d \n",ebMarkers(i,j,k,0),ebMarkers(i,j,k,1));
        printf(" flag_arr(i,j,k).isCovered() = %d \n",flag_arr(i,j,k).isCovered());
        printf(" flag_arr(i,j,k).isSingleValued() = %d \n",flag_arr(i,j,k).isSingleValued());  

        printf(" lev = %d \n",lev);

        printf(" vol_centroid = %f %f \n",vol_centroid(i,j,k,0),vol_centroid(i,j,k,1));
        printf(" bc_centroid = %f %f \n",bc_centroid(i,j,k,0),bc_centroid(i,j,k,1));
        printf("  norm(%f ,%f , %f)     \n",norm_wall[0],norm_wall[1],norm_wall[2]);
        printf(" sumnorm = %e sumError=%e\n",sumnorm,sumError);
        
        
        printf(" dxinv = %f dyinv=%f dzinv=%f inv_hvfrac=%e\n",dxinv[0],dxinv[1],dxinv[2],inv_hvfrac);

        printf(" volfrac = %f \n",vfrac(i,j,k));
        printf(" area frac  ax(i)=%f ax(i+1)=%f \n",apx(i,j,k),apx(i+1,j,k));
        printf("            ay(j)=%f ay(j+1)=%f \n",apy(i,j,k),apy(i,j+1,k));
        printf("            az(k)=%f az(k+1)=%f \n",apz(i,j,k),apz(i,j,k+1));
        printf(" bc area = %f \n",bcarea(i,j,k,0));

        printf(" ============================================================\n");
        printf("  norm(%f ,%f , %f)     \n",norm_wall[0],norm_wall[1],norm_wall[2]);
        printf(">> OLD Error RESx=%e Resy=%e Resz=%e\n", Err_A[0],Err_A[1],Err_A[2]);      
        printf(">> OLD ERROR DIV  = %e \n", sumError*inv_hvfrac);


        Real nmag = 0.0_rt;
        for (int n = 0; n < AMREX_SPACEDIM; n++) { nmag += norm_wall[n]*norm_wall [n];}
        nmag = std::sqrt(nmag);
        //.. check error 
        Real Err_dir[AMREX_SPACEDIM]= {0.0};
        Err_dir[0] = apx(i,j,k) - apx(i+1,j,k) + areaw*norm_wall[0];
        Err_dir[1] = apy(i,j,k) - apy(i,j+1,k) + areaw*norm_wall[1];
        Err_dir[2] = apz(i,j,k) - apz(i,j,k+1) + areaw*norm_wall[2];
        printf(" ------------------------------------------------------------\n");
        printf(" norm wall = %f %f %f \n", norm_wall[0],norm_wall[1],norm_wall[2]);  
        printf(" nmmag norm wall = %f \n", nmag);      
        printf(" areaw = %f \n",areaw);
        printf(" Error CHECK (x,y,z) = (%e %e %e)\n", Err_dir[0], Err_dir[1], Err_dir[2]);      
        printf(">> ERROR DIV CHECK = %e \n", (Err_dir[0]+Err_dir[1]+Err_dir[2])*inv_hvfrac);


        // norm correction based on residual (snm new)  ***
        // normnew = norm_wall + ncorr
        // ncorr =  Err/area    
        Real ncorr[AMREX_SPACEDIM]= {0.0}, normnew[AMREX_SPACEDIM]= {0.0};
    
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
          ncorr[n]   = -Err_dir[n]/areaw;
          normnew[n] = norm_wall[n] + ncorr[n];
        }
        nmag = 0.0_rt;
        for (int n = 0; n < AMREX_SPACEDIM; n++) { nmag += normnew[n]*normnew[n];}
        nmag = std::sqrt(nmag);
                     
        // //.. check new error 
        // Err_dir[0] = apx(i,j,k) - apx(i+1,j,k) + areaw*normnew[0];
        // Err_dir[1] = apy(i,j,k) - apy(i,j+1,k) + areaw*normnew[1];
        // Err_dir[2] = apz(i,j,k) - apz(i,j,k+1) + areaw*normnew[2];
        // printf(" ------------------------------------------------------------\n");
        // printf("  (normal correction only  ** no normalised**)                \n");        
        // printf(" norm new = %f %f %f \n", normnew[0],normnew[1],normnew[2]);       
        // printf(" ncorr = %f %f %f \n", ncorr[0],ncorr[1],ncorr[2]);               
        // printf(" nmmag new = %f \n", nmag);
        // printf(" areaw = %f \n",areaw);
        // printf(" Error (x,y,z) = (%e %e %e)\n", Err_dir[0], Err_dir[1], Err_dir[2]);      
        // printf(">> ERROR DIV  = %e \n", (Err_dir[0]+Err_dir[1]+Err_dir[2])*inv_hvfrac);

        
        // Area correction based on residual (snm new)  ***

        Real area0 = bcarea(i,j,k,0);
        Err_dir[0] = apx(i,j,k) - apx(i+1,j,k) + area0*norm_wall[0];
        Err_dir[1] = apy(i,j,k) - apy(i,j+1,k) + area0*norm_wall[1];
        Err_dir[2] = apz(i,j,k) - apz(i,j,k+1) + area0*norm_wall[2];
      //  printf(" orig Error (x,y,z) = (%e %e %e)\n", Err_dir[0], Err_dir[1], Err_dir[2]);    
        Real Areai[AMREX_SPACEDIM]= {0.0}; 

        Areai[0] = area0*norm_wall[0] - Err_dir[0] ;
        Areai[1] = area0*norm_wall[1] - Err_dir[1] ;
        Areai[2] = area0*norm_wall[2] - Err_dir[2];

        Real areanew = Areai[0]*Areai[0] + Areai[1]*Areai[1] + Areai[2]*Areai[2];
        areanew = std::sqrt(areanew);
        Real normnew2[AMREX_SPACEDIM]= {0.0};
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
          normnew2[n] = Areai[n]/areanew;
        }
        nmag = 0.0_rt;
        for (int n = 0; n < AMREX_SPACEDIM; n++) { nmag += normnew2[n]*normnew2[n];}
        nmag = std::sqrt(nmag);

        //.. check new error 
        Err_dir[0] = apx(i,j,k) - apx(i+1,j,k) + areanew*normnew2[0];
        Err_dir[1] = apy(i,j,k) - apy(i,j+1,k) + areanew*normnew2[1];
        Err_dir[2] = apz(i,j,k) - apz(i,j,k+1) + areanew*normnew2[2];
        printf(" ------------------------------------------------------------\n");
        printf("  (NEW **)                \n");        
        printf(" norm new = %f %f %f \n", normnew2[0],normnew2[1],normnew2[2]);       
        printf(" nmmag new = %f \n", nmag);
        printf(" areaw = %f \n",areanew);
        printf(" Error NEW(x,y,z) = (%e %e %e)\n", Err_dir[0], Err_dir[1], Err_dir[2]);      
        printf(">> ERROR DIV NEW = %e \n", (Err_dir[0]+Err_dir[1]+Err_dir[2])*inv_hvfrac);




        // Real delta_2 = -( sumError)/(norm_wall[0]+norm_wall[1]+norm_wall[2] + 1.e-30);
        // Real area2 = bcarea(i,j,k,0) + delta_2;

        //.. check error 
        // Err_dir[0] = apx(i,j,k) - apx(i+1,j,k) + area2*norm_wall[0];
        // Err_dir[1] = apy(i,j,k) - apy(i,j+1,k) + area2*norm_wall[1];
        // Err_dir[2] = apz(i,j,k) - apz(i,j,k+1) + area2*norm_wall[2];
        // nmag = 0.0_rt;
        // for (int n = 0; n < AMREX_SPACEDIM; n++) { nmag += norm_wall[n]*norm_wall [n];}
        // nmag = std::sqrt(nmag);

        // printf(" ------------------------------------------------------------\n");
        // printf("  (area correction only)                                   \n");        
        // printf(" norm wall = %f %f %f \n", norm_wall[0],norm_wall[1],norm_wall[2]);  
        // printf(" nmag norm wall = %f \n", nmag);      
        // printf(" area2 = %f \n",area2);
        // printf(" Error (x,y,z) = (%e %e %e)\n", Err_dir[0], Err_dir[1], Err_dir[2]);      
        // printf(" >>> NEW ERROR DIV  (delta2) =%e \n", (Err_dir[0]+Err_dir[1]+Err_dir[2])*inv_hvfrac);

        printf(" ============================================================\n");

        printf(" i-1 i i+1  PRIMSWALL \n");
       for (int n = 0; n < 6; n++) {
        printf(" qwall(%d) q=%e \n",n,prim_wall[n]);         
       }
       
       printf(" i-1 i i+1  PRIMS \n");
       printf(" VFRAC %f %f %f \n",vfrac(i-1,j,k),vfrac(i,j,k),vfrac(i+1,j,k));
       printf(" MARKER0 %d %d %d \n",ebMarkers(i-1,j,k,0),ebMarkers(i,j,k,0),ebMarkers(i+1,j,k,0));
       printf(" MARKER1 %d %d %d \n",ebMarkers(i-1,j,k,1),ebMarkers(i,j,k,1),ebMarkers(i+1,j,k,1));
       
       
       for (int n = 0; n < 6; n++) {
        printf(" q(%d) : %e  %e  %e \n",n,prims(i-1,j,k,n),prims(i,j,k,n),prims(i+1,j,k,n));
       }

       printf(" j-1 j j+1  PRIMS \n");
       printf(" VFRAC %f %f %f \n",vfrac(i,j-1,k),vfrac(i,j,k),vfrac(i,j+1,k));
       printf(" MARKER0 %d %d %d \n",ebMarkers(i,j-1,k,0),ebMarkers(i,j,k,0),ebMarkers(i,j+1,k,0));
       printf(" MARKER1 %d %d %d \n",ebMarkers(i,j-1,k,1),ebMarkers(i,j,k,1),ebMarkers(i,j+1,k,1));
       
       for (int n = 0; n < 6; n++) {
        printf(" q(%d) : %e  %e  %e \n",n,prims(i,j-1,k,n),prims(i,j,k,n),prims(i,j+1,k,n));
       }

       printf(" k-1 k k+1  PRIMS \n");
       printf(" VFRAC %f %f %f \n",vfrac(i,j,k-1),vfrac(i,j,k),vfrac(i,j,k+1));
       printf(" MARKER0 %d %d %d \n",ebMarkers(i,j,k-1,0),ebMarkers(i,j,k,0),ebMarkers(i,j,k+1,0));
       printf(" MARKER1 %d %d %d \n",ebMarkers(i,j,k-1,1),ebMarkers(i,j,k,1),ebMarkers(i,j,k+1,1));
       
       for (int n = 0; n < 6; n++) {
        printf("  q(%d) : %e  %e  %e \n",n,prims(i,j,k-1,n),prims(i,j,k,n),prims(i,j,k+1,n));
       }
              
       printf(" FLUX WALL\n");
       for (int n = 0; n < 5; n++) {
        printf(" %d fluxw=%e \n",n,flux_wall[n]);  
       }
       printf(" FLUX X i i+1\n");
       for (int n = 0; n < 5; n++) {
        printf(" %d f(i)=%e f(i+1)=%f Sum(fi)=%e  Sum(fi)/dxcell=%e \n",n,flx_x(i,j,k,n),flx_x(i+1,j,k,n),  
            apx(i + 1, j, k) * flx_x(i+1,j,k,n) - apx(i, j, k) * flx_x(i,j,k,n), 
            (apx(i + 1, j, k) * flx_x(i+1,j,k,n) - apx(i, j, k) * flx_x(i,j,k,n))*inv_hvfrac);  
       }

       printf(" FLUX Y j j+1\n");
       for (int n = 0; n < 5; n++) {
        printf(" %d f(j)=%e f(j+1)=%e Sum(fj)=%e  Sum(fj)/dycell=%e \n",n,flx_y(i,j,k,n),flx_y(i,j+1,k,n), 
            apy(i, j, k + 1) * flx_y(i,j+1,k,n) - apy(i, j, k) * flx_y(i,j,k,n), 
            (apy(i, j, k + 1) * flx_y(i,j+1,k,n) - apy(i, j, k) * flx_y(i,j,k,n))*inv_hvfrac);  
       }

      
       printf(" FLUX Z k k+1\n");
       for (int n = 0; n < 5; n++) {
        printf(" %d f(k)=%e f(k+1)=%e Sum(fk)=%e  Sum(fk)/dzcell=%e\n",n,flx_z(i,j,k,n),flx_z(i,j,k+1,n),
            apz(i, j, k + 1) * flx_z(i,j,k+1,n) - apz(i, j, k) * flx_z(i,j,k,n),  
            (apz(i, j, k + 1) * flx_z(i,j,k+1,n) - apz(i, j, k) * flx_z(i,j,k,n))*inv_hvfrac);  
       }


       printf(" RHS  *******\n");
       for (int n = 0; n < 5; n++) {
          printf(" %d drhs +=%e   RHSTOT=%e \n",n,flux_wall[n]*inv_hvfrac*areaw,rhs(i,j,k,n)); 
        }



      } // end of cell to debug   
            
      //->  
