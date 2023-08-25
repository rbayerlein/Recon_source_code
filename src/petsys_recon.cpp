#include <petsys_recon.h>

using namespace std;

// Image<int>* gb_pnetknn(Image<float>& x)
// {
        
//     int psize = 3;  // patch size
//     int nsize = 7;  // neighbor search range
//     int nneighbor = 50; // Select nneighbor patches from nsize^3 neighbor patches
//     // By Guobao's default: psize = 3, nsize = 7, nneighbor = 50;
    
//     SystemLog::write("searching range = %d x %d x %d\n",nsize,nsize,nsize);
//     SystemLog::write("# of non-local patch = %d",nneighbor);
    
//     int tnsize = nsize*nsize*nsize;
//     int totalindex;
//     int ni = x.getDimI();
//     int nj = x.getDimJ();
//     int nk = x.getDimK();
    
//     int pt, p1, p2, p3;
    
//     int ns = floor(nsize/2);    //neighbor search range: -ns:ns
//     int ps = floor(psize/2);
    
//     Image<float> cpatch(psize,psize,psize); cpatch.clear(); // center patch
//     Image<float> npatch(psize,psize,psize); npatch.clear(); // neighbor patches
    
    
//     float tempdiff ,tempnp;
//     Image<float> dpatch(tnsize,2);  dpatch.clear();
//     // temp memory for neighbors of single voxel
//     // Top row, difference, bottom row, index.
    
    
//     int fdi;
//     Image<int>* fd = new Image<int>(ni*nj*nk,nneighbor,3);  
//     // Final result. contain search result of nearest voxel index for each voxel.
    
//     int i,j,k,i0,j0,k0;
//     // define i0,j0,k0 as neighbor voxels, also center of neighbor patches.

// //     for (k = 36; k < 37; k++)
// //         for (j = 150; j < 151; j++)
// //             for (i = 150; i < 151; i++)
    
    
// //     SystemLog::write("sb %d %d %d\n",ni,nj,nk);
        
//     for (k = 0; k < nk; k++){
//         for (j = 0; j < nj; j++){
//             for (i = 0; i < ni; i++){ // go through all image voxels
                
// //                 SystemLog::write("%d %d %d\n",i,j,k);
                
                
// //                 printf("%d %d %d\n",i,j,k);
//                 totalindex = k*ni*nj + j*ni + i;
// //                 printf("%d\n",totalindex);
//                 // center patch
//                 for(int c1=0;c1<psize;c1++){
//                     for(int c2=0;c2<psize;c2++){
//                         for(int c3=0;c3<psize;c3++){
                            
//                             i0 = i-ps+c1;
//                             j0 = j-ps+c2;
//                             k0 = k-ps+c3;
                            
//                             cpatch(c1,c2,c3)= ( ( (i0<0) || (i0>=ni) || 
//                                                   (j0<0) || (j0>=nj) ||  
//                                                   (k0<0) || (k0>=nk) )) ? 
//                                               0.0 : x(i0,j0,k0);               
                                              
//                         }
//                     }
//                 }
//                 // center patch obtained.
                
//                 // reset index to default
//                 dpatch.clear();
//                 for (int dp = 0; dp < tnsize; dp++)
//                 {
//                     dpatch(dp,1)=dp;
//                 }
                
                
//                 // neighbor patch
                 
//                  for(int c3=0;c3<nsize;c3++){
//                     for(int c2=0;c2<nsize;c2++){
//                         for(int c1=0;c1<nsize;c1++){
                            
                            
                            
//                             // now i0, j0 and k0 are center voxel of neighbor patches
//                             i0 = i-ns+c1;
//                             j0 = j-ns+c2;
//                             k0 = k-ns+c3;
                            
//                             int neighborindex = c1+c2*nsize+c3*nsize*nsize;
                                                        
//                             tempdiff = 0;
//                                 // traverse each voxel of neighbor patch
//                                 for(int d1=0;d1<psize;d1++){
//                                     for(int d2=0;d2<psize;d2++){
//                                         for(int d3=0;d3<psize;d3++){
                                        
//                                             int ii0 = i0-ps+d1;
//                                             int jj0 = j0-ps+d2;
//                                             int kk0 = k0-ps+d3;
//                                             tempnp = ( ( (ii0<0) || (ii0>=ni) ||
//                                                         (jj0<0) || (jj0>=nj) ||
//                                                         (kk0<0) || (kk0>=nk) )) ? 
//                                                     0.0 : x(ii0,jj0,kk0); 
                                                                                            
//                                             tempdiff += (tempnp-cpatch(d1,d2,d3))*(tempnp-cpatch(d1,d2,d3));
//                                 }}}
//                            dpatch(neighborindex,0) = tempdiff;   
//                            // Remember that dpatch(:,344/2) = 0 since it is the same as center patch
                
//                         }
//                     }
//                 }   // finish nsize^3 neighbor of one voxel
//                 // Now the first row of dpatch has the difference between center patch and neighbor patch
//                 // The second row is the index. Currently are 1:nsize^3
                

                 
//                 // Start bubble sort
//                 float tempdp;
//                 #pragma omp for
//                 for (int bb1 = 0; bb1 < tnsize; bb1++)
//                 {
//                     for (int bb2 = tnsize-1; bb2 > bb1; bb2--)
//                     {
//                         if ( dpatch(bb2,0) < dpatch(bb2-1,0) )
//                         {
//                             tempdp = dpatch(bb2,0);
//                             dpatch(bb2,0) = dpatch(bb2-1,0);
//                             dpatch(bb2-1,0) = tempdp;
                            
//                             tempdp = dpatch(bb2,1);
//                             dpatch(bb2,1) = dpatch(bb2-1,1);
//                             dpatch(bb2-1,1) = tempdp;
//                         }
//                     }
//                 }
                
//                 for (fdi = 1; fdi <= nneighbor; fdi++)
//                 // # of neighbor patch is nsize^3, which including center patch itself.
//                 // Hence, fd index start from 1. 
//                 {
//                     pt = int(dpatch(fdi,1));
//                     p3 = pt / (nsize*nsize);
//                     p2 = (pt - p3*nsize*nsize)/nsize;
//                     p1 = pt-p3*nsize*nsize-p2*nsize;
                    
//                     (*fd)(totalindex,fdi-1,0) = p1-ns;
//                     (*fd)(totalindex,fdi-1,1) = p2-ns;
//                     (*fd)(totalindex,fdi-1,2) = p3-ns;
//                 }

//     }}}// go through all image voxels
//     printf("temp OK ");  
    
//     (*fd).write("testfd");
    
//     printf("save OK");  
    
//     return fd;
// }


// void CIPEMBasedPLReconstructor::run(Projector* projector,
//                                  Regularizer* regularizer,
//                                  const int number_of_iterations,
//                                  const int stepsize_for_intermediate_result,
//                                  const float beta, 
//                                  const char* recon_output_folder,
//                                  const char* recon_output_filename_prefix,
//                                  Image<float>& x,
//                                  const float gb_alpha,
//                                  Image<float>& gb_cip)
// {
//     SystemLog::write("EM-based algorithm ...\n");

//     // get the number of subsets
//     SystemLog::write("getting number of subsets ... ");
//     int number_of_subsets = projector->getNumberOfSubsets();
//     SystemLog::write("%d ...\n", number_of_subsets);

//     // initialization
//     int ni = x.getDimI();
//     int nj = x.getDimJ();
//     int nk = x.getDimK();
//     SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
//     Image<float> xem(ni, nj, nk); // temporary image

//     // @Guobao
//     // temp image for guobao penalty, has to be outside if, otherwise compiler does not work...
//     Image<float> ztemp1(ni, nj, nk);    
//     Image<float> ztemp2(ni, nj, nk);
//     Image<int>* fd;
//     float gb_beta = 1 - gb_alpha;
    
//     if (gb_alpha > 0)
//     {
//         SystemLog::write("prepare CIP algorithm\n");
//         fd = gb_pnetknn(gb_cip);
//         (*fd).write("testfd");
        
//         SystemLog::write("Finish KNN for CIP\n");
// //         abort();
//     }


//     // temp image for ordinary penalty
//     Image<float> xtemp1(ni, nj, nk); // temporary image
//     Image<float> xtemp2(ni, nj, nk); // temporary image
    
//     // initialize a vector for projection data
//     SystemLog::write("allocating temporary projections, ... ");
//     size_t proj_buff_size = projector->allocateProjection(0)->getSize();
//     SystemLog::write("length = %lu\n", proj_buff_size);
//     Image<float> y(proj_buff_size);

//     SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
//     Timer timer;
//     timer.start();


//     char lk_fname[512];
//     sprintf(lk_fname, "%s/%s.likelihood",  // _intermediate_it%d.sub%d", 
//                       recon_output_folder,
//                       recon_output_filename_prefix);  //, n+1, subset_id); 
//     std::ofstream f_lk_log(lk_fname);
//     f_lk_log.precision(20);


//     // char lk2_fname[512];
//     // sprintf(lk2_fname, "%s/%s.likelihood_osem",  // _intermediate_it%d.sub%d", 
//     //                   recon_output_folder,
//     //                   recon_output_filename_prefix);  //, n+1, subset_id); 
//     // std::ofstream f2_lk_log(lk2_fname);
//     // f2_lk_log.precision(20);



//     // float cut_fov_radius = (float)(ceil(ni/2) + 10);
//     // float cut_fov_radius2 = cut_fov_radius * cut_fov_radius;

//     // printf("cut_fov_radius=%f, cut_fov_radius2=%f\n", cut_fov_radius, cut_fov_radius2);

//     // float center_fov_x = (ni + 1.0)/2.0;
//     // float center_fov_y = (nj + 1.0)/2.0;

//     // printf("center_fov_x=%f, center_fov_y=%f\n", center_fov_x, center_fov_y);



//     // int ix, iy, iz;



//     for (int n = 0; n < number_of_iterations; n ++) {   // 10000; n ++) {
    
//         SystemLog::write("***** processing iteration #%d (%d in total) *****\n", 
//             n + 1, number_of_iterations);
    
//             // double l1 = 0.0;
//             // for (int k = 0; k < x.getSize(); k ++) {
//             //     l1 += x[k] * (*projector->getSensitivity(subset_id))[k];
//             // }    

//         for (int subset_id = 0; subset_id < number_of_subsets; subset_id ++) {
        
//             SystemLog::write("processing subset #%d (%d in total) ...\n",
//                 subset_id + 1, number_of_subsets);


//             double l1 = 0.0;            
//             double l1_temp = 0.0;
//             for (int k = 0; k < x.getSize(); k ++) {
//                 l1_temp = x[k] * (*projector->getSensitivity(subset_id))[k];
//                 if (isnan(l1_temp)){
//                     l1_temp = 1e-20;
//                 } 
//                 else if (isinf(l1_temp)){
//                     l1_temp = 1e-20;
//                 }
//                 l1 += l1_temp;
//             }    

            
//             //
//             // Likelihood part: EM surrogate
//             //
//             // P * x
//             y.clear();

//             // printf("Debug...\n");

//             projector->doForwardProj(y, x, subset_id); // projections in y


//             // char str0[512];
//             // sprintf(str0, "%s/%s.x_intermediate_it%d.sub%d", 
//             //         recon_output_folder,
//             //         recon_output_filename_prefix, n+1, subset_id);
//             // x.write(str0);



//             // char str1[512];
//             // sprintf(str1, "%s/%s.y_intermediate_it%d.sub%d", 
//             //         recon_output_folder,
//             //         recon_output_filename_prefix, n+1, subset_id);
//             // y.write(str1);




//             // get difference
//             double l2 = 0.0;
//             double t = 0.0;  //, d = 0.0;

//             // calculate the ratio
//             for (size_t i = 0; i < y.getSize(); i ++) {
//                 //
//                 float ytemp = y[i] + (projector->getAdditiveFactor(subset_id))[i];
//                 y[i] = (ytemp > 0.0) ? (projector->getMultiplicativeFactor(subset_id))[i] / ytemp : 0.0;
//                 // y[i] = (ytemp > 0.0) ? 1.0 / ytemp : 0.0;
            
//                 // float yp = y[i] + 1e-20;  // bug
//                 float yp = ytemp + 1e-20;
//                 // // d += (1.0 / (projector->getMultiplicativeFactor(subset_id))[i]);
//                 // d += (1.0 / 1.0);

//                 if (isnan(yp)){
//                     yp = 1e-20;
//                 } 
//                 else if (isinf(yp)){
//                     yp = 1e-20;
//                 }

//                 l2 += (yp > 0.0) ? log(yp) : 0.0f;
// //              l2 += log(yp);
//                 t += yp;

// //                 printf("event.id=%d, ytemp=%f, y[i]=%f\n", i, ytemp, y[i]);
//             }




//             // char str2[512];
//             // sprintf(str2, "%s/%s.y1_intermediate_it%d.sub%d", 
//             //         recon_output_folder,
//             //         recon_output_filename_prefix, n+1, subset_id);
//             // y.write(str2);



//             // printf("lk = %.8f, l1 = %.8f, l2 = %.8f, [t=%f, d=%f]\n", -l1 + l2, l1, l2, t, d);
//             // f_lk_log << l1 << " " << l2 << " " << -l1 + l2 << std::endl;
//             printf("lk = %.8f, l1 = %.8f, l2 = %.8f, [t=%f]\n", -l1 + l2, l1, l2, t);
//             f_lk_log << l1 << " " << l2 << " " << -l1 + l2 << std::endl;



//             //
//             // if necessary check the objective function value
//             //

//             // P' * y
//             xem.clear();
//             projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1

//             // for update
//             for (int j = 0; j < x.getSize(); j ++) {
//                 xem[j] *= x[j];                
//                 // printf("x[j]=%f\n", x[j]);
//                 // printf("xem[j]=%f\n", xem[j]);
//                 // printf("%d, xem[j]=%f, x[j]=%f\n", j, xem[j], x[j]=%f);
//             }


//             // for (iz=0; iz<nk; iz++){
//             //     for (ix=0; ix<ni; ix++){
//             //         for (iy=0; iy<nj; iy++){
//             //             if( (ix - center_fov_x)*(ix - center_fov_x) + (iy - center_fov_y)*(iy - center_fov_y) > cut_fov_radius2){
//             //                 x(ix, iy, iz) = 0.0;
//             //                 xem(ix, iy, iz) = 0.0;
//             //             }
//             //         }
//             //     }
//             // }


//             // char str3[512];
//             // sprintf(str3, "%s/%s.xem_intermediate_it%d.sub%d", 
//             //         recon_output_folder,
//             //         recon_output_filename_prefix, n+1, subset_id);
//             // xem.write(str3);


//             // Regularization part: optimization transfer (half-quadratic approx + De Pierro's surrogate)
//             // The gradient of most surrogate functions has the general form below
//             //
//             // grad U(x; x^n) = diag(d) * (x - x^n) + D*x^n
//             //
//             // So we only need d and D*x^n
//             // define t1 as d, define t2 as D * x^n,
//             // the regularizer only needs to calculate t1 and t2 based on current image x^n
//             // Combining the EM surrogate, we then need to solve a quadratic function of the form
//             //
//             // -sensitivity * x + xem - beta*(diag(t1) * (x - x^n) + t2) * x = 0
//             // -beta*t1 * x * x + (beta * (t1*x^n - t2) - sensitivity) * x + xem = 0
//             // beta*t1 * x * x - (beta * (t1*x^n - t2) - sensitivity) * x - xem = 0
//             //
//             // i.e., a * x * x + b * x + c = 0
//             //
//             if (beta > 0) {

//                 // here t1 in xtemp1, t2 in xtemp2
//                 xtemp1.clear();
//                 xtemp2.clear();
//                 regularizer->calculateGradientBasedOnOptimizationTransfer(x, xtemp1, xtemp2);

//                 if (gb_alpha == 0)
//                 {
//                     // here t1 in xtemp1, t2 in xtemp2

//                     // then find the positive root of the quadratic function
//                     for (int j = 0; j < x.getSize(); j ++) {
//                         double a = beta * xtemp1[j];
//                         double b = (-beta * (xtemp1[j] * x[j] - xtemp2[j]) +
//                                   (*projector->getSensitivity(subset_id))[j]);
//                         double c = +xem[j] * double(number_of_subsets); // + 
                        
//                         //
//                         double det = sqrt(b * b + 4.0 * a * c);
//                         double d = ((b + det) > 0.0) ? (2.0*c) / (b + det) : 0.0;
//                         x[j] = (a == 0.0) ? (b != 0.0 ? c / b : 0.0) : float(d);
//                     }
//                 }
//                 else    // @ guobao
//                 {
//                     ztemp1.clear();
//                     ztemp2.clear();
//                     regularizer->gbcalg(x, ztemp1, ztemp2, gb_cip, fd);
                    
//                     for (int j = 0; j < x.getSize(); j ++) 
//                     {
//                         double a = beta * (gb_beta*xtemp1[j]+gb_alpha*ztemp1[j]);
//                         double b = (-beta * gb_beta  * (xtemp1[j] * x[j] - xtemp2[j])) +
//                                    (-beta * gb_alpha * (ztemp1[j] * x[j] - ztemp2[j])) +
//                                    (*projector->getSensitivity(subset_id))[j];
//                         double c = +xem[j] * double(number_of_subsets); // + 
                    
//                         double det = sqrt(b * b + 4.0 * a * c);
//                         double d = ((b + det) > 0.0) ? (2.0*c) / (b + det) : 0.0;
//                         x[j] = (a == 0.0) ? (b != 0.0 ? c / b : 0.0) : float(d);
//                     }                      
//                 }


//             } else {

//                 // ML-EM update
//                 for (int j = 0; j < x.getSize(); j ++) {
//                     x[j] = ((*projector->getSensitivity(subset_id))[j] > 0.0) ? 
//                            (xem[j] / ((*projector->getSensitivity(subset_id))[j] / 
//                                      double(number_of_subsets))) : 0.0;
//                 }
                
//             }

//             // for (iz=0; iz<nk; iz++){
//             //     for (ix=0; ix<ni; ix++){
//             //         for (iy=0; iy<nj; iy++){
//             //             if( (ix - center_fov_x)*(ix - center_fov_x) + (iy - center_fov_y)*(iy - center_fov_y) > cut_fov_radius2){
//             //                 x(ix, iy, iz) = 0.0;
//             //                 xem(ix, iy, iz) = 0.0;
//             //             }
//             //         }
//             //     }
//             // }

//         }

        

//         // save intermediate result if necessary
//         if ((((n+1) % stepsize_for_intermediate_result) == 0) && 
//             (stepsize_for_intermediate_result > 0) && 
//             (recon_output_folder != 0) &&
//             (recon_output_filename_prefix != 0)) {
//             char str[512];
//             sprintf(str, "%s/%s.intermediate.%d", 
//                     recon_output_folder,
//                     recon_output_filename_prefix, n+1);
//             x.write(str);
//         }

//     }

//     timer.stop();
// }





void EMBasedPLReconstructor::run(Projector* projector,
			   					 Regularizer* regularizer,
               					 const int number_of_iterations,
               					 const int stepsize_for_intermediate_result,
               					 const float beta, 
               					 const char* recon_output_folder,
               					 const char* recon_output_filename_prefix,
               					 Image<float>& x)
{
	SystemLog::write("EM-based algorithm ...\n");

    // get the number of subsets
    SystemLog::write("getting number of subsets ... ");
    int number_of_subsets = projector->getNumberOfSubsets();
    SystemLog::write("%d ...\n", number_of_subsets);

    // initialization
    int ni = x.getDimI();
    int nj = x.getDimJ();
    int nk = x.getDimK();
    SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
    Image<float> xem(ni, nj, nk); // temporary image
    Image<float> xtemp1(ni, nj, nk); // temporary image
    Image<float> xtemp2(ni, nj, nk); // temporary image
    
    // initialize a vector for projection data
    SystemLog::write("allocating temporary projections, ... ");
    size_t proj_buff_size = projector->allocateProjection(0)->getSize();
    SystemLog::write("length = %lu\n", proj_buff_size);
    Image<float> y(proj_buff_size);		// Px = forward projected image data
    Image<float> Px_r_s(proj_buff_size);	// fp image data plus randoms plus scatter (basically plus add_fac)

	SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
	Timer timer;
	timer.start();

    char lk_fname[512];
    sprintf(lk_fname, "%s/%s.likelihood",  // _intermediate_it%d.sub%d", 
                      recon_output_folder,
                      recon_output_filename_prefix);  //, n+1, subset_id); 
    std::ofstream f_lk_log(lk_fname);
    f_lk_log.precision(20);


    // char lk2_fname[512];
    // sprintf(lk2_fname, "%s/%s.likelihood_osem",  // _intermediate_it%d.sub%d", 
    //                   recon_output_folder,
    //                   recon_output_filename_prefix);  //, n+1, subset_id); 
    // std::ofstream f2_lk_log(lk2_fname);
    // f2_lk_log.precision(20);



    // float cut_fov_radius = (float)(ceil(ni/2) + 10);
    // float cut_fov_radius2 = cut_fov_radius * cut_fov_radius;

    // printf("cut_fov_radius=%f, cut_fov_radius2=%f\n", cut_fov_radius, cut_fov_radius2);

    // float center_fov_x = (ni + 1.0)/2.0;
    // float center_fov_y = (nj + 1.0)/2.0;

    // printf("center_fov_x=%f, center_fov_y=%f\n", center_fov_x, center_fov_y);



    // int ix, iy, iz;



    for (int n = 0; n < number_of_iterations; n ++) {   // 10000; n ++) {
    
    	SystemLog::write("***** processing iteration #%d (%d in total) *****\n", 
    		n + 1, number_of_iterations);
    


            // double l1 = 0.0;
            // for (int k = 0; k < x.getSize(); k ++) {
            //     l1 += x[k] * (*projector->getSensitivity(subset_id))[k];
            // }    


        for (int subset_id = 0; subset_id < number_of_subsets; subset_id ++) {
        
        	SystemLog::write("processing subset #%d (%d in total) ...\n",
        		subset_id + 1, number_of_subsets);


            double l1 = 0.0;            
            double l1_temp = 0.0;
            for (int k = 0; k < x.getSize(); k ++) {
                l1_temp = x[k] * (*projector->getSensitivity(subset_id))[k];
                if (isnan(l1_temp)){
                    l1_temp = 1e-20;
                } 
                else if (isinf(l1_temp)){
                    l1_temp = 1e-20;
                }
                l1 += l1_temp;
            }    

        	
            //
            // Likelihood part: EM surrogate
            //
            // P * x
			y.clear();
			Px_r_s.clear();

            // printf("Debug...\n");

            projector->doForwardProj(y, x, subset_id); // projections in y


            // char str0[512];
            // sprintf(str0, "%s/%s.x_intermediate_it%d.sub%d", 
            //         recon_output_folder,
            //         recon_output_filename_prefix, n+1, subset_id);
            // x.write(str0);

        	if((recon_output_folder != 0) &&// (n+1 == number_of_iterations) &&
        	(recon_output_filename_prefix != 0)) {
	            char str1[512];
	            sprintf(str1, "%s/%s.Px_it%d.sub%d", 
	                    recon_output_folder,
	                    recon_output_filename_prefix, n+1, subset_id);
	            y.write(str1);
	        }


            // get difference
            double l2 = 0.0;
            double t = 0.0;  //, d = 0.0;

            // calculate the ratio
            for (size_t i = 0; i < y.getSize(); i ++) {
            	//
            	float ytemp = y[i] + (projector->getAdditiveFactor(subset_id))[i];
            	Px_r_s[i] = ytemp;
                y[i] = (ytemp > 0.0) ? (projector->getMultiplicativeFactor(subset_id))[i] / ytemp : 0.0;
                // y[i] = (ytemp > 0.0) ? 1.0 / ytemp : 0.0;
            
                // float yp = y[i] + 1e-20;  // bug
                float yp = ytemp + 1e-20;
                // // d += (1.0 / (projector->getMultiplicativeFactor(subset_id))[i]);
                // d += (1.0 / 1.0);

                if (isnan(yp)){
                    yp = 1e-20;
                }
                else if (isinf(yp)){
                    yp = 1e-20;
                }

                l2 += (yp > 0.0) ? log(yp) : 0.0f;
//              l2 += log(yp);
                t += yp;

//                 printf("event.id=%d, ytemp=%f, y[i]=%f\n", i, ytemp, y[i]);
            }

        	if((recon_output_folder != 0) &&// (n+1 == number_of_iterations) &&
        	(recon_output_filename_prefix != 0)) {
	            char str1[512];
	            sprintf(str1, "%s/%s.Px_r_s_it%d.sub%d", 
	                    recon_output_folder,
	                    recon_output_filename_prefix, n+1, subset_id);
	            Px_r_s.write(str1);
	        }



            // char str2[512];
            // sprintf(str2, "%s/%s.y1_intermediate_it%d.sub%d", 
            //         recon_output_folder,
            //         recon_output_filename_prefix, n+1, subset_id);
            // y.write(str2);



            // printf("lk = %.8f, l1 = %.8f, l2 = %.8f, [t=%f, d=%f]\n", -l1 + l2, l1, l2, t, d);
            // f_lk_log << l1 << " " << l2 << " " << -l1 + l2 << std::endl;
            printf("lk = %.8f, l1 = %.8f, l2 = %.8f, [t=%f]\n", -l1 + l2, l1, l2, t);
            f_lk_log << l1 << " " << l2 << " " << -l1 + l2 << std::endl;



            //
            // if necessary check the objective function value
            //

            // P' * y
            xem.clear();
            projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1

            // for update
            for (int j = 0; j < x.getSize(); j ++) {
                xem[j] *= x[j];                
                // printf("x[j]=%f\n", x[j]);
                // printf("xem[j]=%f\n", xem[j]);
                // printf("%d, xem[j]=%f, x[j]=%f\n", j, xem[j], x[j]=%f);
            }


            // for (iz=0; iz<nk; iz++){
            //     for (ix=0; ix<ni; ix++){
            //         for (iy=0; iy<nj; iy++){
            //             if( (ix - center_fov_x)*(ix - center_fov_x) + (iy - center_fov_y)*(iy - center_fov_y) > cut_fov_radius2){
            //                 x(ix, iy, iz) = 0.0;
            //                 xem(ix, iy, iz) = 0.0;
            //             }
            //         }
            //     }
            // }


            // char str3[512];
            // sprintf(str3, "%s/%s.xem_intermediate_it%d.sub%d", 
            //         recon_output_folder,
            //         recon_output_filename_prefix, n+1, subset_id);
            // xem.write(str3);


            // Regularization part: optimization transfer (half-quadratic approx + De Pierro's surrogate)
            // The gradient of most surrogate functions has the general form below
            //
            // grad U(x; x^n) = diag(d) * (x - x^n) + D*x^n
            //
            // So we only need d and D*x^n
            // define t1 as d, define t2 as D * x^n,
            // the regularizer only needs to calculate t1 and t2 based on current image x^n
            // Combining the EM surrogate, we then need to solve a quadratic function of the form
            //
            // -sensitivity * x + xem - beta*(diag(t1) * (x - x^n) + t2) * x = 0
            // -beta*t1 * x * x + (beta * (t1*x^n - t2) - sensitivity) * x + xem = 0
            // beta*t1 * x * x - (beta * (t1*x^n - t2) - sensitivity) * x - xem = 0
            //
            // i.e., a * x * x + b * x + c = 0
            //
            if (beta > 0) {

                // here t1 in xtemp1, t2 in xtemp2
                xtemp1.clear();
                xtemp2.clear();
                regularizer->calculateGradientBasedOnOptimizationTransfer(x, xtemp1, xtemp2);

                // then find the positive root of the quadratic function
                for (int j = 0; j < x.getSize(); j ++) {
                    double a = beta * xtemp1[j];
                    double b = (-beta * (xtemp1[j] * x[j] - xtemp2[j]) +
                              (*projector->getSensitivity(subset_id))[j]);
                    double c = +xem[j] * double(number_of_subsets); // + 
                    
                    //
                    double det = sqrt(b * b + 4.0 * a * c);
                    double d = ((b + det) > 0.0) ? (2.0*c) / (b + det) : 0.0;
                    x[j] = (a == 0.0) ? (b != 0.0 ? c / b : 0.0) : float(d);
                }

            } else {

                // ML-EM update
                for (int j = 0; j < x.getSize(); j ++) {
                    x[j] = ((*projector->getSensitivity(subset_id))[j] > 0.0) ? 
                    	   (xem[j] / ((*projector->getSensitivity(subset_id))[j] / 
                    	   			 double(number_of_subsets))) : 0.0;
                }
                
            }

            // for (iz=0; iz<nk; iz++){
            //     for (ix=0; ix<ni; ix++){
            //         for (iy=0; iy<nj; iy++){
            //             if( (ix - center_fov_x)*(ix - center_fov_x) + (iy - center_fov_y)*(iy - center_fov_y) > cut_fov_radius2){
            //                 x(ix, iy, iz) = 0.0;
            //                 xem(ix, iy, iz) = 0.0;
            //             }
            //         }
            //     }
            // }

        }

        

        // save intermediate result if necessary
        if ((((n+1) % stepsize_for_intermediate_result) == 0) && 
        	(stepsize_for_intermediate_result > 0) && 
        	(recon_output_folder != 0) &&
        	(recon_output_filename_prefix != 0)) {
        	char str[512];
        	sprintf(str, "%s/%s.intermediate.%d", 
        			recon_output_folder,
        			recon_output_filename_prefix, n+1);
        	x.write(str);
        }

    }

	timer.stop();
}














void KEMBasedPLReconstructor::run(Projector* projector,
                                 Regularizer* regularizer,
                                 const int number_of_iterations,
                                 const int stepsize_for_intermediate_result,
                                 const float beta, 
                                 const char* recon_output_folder,
                                 const char* recon_output_filename_prefix,
                                 Image<float>& x)
{
    SystemLog::write("KEM-based algorithm ...\n");


    SystemLog::write("Loading Kernel Materix ... ");

    // const char* argv5 = argv[5];
    // char kmat_fname[512];
    // sprintf(kmat_fname, "%s", argv5);                     
    // std::stringstream ostr;
    // ostr.str("");      // clear
    // ostr << kmat_fname;

    std::stringstream ostr;
    ostr.str(""); // clear

    ostr << "./*.sysmat.239x239x679.sym.pmat";    

    printf("%s\n", ostr.str().c_str());


    PSpMatrix* kmat;    
    std::size_t total_sz = 0;  
    kmat = new PSpMatrix();
    kmat->read(ostr.str().c_str()); 
    total_sz += kmat->getNNZ()*sizeof(SSpMatrix::PACKAGE) +
                (kmat->getRowNum()+1)*sizeof(int);

    printf("OK! total memory space occupied: %lu (bytes)\n", total_sz);
    printf("kmat->getNNZ(): %lu (bytes)\n", kmat->getNNZ());
    printf("sizeof(SSpMatrix::PACKAGE): %lu (bytes)\n", sizeof(SSpMatrix::PACKAGE));
    printf("(kmat->getRowNum()+1): %lu (bytes)\n", (kmat->getRowNum()+1));
    printf("(kmat->getRowNum()+1)*sizeof(int): %lu (bytes)\n", (kmat->getRowNum()+1)*sizeof(int));






    // get the number of subsets
    SystemLog::write("getting number of subsets ... ");
    int number_of_subsets = projector->getNumberOfSubsets();
    SystemLog::write("%d ...\n", number_of_subsets);

    // initialization
    int ni = x.getDimI();
    int nj = x.getDimJ();
    int nk = x.getDimK();
    SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
    Image<float> xem(ni, nj, nk); // temporary image
    Image<float> xtemp1(ni, nj, nk); // temporary image
    Image<float> xtemp2(ni, nj, nk); // temporary image

    Image<float> alpha(ni, nj, nk); // temporary image
    Image<float> alpha_bk(ni, nj, nk); // temporary image

    alpha.set(1.0f);
    alpha_bk.set(1.0f);

    // for (int i = 0; i < x.getSize(); i ++) {
    //     x[i] = 1.0;
    //     alpha[i] = (*projector->getSensitivity(0)[i]>0.0) ? 1.0 : 0.0;
    //     alpha_bk[i] = (*projector->getSensitivity(0)[i]>0.0) ? 1.0 : 0.0;
    // }

    // initialize a vector for projection data
    SystemLog::write("allocating temporary projections, ... ");
    size_t proj_buff_size = projector->allocateProjection(0)->getSize();
    SystemLog::write("length = %lu\n", proj_buff_size);
    Image<float> y(proj_buff_size);

    SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
    Timer timer;
    timer.start();




    for (int n = 0; n < number_of_iterations; n ++) {   // 10000; n ++) {
    
        SystemLog::write("***** processing iteration #%d (%d in total) *****\n", 
            n + 1, number_of_iterations);
    
        for (int subset_id = 0; subset_id < number_of_subsets; subset_id ++) {
        
            SystemLog::write("processing subset #%d (%d in total) ...\n",
                subset_id + 1, number_of_subsets);



            y.clear();

            projector->doForwardProj(y, x, subset_id); // projections in y

            // calculate the ratio
            for (size_t i = 0; i < y.getSize(); i ++) {
                //
                float ytemp = y[i] + (projector->getAdditiveFactor(subset_id))[i];
                y[i] = (ytemp > 0.0) ? (projector->getMultiplicativeFactor(subset_id))[i] / ytemp : 0.0;
            }


            // P' * y
            xem.clear();
            projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1
            

            alpha_bk.clear();
            for (int m = 0; m < kmat->getRowNum(); m++) {
                float out = 0.0f;
                for (int n = kmat->getRowPtr()[m]; n < kmat->getRowPtr()[m+1]; n ++) {
                    int col = kmat->getCol()[n] - 1;
                    float val = kmat->getVal()[n];
                    out += xem[col] * val;
                }
                alpha_bk[m] = out;
            }


            // // for update
            // for (int j = 0; j < x.getSize(); j ++) {
            //     xem[j] *= x[j];                
            // }


            // ML-EM update
            for (int j = 0; j < x.getSize(); j ++) {
                alpha[j] = ((*projector->getSensitivity(subset_id))[j] > 0.0) ? 
                       (alpha[j] * alpha_bk[j] / ((*projector->getSensitivity(subset_id))[j] / 
                                 double(number_of_subsets))) : 0.0;
            }


            for (int m = 0; m < kmat->getRowNum(); m++) {
                float out = 0.0f;
                for (int n = kmat->getRowPtr()[m]; n < kmat->getRowPtr()[m+1]; n ++) {
                    int col = kmat->getCol()[n] - 1;
                    float val = kmat->getVal()[n];
                    out += alpha[col] * val;
                }
                x[m] = out;
            }


        }

        // save intermediate result if necessary
        if ((((n+1) % stepsize_for_intermediate_result) == 0) && 
            (stepsize_for_intermediate_result > 0) && 
            (recon_output_folder != 0) &&
            (recon_output_filename_prefix != 0)) {
            char str[512];
            sprintf(str, "%s/%s.intermediate.%d", 
                    recon_output_folder,
                    recon_output_filename_prefix, n+1);
            x.write(str);
        }


        // save intermediate result if necessary
        if ((((n+1) % stepsize_for_intermediate_result) == 0) && 
            (stepsize_for_intermediate_result > 0) && 
            (recon_output_folder != 0) &&
            (recon_output_filename_prefix != 0)) {
            char str[512];
            sprintf(str, "%s/%s.intermediate.%d_alpha", 
                    recon_output_folder,
                    recon_output_filename_prefix, n+1);
            alpha.write(str);
        }

    }

    timer.stop();
}















void ExpLineForwardProjector::run(Projector* projector,
                                 Regularizer* regularizer,
                                 const int number_of_iterations,
                                 const int stepsize_for_intermediate_result,
                                 const float beta, 
                                 const char* recon_output_folder,
                                 const char* recon_output_filename_prefix,
                                 Image<float>& x)
{
    SystemLog::write("Exp(-LineForwardProjector) algorithm ...\n");

    // get the number of subsets
    SystemLog::write("getting number of subsets ... ");
    int number_of_subsets = projector->getNumberOfSubsets();
    SystemLog::write("%d ...\n", number_of_subsets);

    // initialization
    int ni = x.getDimI();
    int nj = x.getDimJ();
    int nk = x.getDimK();
    SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
    Image<float> xem(ni, nj, nk); // temporary image
    Image<float> xtemp1(ni, nj, nk); // temporary image
    Image<float> xtemp2(ni, nj, nk); // temporary image
    
    // initialize a vector for projection data
    SystemLog::write("allocating temporary projections, ... ");
    size_t proj_buff_size = projector->allocateProjection(0)->getSize();
    SystemLog::write("length = %lu\n", proj_buff_size);
    Image<float> y(proj_buff_size);
    Image<float> y_exp(proj_buff_size);

    SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
    Timer timer;
    timer.start();


    // P * x
    y.clear();

    // printf("Debug...\n");
    int subset_id = 0; 
    projector->doForwardProj(y, x, subset_id); // projections in y

    for (size_t i = 0; i < y.getSize(); i ++) {
        y_exp[i] = exp(-y[i]);
    }


    // xem.clear();
    // projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1


    char str[512];
    sprintf(str, "%s/%s", 
            recon_output_folder,
            recon_output_filename_prefix);
    // y.write(str);
    y_exp.write(str);


    timer.stop();
}






void LineForwardProjector::run(Projector* projector,
                                 Regularizer* regularizer,
                                 const int number_of_iterations,
                                 const int stepsize_for_intermediate_result,
                                 const float beta, 
                                 const char* recon_output_folder,
                                 const char* recon_output_filename_prefix,
                                 Image<float>& x)
{
    SystemLog::write("LineForwardProjector algorithm ...\n");

    // get the number of subsets
    SystemLog::write("getting number of subsets ... ");
    int number_of_subsets = projector->getNumberOfSubsets();
    SystemLog::write("%d ...\n", number_of_subsets);

    // initialization
    int ni = x.getDimI();
    int nj = x.getDimJ();
    int nk = x.getDimK();
    SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
    Image<float> xem(ni, nj, nk); // temporary image
    Image<float> xtemp1(ni, nj, nk); // temporary image
    Image<float> xtemp2(ni, nj, nk); // temporary image
    
    // initialize a vector for projection data
    SystemLog::write("allocating temporary projections, ... ");
    size_t proj_buff_size = projector->allocateProjection(0)->getSize();
    SystemLog::write("length = %lu\n", proj_buff_size);
    Image<float> y(proj_buff_size);
    Image<float> y_exp(proj_buff_size);

    SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
    Timer timer;
    timer.start();


    // P * x
    y.clear();

    // printf("Debug...\n");
    int subset_id = 0; 
    projector->doForwardProj(y, x, subset_id); // projections in y

    // for (size_t i = 0; i < y.getSize(); i ++) {
    //     y_exp[i] = exp(-y[i]);
    // }


    // xem.clear();
    // projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1


    char str[512];
    sprintf(str, "%s/%s", 
            recon_output_folder,
            recon_output_filename_prefix);
    y.write(str);
    // y_exp.write(str);


    timer.stop();
}










///////////////////////////////////////////////////////////////////////////////////////////////////
float ADMMBasedPLReconstructor::mu = 1e6;

/*
	does not work very well, because of the choice of mu
*/
void ADMMBasedPLReconstructor::run(Projector* projector,
			   					   Regularizer* regularizer,
               					   const int number_of_iterations,
               					   const int stepsize_for_intermediate_result,
               					   const float beta, 
               					   const char* recon_output_folder,
               					   const char* recon_output_filename_prefix,
               					   Image<float>& x)
{
	int ni = x.getDimI();
	int nj = x.getDimJ();
	int nk = x.getDimK();
	
	Image<float> u(ni, nj, nk);
	Image<float> d(ni, nj, nk);
	Image<float> m(ni, nj, nk);
	
	Image<float> t1(ni, nj, nk);
	Image<float> t2(ni, nj, nk);

	int nite_sub = 10;	
	MeanBasedRegularizer* mean_reg = new MeanBasedRegularizer();
	
	memcpy(u.getPtr(), x.getPtr(), sizeof(float)*x.getSize());

	for (int n = 0; n < number_of_iterations; n ++) {

		// step 1: update u
		SystemLog::write("updating u ...\n");
		for (int j = 0; j < x.getSize(); j ++) {
			m[j] = x[j] - d[j];
		}
	
		for (int k = 0; k < nite_sub; k ++) {
			SystemLog::write("%d,", k);
			t1.clear();
		    t2.clear();
		    regularizer->calculateGradientBasedOnOptimizationTransfer(u, t1, t2);	
	
			for (int j = 0; j < x.getSize(); j ++) {
				u[j] = (mu + beta * t1[j] != 0) ? 
					(mu * m[j] + beta * t1[j] * u[j] - beta * t2[j]) / 
					(mu + beta * t1[j]) : 0.0;
			}	
		}
		SystemLog::write("\n");
	
		// step 2: update x
		SystemLog::write("updating x ...\n");
		for (int j = 0; j < x.getSize(); j ++) {
			m[j] = u[j] + d[j];
		}
		mean_reg->setMeanImage(m);
		EMBasedPLReconstructor::run(projector, mean_reg, 
									1, 0, mu * 0.5, 0, 0, x);
	
		// save intermediate result if necessary
        if ((((n+1) % stepsize_for_intermediate_result) == 0) && 
        	(stepsize_for_intermediate_result > 0) && 
        	(recon_output_folder != 0) &&
        	(recon_output_filename_prefix != 0)) {
        	char str[512];
        	sprintf(str, "%s/%s.intermediate.%d", 
        			recon_output_folder,
        			recon_output_filename_prefix, n+1);
        	x.write(str);
        }

		// step 3: update d
		SystemLog::write("updating d ...\n");
		for (int j = 0; j < x.getSize(); j ++) {
			d[j] = d[j] - (x[j] - u[j]);
		}
	}
}               					 

