#include <petsys_recon.h>

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

    //****************************************************************
    //      get the number of subsets
    //****************************************************************
    SystemLog::write("getting number of subsets ... ");
    int number_of_subsets = projector->getNumberOfSubsets();
    SystemLog::write("%d ...\n", number_of_subsets);

    //****************************************************************
    //      initialization
    //****************************************************************
    int ni = x.getDimI();
    int nj = x.getDimJ();
    int nk = x.getDimK();
    SystemLog::write("allocating temporary images (size: %d,%d,%d) ...\n", ni, nj, nk);
    Image<float> xem(ni, nj, nk); // temporary image
    Image<float> xtemp1(ni, nj, nk); // temporary image
    Image<float> xtemp2(ni, nj, nk); // temporary image
    
    //****************************************************************
    //      initialize a vector for projection data
    //****************************************************************
    SystemLog::write("allocating temporary projections, ... ");
    int proj_buff_size = projector->allocateProjection(0)->getSize();
    SystemLog::write("length = %d\n", proj_buff_size);
    Image<float> y(proj_buff_size);

    //****************************************************************
    //      Start the reconstruction
    //****************************************************************
    SystemLog::write("ready, set, go ... (beta=%.2e) \n", beta);
    Timer timer;
    timer.start();

    for (int n = 0; n < number_of_iterations; n ++) {

        SystemLog::write("***** processing iteration #%d (%d in total) *****\n",
                         n + 1, number_of_iterations);

        for (int subset_id = 0; subset_id < number_of_subsets; subset_id ++) {

            SystemLog::write("processing subset #%d (%d in total) ...\n",
                             subset_id + 1, number_of_subsets);


            //****************************************************************
            //      Forward-project the current image x into y
            //****************************************************************
            // Likelihood part: EM surrogate
            //
            // P * x
            y.clear();
            projector->doForwardProj(y, x, subset_id); // projections in y


            //****************************************************************
            //      Calculate the ratio of measured sinogram to projected sinorgam (update term)
            //****************************************************************
            // calculate the ratio
            for (int i = 0; i < y.getSize(); i ++) {
                //
                float ytemp = y[i] + (projector->getAdditiveFactor(subset_id))[i];                          // add scatter estimation to the projection
                y[i] = (ytemp > 0.0) ? (projector->getMultiplicativeFactor(subset_id))[i] / ytemp : 0.0;    // divide measurement value by (projection+scatter)
            }

            //
            // if necessary check the objective function value
            //

            //****************************************************************
            //      Backproject the update term y into xem
            //****************************************************************
            // P' * y
            xem.clear();
            projector->doBackProj(xem, y, subset_id); // backprojected image in xtemp1

            //****************************************************************
            //      Multiply the current image (x) by the backprojected update term (xem)
            //****************************************************************
            // for update
            for (int j = 0; j < x.getSize(); j ++) {
                xem[j] *= x[j];
            }

            //****************************************************************
            //      Add regularizer terms
            //****************************************************************
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

                //****************************************************************
                //      Divide the image by the sensitivity map (if sensitivity is zero, set to zero)
                //****************************************************************
                // ML-EM update
                for (int j = 0; j < x.getSize(); j ++) {
                    x[j] = ((*projector->getSensitivity(subset_id))[j] > 0.0) ?
                                (xem[j] / ((*projector->getSensitivity(subset_id))[j] /
                                           double(number_of_subsets))) : 0.0;
                }
                
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

    }

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

