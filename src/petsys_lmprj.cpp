#include <petsys_lmprj.h>

//#define CFOV_ENABLED
#define DIST_NO_SQRT
#define BP_ATOMIC 1
const float DrawLine::PI = 3.1415926535898f;
const int DrawLine::GWSAMPLESIZE = 512;
const float DrawLine::FLT_MIN = 1e-30;

DrawLine::DrawLine(const std::size_t img_size_i,
                   const std::size_t img_size_j,
                   const std::size_t img_size_k,
                   const float vox_size_i,
                   const float vox_size_j,
                   const float vox_size_k,
                   const float tw_sigma,
                   const float tw_spacing) :
    m_img_size_i(img_size_i),
    m_img_size_j(img_size_j),
    m_img_size_k(img_size_k),
    m_vox_size_i(vox_size_i),
    m_vox_size_j(vox_size_j),
    m_vox_size_k(vox_size_k),
    m_tw_sigma(tw_sigma),
    m_tw_spacing(tw_spacing)
{
    // precompute image cuts
    for (std::size_t i = 0; i < img_size_i; i ++) {
        m_x_cuts.push_back((-(img_size_i * 0.5f) + i + 0.5f) * m_vox_size_i);
    }

    for (std::size_t j = 0; j < img_size_j; j ++) {
        m_y_cuts.push_back((-(img_size_j * 0.5f) + j + 0.5f) * m_vox_size_j);
    }

    for (std::size_t k = 0; k < img_size_k; k ++) {
        m_z_cuts.push_back((-(img_size_k * 0.5f) + k + 0.5f) * m_vox_size_k);
    }

    // image (negative) boundaries
    m_bd_min_x = -(vox_size_j * img_size_j * 0.5f);
    m_bd_min_y = -(vox_size_i * img_size_i * 0.5f);
    m_bd_min_z = -(vox_size_k * img_size_k * 0.5f);

#ifdef USE_TOF

    if (m_tw_sigma < 0.0f || m_tw_spacing < 0.0f) {
        SystemLog::write("invalid TOF parameters found: %f, %f", m_tw_sigma, m_tw_spacing);
        abort();
    }

#endif

    // precompute gaussian window lookuptable
    // +/- 3-sigma truncation
    SystemLog::write("precomputing timing window lookup table ...\n");
	float c0 = m_tw_spacing * 0.5;
	float c1 = m_tw_sigma*sqrt(2.0);
	float gw_stepsize = 3.0 * m_tw_sigma / GWSAMPLESIZE;
	m_gs_inv = 1.0 / gw_stepsize;
	float s0 = 0;
	for (int i = 0; i < GWSAMPLESIZE; i ++) {
		m_gw_lut.push_back((erf((i*gw_stepsize + c0)/c1) - erf((i*gw_stepsize - c0)/c1)) / 2.0);
	    s0 += m_gw_lut.back();
#if 0        
	    SystemLog::write("%f\n", m_gw_lut.back());
#endif        
	}

    // FOV radius
    if (m_img_size_i != m_img_size_j) {
        SystemLog::write("error, raytracer only support square image!");
        abort();
    }

    // use the radius of the largest circle that fits the square image
    // shrinkage a little bit to avoid boundary check error
    m_fov_radius = (m_img_size_i - 1) * m_vox_size_i * 0.5f;
}

DrawLine::~DrawLine()
{
}

inline bool DrawLine::hitCylFOV(float& tmin, float& tmax,
                                const float fov_radius,
                                const float p0x, const float dx,
                                const float p0y, const float dy)
{
    /*
    	cylindrical surface function
    		x^2 + y^2 = r^2
    	points on a ray
    		x = x0 + t*dx
    		y = y0 + t*dy
    	substitute them into the surface fucntion and
    	solve the quadratic equation to get t
     */
    float dx2dy2 = dx * dx + dy * dy;
    float det = dx2dy2 * fov_radius * fov_radius -
                (dx * p0y - dy * p0x) * (dx * p0y - dy * p0x);

    if (det > 0.0f) {
        float sqr_det = sqrtf(det);
        float xdx_ydy = -(dx * p0x + dy * p0y);
        tmin = (xdx_ydy - sqr_det) / dx2dy2;
        tmax = (xdx_ydy + sqr_det) / dx2dy2;
        return true;
    } else {
        return false;
    }
}



// v0: Jian's code for object's length longer than scanner
inline bool DrawLine::hitBoxFOV(float& tmin, float& tmax,
                                const float p0x, const float dx,
                                const float p0y, const float dy)
{
    //not use true boundaries but cuts, to avoid boundary check error
    float tx0 = (m_x_cuts[0] - p0x) / dx;
    float tx1 = (m_x_cuts[m_img_size_j - 1] - p0x) / dx;

    if (dx < 0.0f) {
        std::swap(tx0, tx1);
    }

    tx0 = std::max(tx0, 0.0f);
    tx1 = std::min(tx1, 1.0f);
    float ty0 = (m_y_cuts[0] - p0y) / dy;
    float ty1 = (m_y_cuts[m_img_size_i - 1] - p0y) / dy;

    if (dy < 0.0f) {
        std::swap(ty0, ty1);
    }

    ty0 = std::max(ty0, 0.0f);
    ty1 = std::min(ty1, 1.0f);

    tmin = std::max(tx0, ty0);
    tmax = std::min(tx1, ty1);

    return (tmin < tmax);
}




typedef int INT;
typedef float REAL;

//Xuezhu's modified code for EXPLORER (object's length maybe shorter than scanner)
// inline bool ImageRayTracer::hitCheck(const REAL p0x, const REAL p0y, const REAL p0z,
//                                      const REAL dx, const REAL dy, const REAL dz,
//                                      REAL& t_min, REAL& t_max)
inline bool DrawLine::hitCheck(const REAL p0x, const REAL p0y, const REAL p0z,
                               const REAL dx,  const REAL dy,  const REAL dz,
                               REAL& t_min, REAL& t_max)
{

    INT m_vdim_i = m_img_size_i;
    INT m_vdim_j = m_img_size_j;
    INT m_vdim_k = m_img_size_k;
    
    // REAL m_vox_size_i = m_vox_size_i;
    // REAL m_vox_size_j = m_vox_size_j;
    // REAL m_vox_size_k = m_vox_size_k;

    REAL m_vbd_x0;
    REAL m_vbd_x1;
    REAL m_vbd_y0;
    REAL m_vbd_y1;
    REAL m_vbd_z0;
    REAL m_vbd_z1;

    // m_vbd_y0 = -(vdim_i * vox_size_i) * 0.5;
    // m_vbd_y1 = +(vdim_i * vox_size_i) * 0.5;
    // m_vbd_x0 = -(vdim_j * vox_size_j) * 0.5;
    // m_vbd_x1 = +(vdim_j * vox_size_j) * 0.5;
    // m_vbd_z0 = -(vdim_k * vox_size_k) * 0.5;
    // m_vbd_z1 = +(vdim_k * vox_size_k) * 0.5;


    // m_bd_min_x = -(vox_size_j * img_size_j * 0.5f);
    // m_bd_min_y = -(vox_size_i * img_size_i * 0.5f);
    // m_bd_min_z = -(vox_size_k * img_size_k * 0.5f);


    // m_vbd_y0 =  m_bd_min_y;  // -
    // m_vbd_y1 = -m_bd_min_y;  // +
    // m_vbd_x0 =  m_bd_min_x;  // -
    // m_vbd_x1 = -m_bd_min_x;  // +
    // m_vbd_z0 =  m_bd_min_z;  // -
    // m_vbd_z1 = -m_bd_min_z;  // +


    m_vbd_x0 = m_x_cuts[0];
    m_vbd_x1 = m_x_cuts[m_img_size_j - 1];

    m_vbd_y0 = m_y_cuts[0];
    m_vbd_y1 = m_y_cuts[m_img_size_i - 1];

    m_vbd_z0 = m_z_cuts[0];
    m_vbd_z1 = m_z_cuts[m_img_size_k - 1];


#if CFOV_ENABLED


//    REAL fov_radius = 58 * m_vox_size_i * 0.5;
    REAL fov_radius = (m_img_size_i - 1) * m_vox_size_i * 0.5;
    REAL dx2dy2 = dx * dx + dy * dy;
    REAL det = dx2dy2 * fov_radius * fov_radius -
                 (dx * p0y - dy * p0x) * (dx * p0y - dy * p0x);

    if (det > 0.0) {
        REAL sqr_det = sqrtf(det);
        REAL xdx_ydy = -(dx * p0x + dy * p0y);
        t_min = (xdx_ydy - sqr_det) / dx2dy2;
        t_max = (xdx_ydy + sqr_det) / dx2dy2;
        return true;
    } else {
        return false;
    }



#else    
    // these are tricks to avoid the `divide by zero' error
    REAL ddx = (dx == 0) ? 1e-20 : dx;
    REAL ddy = (dy == 0) ? 1e-20 : dy;
    REAL ddz = (dz == 0) ? 1e-20 : dz;


    // parameter (need change if ray goes from p1 to p0)
#if 0
    REAL tx0 = (ddx > 0) ? ((m_vbd_x0 - m_vox_size_j - p0x) / dx) : ((m_vbd_x1 + m_vox_size_j - p0x) / ddx);
    REAL tx1 = (ddx > 0) ? ((m_vbd_x1 + m_vox_size_j - p0x) / dx) : ((m_vbd_x0 - m_vox_size_j - p0x) / ddx);
    REAL ty0 = (ddy > 0) ? ((m_vbd_y0 + m_vox_size_i - p0y) / dy) : ((m_vbd_y1 - m_vox_size_i - p0y) / ddy);
    REAL ty1 = (ddy > 0) ? ((m_vbd_y1 - m_vox_size_i - p0y) / dy) : ((m_vbd_y0 + m_vox_size_i - p0y) / ddy);
    REAL tz0 = (ddz > 0) ? ((m_vbd_z0 - p0z) / dz) : ((m_vbd_z1 - p0z) / ddz);
    REAL tz1 = (ddz > 0) ? ((m_vbd_z1 - p0z) / dz) : ((m_vbd_z0 - p0z) / ddz);
#else    
    REAL tx0 = (ddx > 0) ? ((m_vbd_x0 - p0x) / dx) : ((m_vbd_x1 - p0x) / ddx);
    REAL tx1 = (ddx > 0) ? ((m_vbd_x1 - p0x) / dx) : ((m_vbd_x0 - p0x) / ddx);
    REAL ty0 = (ddy > 0) ? ((m_vbd_y0 - p0y) / dy) : ((m_vbd_y1 - p0y) / ddy);
    REAL ty1 = (ddy > 0) ? ((m_vbd_y1 - p0y) / dy) : ((m_vbd_y0 - p0y) / ddy);
    REAL tz0 = (ddz > 0) ? ((m_vbd_z0 - p0z) / dz) : ((m_vbd_z1 - p0z) / ddz);
    REAL tz1 = (ddz > 0) ? ((m_vbd_z1 - p0z) / dz) : ((m_vbd_z0 - p0z) / ddz);
#endif
    /*
    if (dx == 0) {
        tx0 = -99.9;
        tx1 = +99.9;
    }

    if (dy == 0) {
        ty0 = -99.9;
        ty1 = +99.9;
    }

    if (dz == 0) {
        tz0 = -99.9;
        tz1 = +99.9;
    }     
     */
    

    // determine min and max
    // t_min = std::max(std::max(tx0, std::max(ty0, tz0)), 0.0); 
    // t_max = std::min(std::min(tx1, std::min(ty1, tz1)), 1.0);

    t_min = std::max(std::max(tx0, std::max(ty0, tz0)), 0.0f); 
    t_max = std::min(std::min(tx1, std::min(ty1, tz1)), 1.0f);   

    // float temp_min1, temp_min2, temp_max1, temp_max2;
    
    // temp_min1 = std::max(ty0, tz0);
    // temp_min2 = std::max(tx0, temp_min1);
    // t_min = std::max(temp_min2, 0.0); 

    // temp_max1 = std::min(ty1, tz1);
    // temp_max2 = std::min(tx1, temp_max1);
    // t_max = std::min(temp_max2, 1.0);

    return (t_min < t_max);


#endif
}






inline void DrawLine::calcTOFLOREndPoints(const int tbin_id,
        float& p0x, float& p0y, float& p0z,
        float& p1x, float& p1y, float& p1z,
        float& nrm_dx, float& nrm_dy, float& nrm_dz)
{
    // ray direction vector
    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;

    // lor center
    float pcx = (p0x + p1x) * 0.5f;
    float pcy = (p0y + p1y) * 0.5f;
    float pcz = (p0z + p1z) * 0.5f;

    // (inverse) length
    float inv_ray_length = 1.0f / sqrtf(dx * dx + dy * dy + dz * dz);

    // normalized dir
    nrm_dx = dx * inv_ray_length;
    nrm_dy = dy * inv_ray_length;
    nrm_dz = dz * inv_ray_length;

    float tbin_offset1 = tbin_id * m_tw_spacing - 3.0f * m_tw_sigma;
    float tbin_offset2 = tbin_id * m_tw_spacing + 3.0f * m_tw_sigma;

    // get new end-points
    p0x = pcx + tbin_offset1 * nrm_dx;
    p0y = pcy + tbin_offset1 * nrm_dy;
    p0z = pcz + tbin_offset1 * nrm_dz;
    p1x = pcx + tbin_offset2 * nrm_dx;
    p1y = pcy + tbin_offset2 * nrm_dy;
    p1z = pcz + tbin_offset2 * nrm_dz;
}

float DrawLine::fpBresenham(float p0x, float p0y, float p0z,
                            float p1x, float p1y, float p1z,
                            const Image<float>& img,
                            float weight, const int tbin_id)
{
    if (fabs(weight) < FLT_MIN) {
        return 0.0f;
    }

#ifdef USE_TOF
    // normalized dir
    float nrm_dx;
    float nrm_dy;
    float nrm_dz;
    calcTOFLOREndPoints(tbin_id, p0x, p0y, p0z, p1x, p1y, p1z, nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    float tbc_x = (p0x + p1x) * 0.5f;
    float tbc_y = (p0y + p1y) * 0.5f;
    float tbc_z = (p0z + p1z) * 0.5f;
#endif

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;
    float abs_dx = fabs(dx);
    float abs_dy = fabs(dy);
    float dt_x = m_vox_size_j / abs_dx;
    float dt_y = m_vox_size_i / abs_dy;
    float out = 0.0f;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
        return 0.0f;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return 0.0f;
    }

#endif

#ifdef USE_TOF

    tmin = std::max(tmin, 0.0f);
    tmax = std::min(tmax, 1.0f);

    // double check if the ray segment is still inside the FOV
    if (tmin > tmax) {
        return 0.0;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    float coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif

    if (abs_dx > abs_dy) { // x-cuts

        float d = abs_dy / abs_dx;
        weight *= sqrtf(1.0f + d * d);

        double x0 = p0x + tmin * dx;
        double x1 = p0x + tmax * dx;
        int j0 = int((x0 - m_bd_min_x) / m_vox_size_j); // init idx
        int j1 = int((x1 - m_bd_min_x) / m_vox_size_j); // end idx
        double t0 = (m_x_cuts[j0] - p0x) / dx; // init t
        double yc = -(m_bd_min_y + p0y + t0 * dy) / m_vox_size_i; // init coord y and z on x-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pyc = dt_x * dy / m_vox_size_i; // stepsize for coord y and z
        double pzc = dt_x * dz / m_vox_size_k;
        int j = j0; // voxel index
        int dj = (j0 > j1) ? -1 : 1;
        int jmin = std::min(j0, j1);
        int jmax = std::max(j0, j1);

        do {
            int yi = int(yc);
            int zi = int(zc);
            //
#ifdef USE_TOF // TOF
#ifndef DIST_NO_SQRT
            float tdx = m_x_cuts[j] - tbc_x;
            float tdy = m_y_cuts[m_img_size_i - yi - 1] - tbc_y;
            float tdz = m_z_cuts[zi] - tbc_z;
            float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
            float vx = m_x_cuts[j];
            float vy = m_y_cuts[m_img_size_i - yi - 1];
            float vz = m_z_cuts[zi];
            float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
            int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

            if (gw_bin < GWSAMPLESIZE) {
                out += img(yi, j, zi) * m_gw_lut[gw_bin];
            }

#else // non-TOF
            out += img(yi, j, zi);
#endif
            j += dj;
            jmin ++;
            yc -= pyc;
            zc += pzc;
        } while (jmin <= jmax);

    } else {

        float d = abs_dx / abs_dy;
        weight *= sqrtf(1.0f + d * d);

        double y0 = p0y + tmin * dy;
        double y1 = p0y + tmax * dy;
        int i0 = int((y0 - m_bd_min_y) / m_vox_size_i); // init idx
        int i1 = int((y1 - m_bd_min_y) / m_vox_size_i); // end idx
        double t0 = (m_y_cuts[i0] - p0y) / dy; // init t
        double xc = (p0x + t0 * dx - m_bd_min_x) / m_vox_size_j; // init coord x and z on y-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pxc = dt_y * dx / m_vox_size_j; // stepsize for coord y and z
        double pzc = dt_y * dz / m_vox_size_k;
        int i = m_img_size_i - i0 - 1; // voxel index (because of image coord)
        int di = (i0 > i1) ? 1 : -1;
        int imin = std::min(i0, i1);
        int imax = std::max(i0, i1);

        do {
            int xi = int(xc);
            int zi = int(zc);

#ifdef USE_TOF // TOF
#ifndef DIST_NO_SQRT
            float tdx = m_x_cuts[xi] - tbc_x;
            float tdy = m_y_cuts[m_img_size_i - i - 1] - tbc_y;
            float tdz = m_z_cuts[zi] - tbc_z;
            float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
            float vx = m_x_cuts[xi];
            float vy = m_y_cuts[m_img_size_i - i - 1];
            float vz = m_z_cuts[zi];
            float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
            int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

            if (gw_bin < GWSAMPLESIZE) {
                out += img(i, xi, zi) * m_gw_lut[gw_bin];
            }

#else // non-TOF
            out += img(i, xi, zi);
#endif

            i += di;
            imin ++;
            xc += pxc;
            zc += pzc;
        } while (imin <= imax);

    }

    return (out * weight);
}

void DrawLine::bpBresenham(float p0x, float p0y, float p0z,
                           float p1x, float p1y, float p1z,
                           Image<float>& img,
                           float weight, const int tbin_id)
{
    if (fabs(weight) < FLT_MIN) {
        return;
    }

#ifdef USE_TOF
    // normalized dir
    float nrm_dx;
    float nrm_dy;
    float nrm_dz;
    calcTOFLOREndPoints(tbin_id, p0x, p0y, p0z, p1x, p1y, p1z, nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    float tbc_x = (p0x + p1x) * 0.5f;
    float tbc_y = (p0y + p1y) * 0.5f;
    float tbc_z = (p0z + p1z) * 0.5f;
#endif

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;
    float abs_dx = fabs(dx);
    float abs_dy = fabs(dy);
    float dt_x = m_vox_size_j / abs_dx;
    float dt_y = m_vox_size_i / abs_dy;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
        return;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return;
    }

#endif

#ifdef USE_TOF

    tmin = std::max(tmin, 0.0f);
    tmax = std::min(tmax, 1.0f);

    // double check if the ray segment is still inside the FOV
    if (tmin > tmax) {
        return;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    float coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif

    /*
    	only check X-major or Y-major dir
    	usually no Z-major dir could happen
     */
    if (abs_dx > abs_dy) { // X-major

        float d = abs_dy / abs_dx;
        weight *= sqrtf(1.0f + d * d);

        double x0 = p0x + tmin * dx;
        double x1 = p0x + tmax * dx;
        int j0 = int((x0 - m_bd_min_x) / m_vox_size_j); // init idx
        int j1 = int((x1 - m_bd_min_x) / m_vox_size_j); // end idx
        double t0 = (m_x_cuts[j0] - p0x) / dx; // init t
        double yc = -(m_bd_min_y + p0y + t0 * dy) / m_vox_size_i; // init coord y and z on x-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pyc = dt_x * dy / m_vox_size_i; // stepsize for coord y and z
        double pzc = dt_x * dz / m_vox_size_k;
        int j = j0; // voxel index
        int dj = (j0 > j1) ? -1 : 1;
        int jmin = std::min(j0, j1);
        int jmax = std::max(j0, j1);

        do {
            int yi = int(yc);
            int zi = int(zc);

#ifdef USE_TOF // TOF
#ifndef DIST_NO_SQRT
            float tdx = m_x_cuts[j] - tbc_x;
            float tdy = m_y_cuts[m_img_size_i - yi - 1] - tbc_y;
            float tdz = m_z_cuts[zi] - tbc_z;
            float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
            float vx = m_x_cuts[j];
            float vy = m_y_cuts[m_img_size_i - yi - 1];
            float vz = m_z_cuts[zi];
            float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
            int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

            if (gw_bin < GWSAMPLESIZE) {
#ifdef USE_OMP //play trick
#if BP_ATOMIC
                #pragma omp atomic
#endif
#endif
                img(yi, j, zi) += weight * m_gw_lut[gw_bin];
                /* for boundary check (be sure never come here)
                if (j<0 || j>=m_img_size_j || yi<0 || yi>=m_img_size_i || zi<0 || zi>=m_img_size_k) {
                	std::printf("!!");
                	abort();
                }*/
            }

#else // non-TOF
//		if (yi >= 0 && yi < m_img_size_i &&
//			zi >= 0 && zi < m_img_size_k) { <<--- boundary check is not necessary!
#ifdef USE_OMP
#if BP_ATOMIC
            #pragma omp atomic
#endif
#endif
            img(yi, j, zi) += weight;
//		}
#endif

            j += dj;
            jmin ++;
            yc -= pyc;
            zc += pzc;
        } while (jmin <= jmax);

    } else { // Y-major

        float d = abs_dx / abs_dy;
        weight *= sqrtf(1.0f + d * d);

        double y0 = p0y + tmin * dy;
        double y1 = p0y + tmax * dy;
        int i0 = int((y0 - m_bd_min_y) / m_vox_size_i); // init idx
        int i1 = int((y1 - m_bd_min_y) / m_vox_size_i); // end idx
        double t0 = (m_y_cuts[i0] - p0y) / dy; // init t
        double xc = (p0x + t0 * dx - m_bd_min_x) / m_vox_size_j; // init coord x and z on y-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pxc = dt_y * dx / m_vox_size_j; // stepsize for coord y and z
        double pzc = dt_y * dz / m_vox_size_k;
        int i = m_img_size_i - i0 - 1; // voxel index (because of image coord)
        int di = (i0 > i1) ? 1 : -1;
        int imin = std::min(i0, i1);
        int imax = std::max(i0, i1);

        do {
            int xi = int(xc);
            int zi = int(zc);

#ifdef USE_TOF // TOF
#ifndef DIST_NO_SQRT
            float tdx = m_x_cuts[xi] - tbc_x;
            float tdy = m_y_cuts[m_img_size_i - i - 1] - tbc_y;
            float tdz = m_z_cuts[zi] - tbc_z;
            float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
            float vx = m_x_cuts[xi];
            float vy = m_y_cuts[m_img_size_i - i - 1];
            float vz = m_z_cuts[zi];
            float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
            int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

            if (gw_bin < GWSAMPLESIZE) {
#ifdef USE_OMP
#if BP_ATOMIC
                #pragma omp atomic
#endif
#endif
                img(i, xi, zi) += weight * m_gw_lut[gw_bin];
                /* for boundary check
                if (i<0 || i>=m_img_size_i || xi<0 || xi>=m_img_size_j || zi<0 || zi>=m_img_size_k) {
                	std::printf("!!");
                	abort();
                }*/
            }

#else // non-TOF
#ifdef USE_OMP
#if BP_ATOMIC
            #pragma omp atomic
#endif
#endif
            img(i, xi, zi) += weight;
#endif

            i += di;
            imin ++;
            xc += pxc;
            zc += pzc;
        } while (imin <= imax);
    }
}

int DrawLine::bpBresenham(float p0x, float p0y, float p0z,
                          float p1x, float p1y, float p1z,
                          Image<float>& img,
                          float weight, int* index_buff)
{
    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;
    float abs_dx = fabs(dx);
    float abs_dy = fabs(dy);
    float dt_x = m_vox_size_j / abs_dx;
    float dt_y = m_vox_size_i / abs_dy;
    int nnz = 0;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
        return 0;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return 0;
    }

#endif

    /*
    	only check X-major or Y-major dir
    	usually no Z-major dir could happen
     */
    if (abs_dx > abs_dy) { // X-major

        float d = abs_dy / abs_dx;
        weight *= sqrtf(1.0f + d * d);

        double x0 = p0x + tmin * dx;
        double x1 = p0x + tmax * dx;
        int j0 = int((x0 - m_bd_min_x) / m_vox_size_j); // init idx
        int j1 = int((x1 - m_bd_min_x) / m_vox_size_j); // end idx
        double t0 = (m_x_cuts[j0] - p0x) / dx; // init t
        double yc = -(m_bd_min_y + p0y + t0 * dy) / m_vox_size_i; // init coord y and z on x-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pyc = dt_x * dy / m_vox_size_i; // stepsize for coord y and z
        double pzc = dt_x * dz / m_vox_size_k;
        int j = j0; // voxel index
        int dj = (j0 > j1) ? -1 : 1;
        int jmin = std::min(j0, j1);
        int jmax = std::max(j0, j1);

        do {
            int yi = int(yc);
            int zi = int(zc);

//		if (yi >= 0 && yi < m_img_size_i &&
//			zi >= 0 && zi < m_img_size_k) { <<--- boundary check is not necessary!
            if (img(yi, j, zi) == 0) {
                img(yi, j, zi) = weight;
                index_buff[nnz] = img.sub2ind(yi, j, zi);
                nnz ++;
            } else {
                img(yi, j, zi) += weight;
            }

//		}
            j += dj;
            jmin ++;
            yc -= pyc;
            zc += pzc;
        } while (jmin <= jmax);

    } else { // Y-major

        float d = abs_dx / abs_dy;
        weight *= sqrtf(1.0f + d * d);

        double y0 = p0y + tmin * dy;
        double y1 = p0y + tmax * dy;
        int i0 = int((y0 - m_bd_min_y) / m_vox_size_i); // init idx
        int i1 = int((y1 - m_bd_min_y) / m_vox_size_i); // end idx
        double t0 = (m_y_cuts[i0] - p0y) / dy; // init t
        double xc = (p0x + t0 * dx - m_bd_min_x) / m_vox_size_j; // init coord x and z on y-cut
        double zc = (p0z + t0 * dz - m_bd_min_z) / m_vox_size_k;
        double pxc = dt_y * dx / m_vox_size_j; // stepsize for coord y and z
        double pzc = dt_y * dz / m_vox_size_k;
        int i = m_img_size_i - i0 - 1; // voxel index (because of image coord)
        int di = (i0 > i1) ? 1 : -1;
        int imin = std::min(i0, i1);
        int imax = std::max(i0, i1);

        do {
            int xi = int(xc);
            int zi = int(zc);

            if (img(i, xi, zi) == 0) {
                img(i, xi, zi) = weight; //tracer("%f %d\n", xc, i);
                index_buff[nnz] = img.sub2ind(i, xi, zi);
                nnz ++;
            } else {
                img(i, xi, zi) += weight;
            }

            i += di;
            imin ++;
            xc += pxc;
            zc += pzc;
        } while (imin <= imax);
    }

    return nnz;
}

float DrawLine::fpSiddon(float p0x, float p0y, float p0z,
                         float p1x, float p1y, float p1z,
                         const Image<float>& img,
                         float weight, const int tbin_id)
{
    if (fabs(weight) < FLT_MIN) {
        return 0.0f;
    }


    // printf("fpSiddon.debug.1: tbin_id=%d\n", tbin_id);

#ifdef USE_TOF
    // normalized dir
    float nrm_dx;
    float nrm_dy;
    float nrm_dz;
    calcTOFLOREndPoints(tbin_id, p0x, p0y, p0z, p1x, p1y, p1z, nrm_dx, nrm_dy, nrm_dz);


    // printf("fpSiddon.debug.2: tbin_id=%d\n", tbin_id);
      
    // timing bin center coord
    float tbc_x = (p0x + p1x) * 0.5f;
    float tbc_y = (p0y + p1y) * 0.5f;
    float tbc_z = (p0z + p1z) * 0.5f;
#endif

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;
    float out = 0.0f;



#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    // printf("!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)\n");

    // // v0: Jian's code for object's length longer than scanner
    // if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
    //     printf("!!!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)\n");
    //     return 0.0f;
    // }

    // v1: Xuezhu's modified code (EXPLORER) for object's length shorter than scanner
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, tmin, tmax)) {
        return 0;
    }


#else // for cylindrical FOV
    float tmin, tmax;

    // printf("!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)\n");


    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return 0.0f;
    }

#endif

#ifdef USE_TOF

    tmin = std::max(tmin, 0.0f);
    tmax = std::min(tmax, 1.0f);

    // double check if the ray segment is still inside the FOV
    if (tmin > tmax) {
        return 0.0;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    float coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif

    // printf("fpSiddon.debug.3: tmax=%f\n", tmax);


    // weight
    float w0 = sqrtf(dx * dx + dy * dy + dz * dz) * weight;

    // calculate the intersect point of ray with volume
    float p_1st_x = p0x + tmin * dx;
    float p_1st_y = p0y + tmin * dy;
    float p_1st_z = p0z + tmin * dz;

    // see the definition of coordination system
    // index of the first voxel
    int j = int((p_1st_x - m_bd_min_x) / m_vox_size_j); // 1: x
//  j = std::min(j, m_img_size_j - 1); // because of the maximum index value is m_vdim[*]-1
    int i = int((-m_bd_min_y - p_1st_y) / m_vox_size_i); // 0: y
//  i = std::min(i, m_img_size_i - 1);
    int k = int((p_1st_z - m_bd_min_z) / m_vox_size_k); // 2 : z
//  k = std::min(k, m_img_size_k - 1);

    // initial boundary
    float bx0 = (dx > 0.0f) ?
                m_bd_min_x + j * m_vox_size_j :
                m_bd_min_x + (j + 1) * m_vox_size_j;
    float by0 = (dy > 0.0f) ?
                m_bd_min_y + (m_img_size_i - i - 1) * m_vox_size_i :
                m_bd_min_y + (m_img_size_i - i) * m_vox_size_i;
    float bz0 = (dz > 0.0f) ?
                m_bd_min_z + k * m_vox_size_k :
                m_bd_min_z + (k + 1) * m_vox_size_k;

    int di = (dy > 0.0f) ? 1 : -1;
    int dj = (dx > 0.0f) ? 1 : -1;
    int dk = (dz > 0.0f) ? 1 : -1;
    float fx = (dx > 0.0f) ? m_vox_size_j : -m_vox_size_j;
    float fy = (dy > 0.0f) ? m_vox_size_i : -m_vox_size_i;
    float fz = (dz > 0.0f) ? m_vox_size_k : -m_vox_size_k;

    // local variables (using register for optimization)
    float tx = (dx == 0.0f) ? (99999999.9) : (bx0 + fx - p0x) / dx;
    float ty = (dy == 0.0f) ? (99999999.9) : (by0 + fy - p0y) / dy;
    float tz = (dz == 0.0f) ? (99999999.9) : (bz0 + fz - p0z) / dz;
    float vtx = (dx == 0.0f) ? (99999999.9) : fx / dx;
    float vty = (dy == 0.0f) ? (99999999.9) : fy / dy;
    float vtz = (dz == 0.0f) ? (99999999.9) : fz / dz;
    float tm0 = tmin, tm1;


    // printf("fpSiddon.debug.4: tbin_id=%d, tmax=%f, tm0=%f, tmax - tm0=%f\n", tbin_id, tmax, tm0, tmax - tm0);

    // ready to traval the image
    while (1) {

        // check if we are out of image
//            if (fabs(tm0 - tmax) < 1e-10) //(tm0 >= t_max) // added! 07-28-2011
//                break;
        if (tmax - tm0 < 1e-20)
            break;

        // boundaries of current voxel
        int i0 = i;
        int j0 = j;
        int k0 = k;

        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += vtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += vty;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        }

        // double check, because of numerical precision
        // the previous check might not be robust
        if (tm1 > tmax) {
            tm1 = tmax;
        }

        float ww = (tm1 - tm0) * w0;

        // put value to image domain
#ifdef USE_TOF
#ifndef DIST_NO_SQRT
        float tdx = m_x_cuts[j0] - tbc_x;
        float tdy = m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        float tdz = m_z_cuts[k0] - tbc_z;
        float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        float vx = m_x_cuts[j0];
        float vy = m_y_cuts[m_img_size_i - i0 - 1];
        float vz = m_z_cuts[k0];
        float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
        int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

        if (gw_bin < GWSAMPLESIZE) {
            out += img(i0, j0, k0) * ww * m_gw_lut[gw_bin];
        }

#else
        out += img(i0, j0, k0) * ww;

#endif

        tm0 = tm1; // save previous result // modified! 07-28-2011

    }


    // printf("out=%f, tbin_id=%d\n", out, tbin_id);


    // integral!
    return out;
}

void DrawLine::bpSiddon(float p0x, float p0y, float p0z,
                        float p1x, float p1y, float p1z,
                        Image<float>& img,
                        float weight, const int tbin_id)
{
    if (fabs(weight) < FLT_MIN) {
        return;
    }

#ifdef USE_TOF
    // normalized dir
    float nrm_dx;
    float nrm_dy;
    float nrm_dz;
    calcTOFLOREndPoints(tbin_id, p0x, p0y, p0z, p1x, p1y, p1z, nrm_dx, nrm_dy, nrm_dz);
    // timing bin center coord
    float tbc_x = (p0x + p1x) * 0.5f;
    float tbc_y = (p0y + p1y) * 0.5f;
    float tbc_z = (p0z + p1z) * 0.5f;
#endif

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    // // v0: Jian's code for object's length longer than scanner
    // if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
    //     return;
    // }

    // v1: Xuezhu's modified code (EXPLORER) for object's length shorter than scanner
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, tmin, tmax)) {
        return;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return;
    }

#endif

#ifdef USE_TOF

    tmin = std::max(tmin, 0.0f);
    tmax = std::min(tmax, 1.0f);

    // double check if the ray segment is still inside the FOV
    if (tmin > tmax) {
        return;
    }

#ifdef DIST_NO_SQRT // required by type-II distance (if enabled)
    float coef_d = -(nrm_dx * tbc_x + nrm_dy * tbc_y + nrm_dz * tbc_z);
#endif
#endif

    // weight
    float w0 = sqrtf(dx * dx + dy * dy + dz * dz) * weight;

    // calculate the intersect point of ray with volume
    float p_1st_x = p0x + tmin * dx;
    float p_1st_y = p0y + tmin * dy;
    float p_1st_z = p0z + tmin * dz;

    // see the definition of coordination system
    // index of the first voxel
    int j = int((p_1st_x - m_bd_min_x) / m_vox_size_j); // 1: x
//  j = std::min(j, m_img_size_j - 1); // because of the maximum index value is m_vdim[*]-1
    int i = int((-m_bd_min_y - p_1st_y) / m_vox_size_i); // 0: y
//  i = std::min(i, m_img_size_i - 1);
    int k = int((p_1st_z - m_bd_min_z) / m_vox_size_k); // 2 : z
//  k = std::min(k, m_img_size_k - 1);

    // initial boundary
    float bx0 = (dx > 0.0f) ?
                m_bd_min_x + j * m_vox_size_j :
                m_bd_min_x + (j + 1) * m_vox_size_j;
    float by0 = (dy > 0.0f) ?
                m_bd_min_y + (m_img_size_i - i - 1) * m_vox_size_i :
                m_bd_min_y + (m_img_size_i - i) * m_vox_size_i;
    float bz0 = (dz > 0.0f) ?
                m_bd_min_z + k * m_vox_size_k :
                m_bd_min_z + (k + 1) * m_vox_size_k;

    int di = (dy > 0.0f) ? 1 : -1;
    int dj = (dx > 0.0f) ? 1 : -1;
    int dk = (dz > 0.0f) ? 1 : -1;
    float fx = (dx > 0.0f) ? m_vox_size_j : -m_vox_size_j;
    float fy = (dy > 0.0f) ? m_vox_size_i : -m_vox_size_i;
    float fz = (dz > 0.0f) ? m_vox_size_k : -m_vox_size_k;

    // local variables (using register for optimization)
    float tx = (dx == 0.0f) ? 99999999.9 : (bx0 + fx - p0x) / dx;
    float ty = (dy == 0.0f) ? 99999999.9 : (by0 + fy - p0y) / dy;
    float tz = (dz == 0.0f) ? 99999999.9 : (bz0 + fz - p0z) / dz;
    float vtx = (dx == 0.0f) ? 99999999.9 : fx / dx;
    float vty = (dy == 0.0f) ? 99999999.9 : fy / dy;
    float vtz = (dz == 0.0f) ? 99999999.9 : fz / dz;
    float tm0 = tmin, tm1;

    // ready to traval the image
    while (1) {

        // check if we are out of image
//        if (fabs(tm0 - tmax) < 1e-10) //(tm0 >= t_max) // added! 07-28-2011
//            break;
        if (tmax - tm0 < 1e-20)
            break;

        int i0 = i;
        int j0 = j;
        int k0 = k;

        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += vtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += vty;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        }

        // double check, because of numerical precision
        // the previous check might not be robust
        if (tm1 > tmax) {
            tm1 = tmax;
        }

        // calc weight for current voxel
        float ww =  (tm1 - tm0) * w0;

        // put value to image domain
#ifdef USE_TOF
////<---
#ifndef DIST_NO_SQRT
        float tdx = m_x_cuts[j0] - tbc_x;
        float tdy = m_y_cuts[m_img_size_i - i0 - 1] - tbc_y;
        float tdz = m_z_cuts[k0] - tbc_z;
        float vox_to_tbc_dist = sqrtf(tdx * tdx + tdy * tdy + tdz * tdz);
#else
        float vx = m_x_cuts[j0];
        float vy = m_y_cuts[m_img_size_i - i0 - 1];
        float vz = m_z_cuts[k0];
        float vox_to_tbc_dist = fabs(nrm_dx * vx + nrm_dy * vy + nrm_dz * vz + coef_d);
#endif
        int gw_bin = int(vox_to_tbc_dist * m_gs_inv);

        if (gw_bin < GWSAMPLESIZE) {
#ifdef USE_OMP
#if BP_ATOMIC
            #pragma omp atomic
#endif
#endif
            img(i0, j0, k0) += ww * m_gw_lut[gw_bin];
        }

#else
////<---
#ifdef USE_OMP
#if BP_ATOMIC
        #pragma omp atomic
#endif
#endif
        img(i0, j0, k0) += ww;
#endif

        tm0 = tm1; // save previous result // modified! 07-28-2011
    };
}

int DrawLine::bpSiddon(float p0x, float p0y, float p0z,
                       float p1x, float p1y, float p1z,
                       Image<float>& img,
                       float weight, int* index_buff)
{
    if (fabs(weight) < FLT_MIN) {
        return 0;
    }

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;
    int nnz = 0;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    // // v0: Jian's code for object's length longer than scanner
    // if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
    //     return 0;
    // }

    // v1: Xuezhu's modified code (EXPLORER) for object's length shorter than scanner
    if (!hitCheck(p0x, p0y, p0z, dx, dy, dz, tmin, tmax)) {
        return 0;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return 0;
    }

#endif

    // weight
    float w0 = sqrtf(dx * dx + dy * dy + dz * dz) * weight;

    // calculate the intersect point of ray with volume
    float p_1st_x = p0x + tmin * dx;
    float p_1st_y = p0y + tmin * dy;
    float p_1st_z = p0z + tmin * dz;

    // see the definition of coordination system
    // index of the first voxel
    int j = int((p_1st_x - m_bd_min_x) / m_vox_size_j); // 1: x
//  j = std::min(j, m_img_size_j - 1); // because of the maximum index value is m_vdim[*]-1
    int i = int((-m_bd_min_y - p_1st_y) / m_vox_size_i); // 0: y
//  i = std::min(i, m_img_size_i - 1);
    int k = int((p_1st_z - m_bd_min_z) / m_vox_size_k); // 2 : z
//  k = std::min(k, m_img_size_k - 1);

    // initial boundary
    float bx0 = (dx > 0.0f) ?
                m_bd_min_x + j * m_vox_size_j :
                m_bd_min_x + (j + 1) * m_vox_size_j;
    float by0 = (dy > 0.0f) ?
                m_bd_min_y + (m_img_size_i - i - 1) * m_vox_size_i :
                m_bd_min_y + (m_img_size_i - i) * m_vox_size_i;
    float bz0 = (dz > 0.0f) ?
                m_bd_min_z + k * m_vox_size_k :
                m_bd_min_z + (k + 1) * m_vox_size_k;

    int di = (dy > 0.0f) ? 1 : -1;
    int dj = (dx > 0.0f) ? 1 : -1;
    int dk = (dz > 0.0f) ? 1 : -1;
    float fx = (dx > 0.0f) ? m_vox_size_j : -m_vox_size_j;
    float fy = (dy > 0.0f) ? m_vox_size_i : -m_vox_size_i;
    float fz = (dz > 0.0f) ? m_vox_size_k : -m_vox_size_k;

    // local variables (using register for optimization)
    float tx = (dx == 0.0f) ? 99999999.9 : (bx0 + fx - p0x) / dx;
    float ty = (dy == 0.0f) ? 99999999.9 : (by0 + fy - p0y) / dy;
    float tz = (dz == 0.0f) ? 99999999.9 : (bz0 + fz - p0z) / dz;
    float vtx = (dx == 0.0f) ? 99999999.9 : fx / dx;
    float vty = (dy == 0.0f) ? 99999999.9 : fy / dy;
    float vtz = (dz == 0.0f) ? 99999999.9 : fz / dz;
    float tm0 = tmin, tm1;

    // ready to traval the image
    while (1) {

        // check if we are out of image
//        if (fabs(tm0 - tmax) < 1e-10) //(tm0 >= t_max) // added! 07-28-2011
//            break;
        if (tmax - tm0 < 1e-20)
            break;

        int i0 = i;
        int j0 = j;
        int k0 = k;

        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += vtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += vty;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        }

        // double check, because of numerical precision
        // the previous check might not be robust
        if (tm1 > tmax) {
            tm1 = tmax;
        }

        // calc weight for current voxel
        float ww =  (tm1 - tm0) * w0;

        if (img(i0, j0, k0) == 0.0) {
            img(i0, j0, k0) = ww;
            index_buff[nnz] = img.sub2ind(i0, j0, k0);
            nnz ++;
        } else {
            img(i0, j0, k0) += ww;
        }

        tm0 = tm1; // save previous result // modified! 07-28-2011
    }

    return nnz;
}

void DrawLine::bpSiddonSquare(float p0x, float p0y, float p0z,
                        	  float p1x, float p1y, float p1z,
                        	  Image<float>& img,
                        	  float weight)
{
    if (fabs(weight) < FLT_MIN) {
        return;
    }

#ifdef USE_TOF
	SystemLog::write("bpSiddonSquare() cannot be used with TOF!\n");
	abort();
#endif

    float dx = p1x - p0x;
    float dy = p1y - p0y;
    float dz = p1z - p0z;

#ifndef CFOV_ENABLED // for square FOV
    float tmin, tmax;

    if (!hitBoxFOV(tmin, tmax, p0x, dx, p0y, dy)) {
        return;
    }

#else // for cylindrical FOV
    float tmin, tmax;

    if (!hitCylFOV(tmin, tmax, m_fov_radius, p0x, dx, p0y, dy)) {
        return;
    }

#endif

    // weight
    float w0 = sqrtf(dx * dx + dy * dy + dz * dz); // * weight;

    // calculate the intersect point of ray with volume
    float p_1st_x = p0x + tmin * dx;
    float p_1st_y = p0y + tmin * dy;
    float p_1st_z = p0z + tmin * dz;

    // see the definition of coordination system
    // index of the first voxel
    int j = int((p_1st_x - m_bd_min_x) / m_vox_size_j); // 1: x
//  j = std::min(j, m_img_size_j - 1); // because of the maximum index value is m_vdim[*]-1
    int i = int((-m_bd_min_y - p_1st_y) / m_vox_size_i); // 0: y
//  i = std::min(i, m_img_size_i - 1);
    int k = int((p_1st_z - m_bd_min_z) / m_vox_size_k); // 2 : z
//  k = std::min(k, m_img_size_k - 1);

    // initial boundary
    float bx0 = (dx > 0.0f) ?
                m_bd_min_x + j * m_vox_size_j :
                m_bd_min_x + (j + 1) * m_vox_size_j;
    float by0 = (dy > 0.0f) ?
                m_bd_min_y + (m_img_size_i - i - 1) * m_vox_size_i :
                m_bd_min_y + (m_img_size_i - i) * m_vox_size_i;
    float bz0 = (dz > 0.0f) ?
                m_bd_min_z + k * m_vox_size_k :
                m_bd_min_z + (k + 1) * m_vox_size_k;

    int di = (dy > 0.0f) ? 1 : -1;
    int dj = (dx > 0.0f) ? 1 : -1;
    int dk = (dz > 0.0f) ? 1 : -1;
    float fx = (dx > 0.0f) ? m_vox_size_j : -m_vox_size_j;
    float fy = (dy > 0.0f) ? m_vox_size_i : -m_vox_size_i;
    float fz = (dz > 0.0f) ? m_vox_size_k : -m_vox_size_k;

    // local variables (using register for optimization)
    float tx = (dx == 0.0f) ? 99999999.9 : (bx0 + fx - p0x) / dx;
    float ty = (dy == 0.0f) ? 99999999.9 : (by0 + fy - p0y) / dy;
    float tz = (dz == 0.0f) ? 99999999.9 : (bz0 + fz - p0z) / dz;
    float vtx = (dx == 0.0f) ? 99999999.9 : fx / dx;
    float vty = (dy == 0.0f) ? 99999999.9 : fy / dy;
    float vtz = (dz == 0.0f) ? 99999999.9 : fz / dz;
    float tm0 = tmin, tm1;

    // ready to traval the image
    while (1) {

        // check if we are out of image
//        if (fabs(tm0 - tmax) < 1e-10) //(tm0 >= t_max) // added! 07-28-2011
//            break;
        if (tmax - tm0 < 1e-20)
            break;

        int i0 = i;
        int j0 = j;
        int k0 = k;

        // check which direction should be updated
        if (tx < ty) {
            if (tx < tz) {
                tm1 = tx;
                j += dj;
                tx += vtx;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        } else {
            if (ty < tz) {
                tm1 = ty;
                i -= di;    // ! see definition of coordinate system
                ty += vty;
            } else {
                tm1 = tz;
                k += dk;
                tz += vtz;
            }
        }

        // double check, because of numerical precision
        // the previous check might not be robust
        if (tm1 > tmax) {
            tm1 = tmax;
        }

        // calc weight for current voxel
        float ww =  (tm1 - tm0) * w0;

        // put value to image domain
////<--- // nonTOF
#ifdef USE_OMP
#if BP_ATOMIC
        #pragma omp atomic
#endif
#endif
        img(i0, j0, k0) += ww * ww * weight;

        tm0 = tm1; // save previous result // modified! 07-28-2011
    };
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ListModeProjector
//
////////////////////////////////////////////////////////////////////////////////////////////////////

int ListModeProjector::raytracer_type = ListModeProjector::SIDDON; // defaul raytracer

ListModeProjector::ListModeProjector() : m_proj(0)
{
}

ListModeProjector::~ListModeProjector()
{
    if (!m_proj) {
        delete m_proj;
    }
}

bool ListModeProjector::initialize(const int num_of_rings,
                              	   const float ring_pitch,
                              	   const int num_of_blocks_t,
                              	   const int num_of_xtals_t,
                              	   const float ring_diameter,
                              	   const float xtal_pitch_t,
                              	   const float xtal_length,
                              	   const int img_size_ij,
                              	   const int img_size_k,
                              	   const float vox_size_ij,
                              	   const float vox_size_k,
                              	   const float tw_resolution,
                              	   const float tw_spacing,
                              	   const bool half_ang)
{
	SystemLog::write("initializing listmode-based projector ...\n");

	SystemLog::write("double checking projector parameters ...\n");
		
    if (num_of_rings <= 0 || ring_pitch <= 0.0 ||
        num_of_blocks_t <= 0 || num_of_xtals_t <= 0 ||
        ring_diameter <= 0.0 || img_size_ij <= 0 ||
        img_size_k <= 0 || vox_size_ij <= 0.0 || vox_size_k <= 0.0) {

        SystemLog::write("error, invalid parameter found: "
        				 "%d %f %d %d %f %d %d %f %f",
                    	 num_of_rings, ring_pitch, 
                    	 num_of_blocks_t, 
                    	 num_of_xtals_t,
                    	 ring_diameter, 
                    	 img_size_ij, img_size_k, 
                    	 vox_size_ij, vox_size_k);

        return false;
    }

//    double xtal_pitch_t = xtal_size_t + gap_size_t;
	SystemLog::write("calculating crystal transaxial positions (centers) ...\n");
	
    double ang_step = 360.0 / double(num_of_blocks_t);
    double ang_init = (half_ang) ? (0.5 * ang_step) : 0.0;
    double bx, by, angle;
    CRYSTALXY xtalxy;

    for (int a = 0; a < num_of_blocks_t; a ++) {
        angle = ang_init + ang_step * a;
        bx = (ring_diameter + xtal_length) * 0.5 * cos(angle * DrawLine::PI / 180.0);
        by = (ring_diameter + xtal_length) * 0.5 * sin(angle * DrawLine::PI / 180.0);
#if 0
		SystemLog::write("dectector module #%d\n", a+1);
#endif		
        for (int i = 0; i < num_of_xtals_t; i ++) {
            double x = xtal_pitch_t * (-num_of_xtals_t * 0.5 + 0.5 + i);
            xtalxy.x = bx + x * cos((angle + 90.0) * DrawLine::PI / 180.0);
            xtalxy.y = by + x * sin((angle + 90.0) * DrawLine::PI / 180.0);
            m_xtal_pos_xy.push_back(xtalxy);
#if 0            
            SystemLog::write("%.8f, %.8f\n", xtalxy.x, xtalxy.y);           
#endif             
        }
#if 0        
        SystemLog::write("\n");
#endif        
    }

    // calc ring axial pos
    SystemLog::write("calculating crystal ring offsets ...\n");
    for (int i = 0; i < num_of_rings; i ++) {
        m_ring_z.push_back((-num_of_rings * 0.5f + i + 0.5f)*ring_pitch);
#if 0        
        SystemLog::write("%.8f\n", m_ring_z.back());
#endif        
    }

    // create projector
    SystemLog::write("creating raytracer ... \n");
    m_proj = new DrawLine(img_size_ij, img_size_ij, img_size_k,
                          vox_size_ij, vox_size_ij, vox_size_k,
                          tw_resolution * 0.15 / 2.355,
                          tw_spacing * 0.15);

    return (m_proj != 0);
}

/*
void ListModeProjector::readRawData(const char* filename)
{
    std::ifstream input(filename);

    if (!input) {
#if DEBUG
        std::printf("open file failed: %s", filename);
#endif
        abort();
    } else {
        input.seekg(0, std::ios::end);
        size_t len_of_file = input.tellg();
        input.seekg(0, std::ios::beg);
        int num_of_events = len_of_file / sizeof(LMEVENT);

        if (num_of_events <= 0) {
#if DEBUG
            std::printf("no event found");
#endif
            abort();
        } else {

            m_raw_lmdata_buff.resize(num_of_events);
            input.read((char*)&m_raw_lmdata_buff[0], sizeof(LMEVENT) * num_of_events);

            if (input.fail()) {
#if DEBUG
                std::printf("read file failed: %s", filename);
#endif
                abort();
            }

        }
    }
}*/
void ListModeProjector::readRawData(const char* filename)
{
	SystemLog::write("reading raw data file ...\n");
	
	char str[512];
	SystemLog::write("checking event file ...\n");
	sprintf(str, "%s.lm", filename); 
	SystemLog::write("file in %s ...\n", str);
	SystemLog::write("calculating the number of events ...\n");
	size_t number_of_events = getFileLength(str) / sizeof(LMEVENT);
	SystemLog::write("number of events found: %lu\n", number_of_events);
	
	if (number_of_events == 0) {
		SystemLog::write("no event found, file may be invalid.\n");
		abort();
	}
	
	SystemLog::write("reading events ...\n");
	FILE* fid = fopen(str, "rb");
	if (fid != NULL) {
		m_raw_lmdata_buff.resize(number_of_events);
		fread((LMEVENT*)&m_raw_lmdata_buff[0], sizeof(LMEVENT), number_of_events, fid);
		fclose(fid);
		SystemLog::write("ok.\n");				
	} else {
		SystemLog::write("open file failed: %s.\n", str);
		abort();
	}	
	
	SystemLog::write("checking additive factor file ...\n");
	sprintf(str, "%s.add_fac", filename);
	SystemLog::write("file in %s ...\n", str);
	fid = fopen(str, "rb");
	if (fid != NULL) {
		m_add_fac_buff.resize(number_of_events);
		fread((LMEVENT*)&m_add_fac_buff[0], sizeof(float), number_of_events, fid);
		fclose(fid);
		SystemLog::write("ok.\n");
	} else {
		m_add_fac_buff.resize(number_of_events);
		SystemLog::write("open file failed: %s.\n", str);
		//abort();
		SystemLog::write("assuming default value ...\n");
		memset(&m_add_fac_buff[0], 0, sizeof(float) * number_of_events);
	}
	


    // assume mul_fac == 1.0
	SystemLog::write("checking multiplicative factor file ...\n");
	sprintf(str, "%s.mul_fac", filename);
	SystemLog::write("file in %s ...\n", str);
	fid = fopen(str, "rb");
	if (fid != NULL) {
		m_mul_fac_buff.resize(number_of_events);
		fread((LMEVENT*)&m_mul_fac_buff[0], sizeof(float), number_of_events, fid);
		fclose(fid);
		SystemLog::write("ok.\n");
	} else {
		SystemLog::write("open file failed: %s.\n", str);
		//abort();
		SystemLog::write("assuming default value ...\n");
		m_mul_fac_buff.resize(number_of_events);
		for (size_t n = 0; n < number_of_events; n ++) {
			m_mul_fac_buff[n] = 1.0;
		}
	}




/*
    SystemLog::write("checing angle index file ...\n");
    sprintf(str, "%s.aid", filename);
    SystemLog::write("file in %s ...\n", str);
    fid = fopen(str, "rb");
    if (fid != NULL) {
        m_angle_index_buff.resize(number_of_events);
        fread((int*)&m_angle_index_buff[0], sizeof(int), number_of_events, fid);
        fclose(fid);
        SystemLog::write("ok.\n");
    } else {
        SystemLog::write("open file failed: %s.\n", str);
        SystemLog::write("Warning: subsets will be made using uniform sampling ...\n");
    }
*/

#if 0 // random shuffle
	// added random shuffle
	SystemLog::write("performing random shuffle ...\n");
	SystemLog::write("generating random index ...\n");
	std::vector<int> random_index(number_of_events);
	for (int i = 0; i < number_of_events; i ++) {
		random_index[i] = i;
	}
	// create random index
	std::random_shuffle(random_index.begin(), random_index.end());

	//
	SystemLog::write("processing data ... ");
    std::vector<LMEVENT> tmp_lm(number_of_events);
    std::vector<float> tmp_add(number_of_events);
    std::vector<float> tmp_mul(number_of_events);

	// reorder listmode data
	for (int i = 0; i < number_of_events; i ++) {
        int j = random_index[i];
		tmp_lm[i] = m_raw_lmdata_buff[j];
		tmp_add[i] = m_add_fac_buff[j];
		tmp_mul[i] = m_mul_fac_buff[j];
	}

    memcpy(&m_raw_lmdata_buff[0], &tmp_lm[0], sizeof(LMEVENT) * number_of_events);
    memcpy(&m_add_fac_buff[0], &tmp_add[0], sizeof(float) * number_of_events);
    memcpy(&m_mul_fac_buff[0], &tmp_mul[0], sizeof(float) * number_of_events);

	SystemLog::write("done!\n");
#endif    
}

void ListModeProjector::makeSubsets(int number_of_subsets)
{
	SystemLog::write("making subsets [%d]...\n", number_of_subsets);
	
	m_number_of_subsets = number_of_subsets;
	
    if (m_raw_lmdata_buff.empty()) {
        std::printf("no data found, can't make subsets\n");
        abort();
    }
    
//    if (m_angle_index_buff.empty()) {

        // should match the subset in FProj and BProj
        SystemLog::write("resizing ... \n");
    //    m_lmdata_subsets.resize(number_of_subsets);
        m_mul_fac_subsets.resize(number_of_subsets);
        m_add_fac_subsets.resize(number_of_subsets);
        
        size_t number_of_events = m_raw_lmdata_buff.size();
    	size_t ne = number_of_events / number_of_subsets;
        for (size_t n = 0; n < number_of_subsets; n ++) {
    		if (n < number_of_events % number_of_subsets) {
    //		m_lmdata_subsets[n].resize(ne + 1);
    		    m_mul_fac_subsets[n].resize(ne + 1);
    		    m_add_fac_subsets[n].resize(ne + 1);
    		} else {
    //			m_lmdata_subsets[n].resize(ne);
    		    m_mul_fac_subsets[n].resize(ne);
    		    m_add_fac_subsets[n].resize(ne);
    		}
    		SystemLog::write("number of events in subset #%lu = %lu\n",
    			n + 1, m_mul_fac_subsets[n].size());
        }
        
        SystemLog::write("grouping ...\n");
        for (size_t n = 0; n < number_of_events; n ++) {
    //    	m_lmdata_subsets[n % number_of_subsets][n / number_of_subsets] = m_raw_lmdata_buff[n];
        	m_mul_fac_subsets[n % number_of_subsets][n / number_of_subsets] = m_mul_fac_buff[n];
        	m_add_fac_subsets[n % number_of_subsets][n / number_of_subsets] = m_add_fac_buff[n];
        }
        
        m_max_subset_size = m_mul_fac_subsets[0].size();
        
        SystemLog::write("ok.\n");

//    } else { // make geometric subsets

//        int number_of_angles = m_xtal_pos_xy.size() / 2; // compute total number of angles
//        if (number_of_angles % number_of_subsets != 0) {
//            SystemLog::write("invalid number of subsets, currently not supported yet!\n");
//           SystemLog::write("number of subsets must be an integer fraction of number of angles [%d]\n", number_of_angles);
//            abort();
//        }

//    }
}

int ListModeProjector::getNumberOfSubsets() const
{
    return m_number_of_subsets;
}

void ListModeProjector::doForwardProj(Image<float>& proj,
                                      const Image<float>& img,
                                      int subset_id) const
{

    // printf("Debug #1...doForwardProj\n");

	int ni = img.getDimI();
	int nj = img.getDimJ();
	int nk = img.getDimK();
	Image<float> tmp_img(ni, nj, nk);

	// r = R * x
	m_ipsf_model->blur(tmp_img, img, IPSFModel::F_BLURRING);
	
	// y = G * r
//    doForwardProj(proj, tmp_img, m_lmdata_subsets[subset_id]);
    doFProj(proj, tmp_img, subset_id);
}

void ListModeProjector::doBackProj(Image<float>& img,
                                   const Image<float>& proj,
                                   int subset_id) const
{
    int ni = img.getDimI();
	int nj = img.getDimJ();
	int nk = img.getDimK();
	Image<float> tmp_img(ni, nj, nk);
	
	// r = G' * y
//    doBackProj(tmp_img, proj, m_lmdata_subsets[subset_id]);
    doBProj(tmp_img, proj, subset_id);
    
    // x = R' * r;
    m_ipsf_model->blur(img, tmp_img, IPSFModel::B_BLURRING);
}

Image<float>* ListModeProjector::allocateProjection(int subset_id) const
{
	if (m_max_subset_size == 0) {
        std::printf("maximum subset size is equal to 0, can't allocate projection data\n");
        abort();
    }

    Image<float>* proj = new Image<float>(m_max_subset_size);
    return proj;
}

const std::vector<float>& ListModeProjector::getMultiplicativeFactor(int subset_id) const
{
    return m_mul_fac_subsets[subset_id];
}

const std::vector<float>& ListModeProjector::getAdditiveFactor(int subset_id) const
{
    return m_add_fac_subsets[subset_id];
}

void ListModeProjector::initializeSensitivityImage(const char* file, const SIZE& image_size)
{
	if (m_sensitivity == 0) {
		
		SystemLog::write("listmode projector: intializing sensitivity image ...\n");
		
		Image<float> sensitivity(image_size);
		if (!sensitivity.read(file)) {
			SystemLog::write("read file '%s' failed. use default sensitivity image.\n", file);
			sensitivity.set(1.0);
		}
		
		// create sensitivity
		m_sensitivity = new Image<float>(image_size);
		
		if (m_ipsf_model != 0) {
			m_ipsf_model->blur(*m_sensitivity, sensitivity, IPSFModel::B_BLURRING);

            SystemLog::write("\n\nm_ipsf_model->blur(*m_sensitivity, sensitivity, IPSFModel::B_BLURRING).\n\n");
            
		} else {
			SystemLog::write("error: iPSFModel not initialized! "
							 "required by sensitivity calculation.\n");
			abort();
		}
	}
}

const Image<float>* ListModeProjector::getSensitivity(int subset_id) const
{
	return m_sensitivity;
}

long long ListModeProjector::getFileLength(const char* filename)
{
    std::ifstream input(filename);

    if (!input) {
#if DEBUG
        std::printf("open file failed: %s", filename);
#endif
        return false;
    } else {
        input.seekg(0, std::ios::end);
        size_t len_of_file = input.tellg();
        input.seekg(0, std::ios::beg);
        input.close();

        return (long long)len_of_file;
    }
}

DrawLine* ListModeProjector::getRayTracer()
{
    return m_proj;
}

const std::vector<ListModeProjector::CRYSTALXY>& ListModeProjector::getXtalPosXY() const
{
    return m_xtal_pos_xy;
}

const std::vector<float>& ListModeProjector::getRingZ() const
{
    return m_ring_z;
}

void ListModeProjector::doForwardProj(Image<float>& proj,
                                      const Image<float>& img,
                                      const std::vector<LMEVENT>& lm) const
{
    
    // printf("Debug #2...doForwardProj\n");

    size_t num_of_events = lm.size();
    LMEVENT e;
    float p0x, p0y, p0z;
    float p1x, p1y, p1z;
    size_t i;

    switch (ListModeProjector::raytracer_type) {
    case SIDDON: {
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_FP)
#endif

        for (i = 0; i < num_of_events; i ++) { 
            e = lm[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];

            // printf("SIDDON: i=%lu, p0x=%f, p0y=%f, p0z=%f, p1x=%f, p1y=%f, p1z=%f, e.tbin_id=%d\n", i, p0x, p0y, p0z, p1x, p1y, p1z, e.tbin_id);

            proj[i] = m_proj->fpSiddon(p0x, p0y, p0z,
                                       p1x, p1y, p1z,
                                       img, 1.0f, e.tbin_id);

        }
    }
    break;

    case BRESENHAM: {
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_FP)
#endif

        for (i = 0; i < num_of_events; i ++) {
            e = lm[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];

            // printf("BRESENHAM: i=%lu, p0x=%f, p0y=%f, p0z=%f, p1x=%f, p1y=%f, p1z=%f, e.tbin_id=%d\n", i, p0x, p0y, p0z, p1x, p1y, p1z, e.tbin_id);

            proj[i] = m_proj->fpBresenham(p0x, p0y, p0z,
                                          p1x, p1y, p1z,
                                          img, 1.0f, e.tbin_id);

        }
    }
    break;

    default:
        std::printf("unknown raytracer\n");
        abort();
        break;
    }

}

void ListModeProjector::doBackProj(Image<float>& img,
                                   const Image<float>& proj,
                                   const std::vector<LMEVENT>& lm) const
{
    size_t num_of_events = lm.size();
    LMEVENT e;
    float p0x, p0y, p0z;
    float p1x, p1y, p1z;
    size_t i;

    switch (ListModeProjector::raytracer_type) {

    case SIDDON: {
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_BP)
#endif

        for (i = 0; i < num_of_events; i ++) {
            e = lm[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];
            m_proj->bpSiddon(p0x, p0y, p0z, p1x, p1y, p1z, img, proj[i], e.tbin_id);
        }
    }
    break;

    case BRESENHAM: {
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_BP)
#endif

        for (i = 0; i < num_of_events; i ++) {
            e = lm[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];
            m_proj->bpBresenham(p0x, p0y, p0z, p1x, p1y, p1z, img, proj[i], e.tbin_id);
        }
    }
    break;

    default:
        std::printf("unknown raytracer\n");
        abort();
        break;
    }
}

void ListModeProjector::doFProj(Image<float>& proj,
                                const Image<float>& img,
                                int subset_id) const
{
    

    size_t num_of_events = m_raw_lmdata_buff.size();

    // printf("Debug #3...doFProj: subset_id=%d, num_of_events=%lu\n", subset_id, num_of_events);

    LMEVENT e;
    float p0x, p0y, p0z;
    float p1x, p1y, p1z;
    size_t i;

    switch (ListModeProjector::raytracer_type) {
    case SIDDON: {

    // printf("switch-case SIDDON....\n");

        SystemLog::write("FP using SIDDON, initialized number of threads: %d (fp) \n", 
            Projector::NUMBER_OF_THREADS_FP
        );
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_FP)
#endif   

        for (i = subset_id; i < num_of_events; i += m_number_of_subsets) { 
            e = m_raw_lmdata_buff[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];   

            // printf("SIDDON: i=%lu, p0x=%f, p0y=%f, p0z=%f, p1x=%f, p1y=%f, p1z=%f, e.tbin_id=%d\n", i, p0x, p0y, p0z, p1x, p1y, p1z, e.tbin_id);

            proj[i / m_number_of_subsets] = m_proj->fpSiddon(p0x, p0y, p0z,
                                       						 p1x, p1y, p1z,
                                       						 img, 1.0f, e.tbin_id);

            // printf("proj[i / m_number_of_subsets]=%f \n", proj[i / m_number_of_subsets]);

        }
    }
    break;

    case BRESENHAM: {
        SystemLog::write("FP using BRESENHAM, initialized number of threads: %d (fp) \n", 
            Projector::NUMBER_OF_THREADS_FP
        );
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_FP)
#endif

        for (i = subset_id; i < num_of_events; i += m_number_of_subsets) {
            e = m_raw_lmdata_buff[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];
            
            // printf("BRESENHAM: i=%lu, p0x=%f, p0y=%f, p0z=%f, p1x=%f, p1y=%f, p1z=%f, e.tbin_id=%d\n", i, p0x, p0y, p0z, p1x, p1y, p1z, e.tbin_id);

            proj[i / m_number_of_subsets] = m_proj->fpBresenham(p0x, p0y, p0z,
                                          						p1x, p1y, p1z,
                                          						img, 1.0f, e.tbin_id);

        }
    }
    break;

    default:
        std::printf("unknown raytracer\n");
        abort();
        break;
    }

}

void ListModeProjector::doBProj(Image<float>& img,
                                const Image<float>& proj,
                                int subset_id) const
{
    size_t num_of_events = m_raw_lmdata_buff.size();
    LMEVENT e;
    float p0x, p0y, p0z;
    float p1x, p1y, p1z;
    size_t i;

    switch (ListModeProjector::raytracer_type) {

    case SIDDON: {
        SystemLog::write("BP using SIDDON, initialized number of threads: %d (bp) \n", 
            Projector::NUMBER_OF_THREADS_BP
        );
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_BP)
#endif

        for (i = subset_id; i < num_of_events; i += m_number_of_subsets) {
            e = m_raw_lmdata_buff[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];
            m_proj->bpSiddon(p0x, p0y, p0z, p1x, p1y, p1z, 
            				 img, proj[i / m_number_of_subsets], e.tbin_id);
        }
    }
    break;

    case BRESENHAM: {
        SystemLog::write("BP using BRESENHAM, initialized number of threads: %d (bp) \n", 
            Projector::NUMBER_OF_THREADS_BP
        );
#if USE_OMP
        #pragma omp parallel for private(i, e, p0x, p0y, p0z, p1x, p1y, p1z) num_threads(Projector::NUMBER_OF_THREADS_BP)
#endif

        for (i = subset_id; i < num_of_events; i += m_number_of_subsets) {
            e = m_raw_lmdata_buff[i];
            p0x = m_xtal_pos_xy[e.xtal_id_1].x;
            p0y = m_xtal_pos_xy[e.xtal_id_1].y;
            p0z = m_ring_z[e.ring_id_1];
            p1x = m_xtal_pos_xy[e.xtal_id_2].x;
            p1y = m_xtal_pos_xy[e.xtal_id_2].y;
            p1z = m_ring_z[e.ring_id_2];
            m_proj->bpBresenham(p0x, p0y, p0z, p1x, p1y, p1z, 
            					img, proj[i / m_number_of_subsets], e.tbin_id);
        }
    }
    break;

    default:
        std::printf("unknown raytracer\n");
        abort();
        break;
    }
}
