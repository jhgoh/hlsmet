#ifndef MET_H
#define MET_H

#include <iostream>
#include <cmath>

// top algs
void met_ref(float in_pt[NPART], float in_phi[NPART], float& out_pt, float& out_phi);
void met_hw(pt_t data_pt[NPART], phi_t data_phi[NPART], pt2_t& res_pt2, phi_t& res_phi);

/*
#define PROJ_TAB_SIZE (1<<(PHI_SIZE-2))
template<class pt_T, class phi_T,class pxy_T>
void ProjXY(pt_T pt, phi_T phi, pxy_T& px, pxy_T& py, bool debug=false)
{
  // Initialize the lookup tables
#ifdef __HLS_SYN__
  bool initialized = false;
  pt_t cos_table[PROJ_TAB_SIZE];
#else 
  static bool initialized = false;
  static pt_t cos_table[PROJ_TAB_SIZE];
#endif
  if (!initialized) {
    init_projx_table(cos_table);
    initialized = true;
  }
    //map phi to first quadrant value: range [0, 2^(PHI_SIZE-2))
    ap_uint<PHI_SIZE-2> phiQ1 = phi;
    if(phi>=(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map 64-128 (0-63) to 63-0
    if(phi<0 && phi>=-(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map -64-1 (0-63) to 63-0

    // get x component and flip sign if necessary
    x = (pt * cos_table[phiQ1]) >> PT_SIZE;
    if(debug) std::cout << pt << "  cos_table[" << phiQ1 << "] = " << cos_table[phiQ1] << "  " << x << std::endl;
    if( phi>=(1<<(PHI_SIZE-2))
        || phi<-(1<<(PHI_SIZE-2)))
        x = -x;

    return;
}


//
// Lookup tables for pt projections to X, Y
//   n.b. store only (0,pi/2) to reduce size
//

#define PROJ_TAB_SIZE (1<<(PHI_SIZE-2))
// use table size of 1024/4=256 for now, so no precision is lost

template<class pt_T>
void init_projx_table(pt_T table_out[PROJ_TAB_SIZE]) {
    // Return table of cos(phi) where phi is in (0,pi/2)
    // multiply result by 2^(PT_SIZE) (equal to 1 in our units)
    for (int i = 0; i < PROJ_TAB_SIZE; i++) {
        //store result, guarding overflow near costheta=1
        pt2_t x = round((1<<PT_SIZE) * cos(float(i)/PROJ_TAB_SIZE * M_PI/2));
        // (using extra precision here (pt2_t, not pt_t) to check the out of bounds condition)
        if(x >= (1<<PT_SIZE)) table_out[i] = (1<<PT_SIZE)-1;
        else table_out[i] = x;
    }
    return;
}

template<class pt_T, class phi_T,class pxy_T>
    void ProjX(pt_T pt, phi_T phi, pxy_T &x, bool debug=false){
    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t cos_table[PROJ_TAB_SIZE];
#else 
    static bool initialized = false;
    static pt_t cos_table[PROJ_TAB_SIZE];
#endif
    if (!initialized) {
        init_projx_table(cos_table);
        initialized = true;
    }
    //map phi to first quadrant value: range [0, 2^(PHI_SIZE-2))
    ap_uint<PHI_SIZE-2> phiQ1 = phi;
    if(phi>=(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map 64-128 (0-63) to 63-0
    if(phi<0 && phi>=-(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map -64-1 (0-63) to 63-0

    // get x component and flip sign if necessary
    x = (pt * cos_table[phiQ1]) >> PT_SIZE;
    if(debug) std::cout << pt << "  cos_table[" << phiQ1 << "] = " << cos_table[phiQ1] << "  " << x << std::endl;
    if( phi>=(1<<(PHI_SIZE-2))
        || phi<-(1<<(PHI_SIZE-2)))
        x = -x;

    return;
}



template<class pt_T>
void init_projy_table(pt_T table_out[PROJ_TAB_SIZE]) {
    // Return table of sin(phi) where phi is in (0,pi/2)
    // multiply result by 1=2^(PT-SIZE)
    // see comments in the above ProjX function
    for (int i = 0; i < PROJ_TAB_SIZE; i++) {
        pt2_t x = round((1<<PT_SIZE) * sin(float(i)/PROJ_TAB_SIZE * M_PI/2));
        if(x >= (1<<PT_SIZE)) table_out[i] = (1<<PT_SIZE)-1;
        else table_out[i] = x;
    }
    return;
}

template<class pt_T, class phi_T,class pxy_T>
    void ProjY(pt_T pt, phi_T phi, pxy_T &y, bool debug=false){
    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t sin_table[PROJ_TAB_SIZE];
#else 
    static bool initialized = false;
    static pt_t sin_table[PROJ_TAB_SIZE];
#endif
    if (!initialized) {
        init_projy_table(sin_table);
        initialized = true;
    }
    //map phi to first quadrant value: range [0, 2^(PHI_SIZE-2))
    ap_uint<PHI_SIZE-2> phiQ1 = phi;
    if(phi>=(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map 64-128 (0-63) to 63-0
    if(phi<0 && phi>=-(1<<(PHI_SIZE-2))) phiQ1 = (1<<(PHI_SIZE-2)) -1 - phiQ1; // map -64-1 (0-63) to 63-0

    // get y component and flip sign if necessary
    y = (pt * sin_table[phiQ1]) >> PT_SIZE;
    if(debug) std::cout << pt << "  sin_table[" << phiQ1 << "] = " << sin_table[phiQ1] << "  " << y << std::endl;
    if( phi<0 ) y = -y;

    return;
}





//
// This may not the most efficient way, I think, since numbers 513-1023 all map to 1 !
//   Currently we mitigate, by simply not storing inverses of large numbers, 
//   but there may be a better way to do this (TODO - study this!)
//
#define DROP_BITS 2
#define INV_TAB_SIZE (1<<(PT_SIZE-DROP_BITS))
// return the inverse of a 'pt_size' bits (1024) number
template<class pt_T>
void init_inv_table(pt_T table_out[INV_TAB_SIZE]) {
    // multiply result by 1=2^(PT-SIZE)
    table_out[0]=(1<<PT_SIZE)-1;
    for (int i = 1; i < INV_TAB_SIZE; i++) {
        table_out[i] = round((1<<PT_SIZE) / float(i));
    }
    return;
}

// sets number of values in atan lookup table.
//   covers 0,pi/4, so naturally use 2^(PHI_SIZE)/8
#define ATAN_SIZE (PHI_SIZE-3)
#define ATAN_TAB_SIZE (1<<ATAN_SIZE)
// Get arctan of a number in (0,1), represented as integers 0 to 2^(pt_size)=1024
// for a 1-1 mapping, we can get away with an eighth of the bits used for full phi
template<class phi_T>
void init_atan_table(phi_T table_out[ATAN_TAB_SIZE]) {
    // multiply result by 1=2^(PT-SIZE) 
    table_out[0]=int(0);
    for (int i = 1; i < ATAN_TAB_SIZE; i++) {
        table_out[i] = int(round(atan(float(i)/ATAN_TAB_SIZE) * (1<<(PHI_SIZE-3)) / (M_PI/4)));
    }
    return;
}

template<class pxy_T, class phi_T>
    void PhiFromXY(pxy_T px, pxy_T py, phi_T &phi){

    // Initialize the lookup tables
#ifdef __HLS_SYN__
    bool initialized = false;
    pt_t inv_table[INV_TAB_SIZE];
    pt_t atan_table[ATAN_TAB_SIZE];
#else 
    static bool initialized = false;
    static pt_t inv_table[INV_TAB_SIZE];
    static pt_t atan_table[ATAN_TAB_SIZE];
#endif
    if (!initialized) {
        init_inv_table(inv_table);
        init_atan_table(atan_table);
        initialized = true;
    }

    if(px==0 && py==0){ phi = 0; return; }

    // get q1 coordinates
    pt_t x =  px; //px>=0 ? px : -px;
    pt_t y =  py; //py>=0 ? py : -py;
    if(px<0) x = -px;
    if(py<0) y = -py;
    // transform so a<b
    pt_t a = x; //x<y ? x : y;
    pt_t b = y; //x<y ? y : x;
    if(a>b){ a = y; b = x; }

    pt_t inv_b;
    if(b>= (1<<(PT_SIZE-DROP_BITS))) inv_b = 1; 
    // don't bother to store these large numbers in the LUT...
    // instead approximate their inverses as 1
    else inv_b = inv_table[b];

    pt_t a_over_b = a * inv_b; // x 2^(PT_SIZE-DROP_BITS)
    ap_uint<ATAN_SIZE> atan_index = a_over_b >> (PT_SIZE-ATAN_SIZE); // keep only most significant bits
    phi = atan_table[atan_index];

    // rotate from (0,pi/4) to full quad1
    if(y>x) phi = (1<<(PHI_SIZE-2)) - phi; //phi = pi/2 - phi
    // other quadrants
    if( px < 0 && py > 0 ) phi = (1<<(PHI_SIZE-1)) - phi;    // Q2 phi = pi - phi
    if( px > 0 && py < 0 ) phi = -phi;                       // Q4 phi = -phi
    if( px < 0 && py < 0 ) phi = -((1<<(PHI_SIZE-1)) - phi); // Q3 composition of both

    return;
}

// General comments:
// division = multiplication and bit shift
// if a, b uint<16>, then a in (0,2^16-1) and 1==2^16
// then 1/b=2^16/b and a/b=a*(2^16/b)
// can convert to decimal by shifting 16 bits
// a/b = a*(2^16/b) >> 16
*/
#endif
