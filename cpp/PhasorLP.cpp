/* PhasorLP.cpp
   Copyright (C) by Bjoern Erlach

   Phasor filter (cascaded complex one pole sections) based lowpass filters.
   Can withstand even the most jumpy audio rate frequency modulation with
   high resonance.  This is the joy of phasor filters.
   shout-out to Max Mathews!
*/

#include <complex.h>
#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>
#include <math.h>
#include "Oversamp.h"
#include <stdio.h>

static InterfaceTable *ft;

#include "OSNonlinearClippers.h"

struct PhasorRLPF2 : public Unit
{
    float m_lastFreq;
    double m_m1r;
    double m_m1i;
    double m_m2r;
    double m_m2i;
};


struct PhasorRLPF4 : public Unit
{
    float m_freq;
    double m_m1r;
    double m_m1i;
    double m_m2r;
    double m_m2i;
    double m_m3r;
    double m_m3i;
    double m_m4r;
    double m_m4i;
};


struct PhasorRLPF4OS4 : public Oversamp2
{
    float m_freq;
    double m_m1r;
    double m_m1i;
    double m_m2r;
    double m_m2i;
    double m_m3r;
    double m_m3i;
    double m_m4r;
    double m_m4i;
};


extern "C"  {

void load(InterfaceTable *inTable);

int api_version(void);
void PhasorRLPF2_Ctor (PhasorRLPF2 *unit);
void PhasorRLPF2_next (PhasorRLPF2 *unit, int inNumSamples);

void PhasorRLPF4_Ctor (PhasorRLPF4 *unit);
void PhasorRLPF4_next (PhasorRLPF4 *unit, int inNumSamples);
void PhasorRLPF4_next_a (PhasorRLPF4 *unit, int inNumSamples);

void PhasorRLPF4OS4_Ctor (PhasorRLPF4OS4 *unit);
void PhasorRLPF4OS4_next (PhasorRLPF4OS4 *unit, int inNumSamples);
void PhasorRLPF4OS4_next_a (PhasorRLPF4OS4 *unit, int inNumSamples);
}


int api_version(void) 
{ 
    return sc_api_version; 
}


void PhasorRLPF2_Ctor (PhasorRLPF2 *unit)
{
    SETCALC(PhasorRLPF2_next);
    unit->m_m1r = 0.0;
    unit->m_m1i = 0.0;
    unit->m_m2r = 0.0;
    unit->m_m2i = 0.0;
    ZOUT0(0) = 0.f;
}


void PhasorRLPF2_next (PhasorRLPF2 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float freq = ZIN0(1);
    float res = ZIN0(2);
    double m1r = unit->m_m1r;
    double m1i = unit->m_m1i;
    double m2r = unit->m_m2r;
    double m2i = unit->m_m2i;

    freq = 2.f*tanf(0.5f*M_PI*freq*SAMPLEDUR);
    // butterworth prototype and freq scale
    double p0r = (-0.70710678 + res) * freq;
    double p0i = 0.70710678 * freq;
    //p0r = res * p0r;
    float gain = sqrtf(p0r*p0r+p0i*p0i);
    gain *= gain;

    double denomr = (2.0-p0r);
    double denomi = -p0i;
    double d = (denomr*denomr+denomi*denomi);
    double g = sqrtf(d);

    p0r = (2.0+p0r);
    double tmp = (p0r*denomr + p0i*denomi) / d;
    p0i = (p0i*denomr - p0r*denomi) / d;
    p0r = tmp;
    gain = gain / (g*g);

    for(int i=0; i<inNumSamples; ++i) {
        float x = ZXP(in);
        double tr, ti;
        double yr, yi;

        tr = x + p0r*m1r - p0i*m1i;
        ti = (p0r*m1i + p0i*m1r);

        yr = tr+m1r;
        yi = ti+m1i;

        m1r = tr;
        m1i = ti;

        tr = yr + p0r*m2r + p0i*m2i;
        ti = yi + (p0r*m2i - p0i*m2r);

        yr = tr+m2r;

        m2r = tr;
        m2i = ti;

        ZXP(out) = yr * gain;

    }

    unit->m_m1r = m1r;
    unit->m_m1i = m1i;
    unit->m_m2r = m2r;
    unit->m_m2i = m2i;
}

//////////////////////////////////////////////////////////////////////

void PhasorRLPF4_Ctor (PhasorRLPF4 *unit)
{
    if (INRATE(1) == calc_FullRate) {
        SETCALC(PhasorRLPF4_next_a);
    } else {
        SETCALC(PhasorRLPF4_next);
    }
    unit->m_m1r = 0.0;
    unit->m_m1i = 0.0;
    unit->m_m2r = 0.0;
    unit->m_m2i = 0.0;
    unit->m_m3r = 0.0;
    unit->m_m3i = 0.0;
    unit->m_m4r = 0.0;
    unit->m_m4i = 0.0;
    unit->m_freq = sc_max(ZIN0(1), 20) * 0.5 * M_PI* SAMPLEDUR;
    ZOUT0(0) = 0.f;
}


void PhasorRLPF4_next (PhasorRLPF4 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float scale = M_PI*SAMPLEDUR;
    float newfreq = ZIN0(1) * scale;
    float res = 1.f-ZIN0(2);
    float clip1 = ZIN0(3);
    float blend1 = ZIN0(4);
    float clip2 = ZIN0(5);
    float blend2 = ZIN0(6);
    float clip3 = ZIN0(7);
    float blend3 = ZIN0(8);
    float clip4 = ZIN0(9);
    float blend4 = ZIN0(10);
    double m1r = unit->m_m1r;
    double m1i = unit->m_m1i;
    double m2r = unit->m_m2r;
    double m2i = unit->m_m2i;
    double m3r = unit->m_m3r;
    double m3i = unit->m_m3i;
    double m4r = unit->m_m4r;
    double m4i = unit->m_m4i;


    float freq = unit->m_freq;
    float freqdiff = (newfreq-freq)/inNumSamples;

    RGET

    for(int i=0; i<inNumSamples; ++i) {

        float w = 2.f*tanf(freq);
        // butterworth prototype and freq scale
        double p0r = (-0.38268343) * w;
        double p0i =  0.92387953 * w;
        double p1r = -0.92387953 * w;
        double p1i =  0.38268343 * w;
        p0r = res * p0r;
        float gain0 = sqrtf(p0r*p0r+p0i*p0i);
        float gain1 = sqrtf(p1r*p1r+p1i*p1i);
        //gain *= gain;
        
        double denomr = (2.0-p0r);
        double denomi = -p0i;
        double d = (denomr*denomr+denomi*denomi);
        //double g = sqrtf(d);
        p0r = (2.0+p0r);
        double tmp = (p0r*denomr + p0i*denomi) / d;
        p0i = (p0i*denomr - p0r*denomi) / d;
        p0r = tmp;
        gain0 /= sqrtf(d);
        
        denomr = (2.0-p1r);
        denomi = -p1i;
        d = (denomr*denomr+denomi*denomi);
        //g = sqrtf(d);
        p1r = (2.0+p1r);
        tmp = (p1r*denomr + p1i*denomi) / d;
        p1i = (p1i*denomr - p1r*denomi) / d;
        p1r = tmp;
        gain1 /= sqrtf(d);

        float x = ZXP(in) * gain0;
        double tr, ti;
        double yr, yi;

        tr = x + p0r*m1r - p0i*m1i;
        ti = (p0r*m1i + p0i*m1r);

        yr = 0.5*(tr+m1r)*gain0;
        yi = 0.5*(ti+m1i)*gain0;

        m1r = tr;
        m1i = ti;

        tr = yr + p0r*m2r + p0i*m2i;
        ti = yi + (p0r*m2i - p0i*m2r);

        yr = 0.5*(tr+m2r)*gain1;
        yi = 0.5*(ti+m2i)*gain1;

        yr += blend2 * (tanh_a(clip2*yr)/clip2 - yr);
        yi += blend1 * (tanh_a(clip1*yi)/clip1 - yi);

        m2r = tr;
        m2i = ti;

        tr = yr + p1r*m3r - p1i*m3i;
        ti = yi + (p1r*m3i + p1i*m3r);

        yr = 0.5*(tr+m3r)*gain1;
        yi = 0.5*(ti+m3i)*gain1;

        yr += blend3 * (tanh_a(clip3*yr)/clip3 - yr);

        m3r = tr;
        m3i = ti;

        tr = yr + p1r*m4r + p1i*m4i;
        ti = yi + (p1r*m4i - p1i*m4r);

        yr = 0.5*(tr+m4r);
        //yi = ti+m2i;
        yr += blend4 * (tanh_a(clip4*yr)/clip4 - yr);

        m4r = tr;
        m4i = ti;
       
        ZXP(out) = yr;
        freq = freq + freqdiff;
    }

    unit->m_freq = freq;
    unit->m_m1r = m1r;
    unit->m_m1i = m1i;
    unit->m_m2r = m2r;
    unit->m_m2i = m2i;
    unit->m_m3r = m3r;
    unit->m_m3i = m3i;
    unit->m_m4r = m4r;
    unit->m_m4i = m4i;
}


void PhasorRLPF4_next_a (PhasorRLPF4 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float scale = M_PI*SAMPLEDUR;
    float *freqIn = ZIN(1);
    float res = 1.f-ZIN0(2);
    float clip1 = ZIN0(3);
    float blend1 = ZIN0(4);
    float clip2 = ZIN0(5);
    float blend2 = ZIN0(6);
    float clip3 = ZIN0(7);
    float blend3 = ZIN0(8);
    float clip4 = ZIN0(9);
    float blend4 = ZIN0(10);
    float noise = ZIN0(11);
    double m1r = unit->m_m1r;
    double m1i = unit->m_m1i;
    double m2r = unit->m_m2r;
    double m2i = unit->m_m2i;
    double m3r = unit->m_m3r;
    double m3i = unit->m_m3i;
    double m4r = unit->m_m4r;
    double m4i = unit->m_m4i;

    RGET

    for(int i=0; i<inNumSamples; ++i) {
        float freq = ZXP(freqIn) * scale;
        float w = 2.f*tanf(freq);
        // butterworth prototype and freq scale
        double p0r = (-0.38268343) * w;
        double p0i =  0.92387953 * w;
        double p1r = -0.92387953 * w;
        double p1i =  0.38268343 * w;
        p0r = res * p0r;
        float gain0 = sqrtf(p0r*p0r+p0i*p0i);
        float gain1 = sqrtf(p1r*p1r+p1i*p1i);
        //gain *= gain;
        
        double denomr = (2.0-p0r);
        double denomi = -p0i;
        double d = (denomr*denomr+denomi*denomi);

        p0r = (2.0+p0r);
        double tmp = (p0r*denomr + p0i*denomi) / d;
        p0i = (p0i*denomr - p0r*denomi) / d;
        p0r = tmp;
        gain0 /= sqrtf(d);
        
        denomr = (2.0-p1r);
        denomi = -p1i;
        d = (denomr*denomr+denomi*denomi);

        p1r = (2.0+p1r);
        tmp = (p1r*denomr + p1i*denomi) / d;
        p1i = (p1i*denomr - p1r*denomi) / d;
        p1r = tmp;
        gain1 /= sqrtf(d);

        float x = ZXP(in) * gain0;
        double tr, ti;
        double yr, yi;

        tr = x + p0r*m1r - p0i*m1i + noise*frand2(s1,s2,s3);
        ti = (p0r*m1i + p0i*m1r) + noise*frand2(s1,s2,s3);

        tr += blend1 * (tanh_a(clip1*tr)/clip1 - tr);
        ti += blend2 * (tanh_a(clip2*ti)/clip2 - ti);

        yr = 0.5*(tr+m1r)*gain0;
        yi = 0.5*(ti+m1i)*gain0;

        m1r = tr;
        m1i = ti;

        tr = yr + p0r*m2r + p0i*m2i + noise*frand2(s1,s2,s3);
        ti = yi + (p0r*m2i - p0i*m2r) + noise*frand2(s1,s2,s3);

        tr += blend3 * (tanh_a(clip3*tr)/clip3 - tr);
        ti += blend4 * (tanh_a(clip4*ti)/clip4 - ti);

        yr = 0.5*(tr+m2r)*gain1;
        yi = 0.5*(ti+m2i)*gain1;

        m2r = tr;
        m2i = ti;

        tr = yr + p1r*m3r - p1i*m3i;
        ti = yi + (p1r*m3i + p1i*m3r);

        yr = 0.5*(tr+m3r)*gain1;
        yi = 0.5*(ti+m3i)*gain1;

        m3r = tr;
        m3i = ti;

        tr = yr + p1r*m4r + p1i*m4i;
        ti = yi + (p1r*m4i - p1i*m4r);

        yr = 0.5*(tr+m4r);

        m4r = tr;
        m4i = ti;
       
        ZXP(out) = yr;
    }

    unit->m_m1r = m1r;
    unit->m_m1i = m1i;
    unit->m_m2r = m2r;
    unit->m_m2i = m2i;
    unit->m_m3r = m3r;
    unit->m_m3i = m3i;
    unit->m_m4r = m4r;
    unit->m_m4i = m4i;
	RPUT
}


//////////////////////////////////////////////////////////////////////

void PhasorRLPF4OS4_Ctor (PhasorRLPF4OS4 *unit)
{
    if (INRATE(1) == calc_FullRate) {
        SETCALC(PhasorRLPF4OS4_next_a);
    } else {
        SETCALC(PhasorRLPF4OS4_next);
    }
    OVERSAMPLE2_INIT;
    unit->m_m1r = 0.0;
    unit->m_m1i = 0.0;
    unit->m_m2r = 0.0;
    unit->m_m2i = 0.0;
    unit->m_m3r = 0.0;
    unit->m_m3i = 0.0;
    unit->m_m4r = 0.0;
    unit->m_m4i = 0.0;
    unit->m_freq = sc_max(ZIN0(1), 20) * 0.5 * M_PI* SAMPLEDUR;
    ZOUT0(0) = 0.f;
}


void PhasorRLPF4OS4_next (PhasorRLPF4OS4 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float scale = M_PI*SAMPLEDUR;
    float newfreq = ZIN0(1) * scale;
    float res = 1.f-ZIN0(2);
    float clip1 = ZIN0(3);
    float blend1 = ZIN0(4);
    float clip2 = ZIN0(5);
    float blend2 = ZIN0(6);
    float clip3 = ZIN0(7);
    float blend3 = ZIN0(8);
    float clip4 = ZIN0(9);
    float blend4 = ZIN0(10);
    float noise = ZIN0(11);
    double m1r = unit->m_m1r;
    double m1i = unit->m_m1i;
    double m2r = unit->m_m2r;
    double m2i = unit->m_m2i;
    double m3r = unit->m_m3r;
    double m3i = unit->m_m3i;
    double m4r = unit->m_m4r;
    double m4i = unit->m_m4i;
    int nsamps = inNumSamples*2;

    float freq = unit->m_freq;
    float freqdiff = (newfreq-freq)/inNumSamples;

    UPSAMPLE2;
    //fprintf(stderr, "%f %f\n", p0r, p0i);
    RGET

    for(int i=0; i<nsamps; ++i) {

        float w = 2.f*tanf(freq);
        // butterworth prototype and freq scale
        double p0r = (-0.38268343) * w;
        double p0i =  0.92387953 * w;
        double p1r = -0.92387953 * w;
        double p1i =  0.38268343 * w;
        p0r = res * p0r;
        float gain0 = sqrtf(p0r*p0r+p0i*p0i);
        float gain1 = sqrtf(p1r*p1r+p1i*p1i);
        //gain *= gain;
        
        double denomr = (2.0-p0r);
        double denomi = -p0i;
        double d = (denomr*denomr+denomi*denomi);
        //double g = sqrtf(d);
        p0r = (2.0+p0r);
        double tmp = (p0r*denomr + p0i*denomi) / d;
        p0i = (p0i*denomr - p0r*denomi) / d;
        p0r = tmp;
        gain0 /= sqrtf(d);
        
        denomr = (2.0-p1r);
        denomi = -p1i;
        d = (denomr*denomr+denomi*denomi);
        //g = sqrtf(d);
        p1r = (2.0+p1r);
        tmp = (p1r*denomr + p1i*denomi) / d;
        p1i = (p1i*denomr - p1r*denomi) / d;
        p1r = tmp;
        gain1 /= sqrtf(d);

        float x = domemoff[i] * gain0;
        double tr, ti;
        double yr, yi;

        tr = x + p0r*m1r - p0i*m1i + noise*frand2(s1,s2,s3);
        ti = (p0r*m1i + p0i*m1r) + noise*frand2(s1,s2,s3);

        tr += blend1 * (tanh_a(clip1*tr)/clip1 - tr);
        ti += blend2 * (tanh_a(clip2*ti)/clip2 - ti);

        yr = 0.5*(tr+m1r)*gain0;
        yi = 0.5*(ti+m1i)*gain0;

        m1r = tr;
        m1i = ti;

        tr = yr + p0r*m2r + p0i*m2i + noise*frand2(s1,s2,s3);
        ti = yi + (p0r*m2i - p0i*m2r) + noise*frand2(s1,s2,s3);

        tr += blend3 * (tanh_a(clip3*tr)/clip3 - tr);
        ti += blend4 * (tanh_a(clip4*ti)/clip4 - ti);

        yr = 0.5*(tr+m2r)*gain1;
        yi = 0.5*(ti+m2i)*gain1;

        yr += blend3 * (tanh_a(clip3*yr)/clip3 - yr);

        m2r = tr;
        m2i = ti;

        tr = yr + p1r*m3r - p1i*m3i;
        ti = yi + (p1r*m3i + p1i*m3r);

        yr = 0.5*(tr+m3r)*gain1;
        yi = 0.5*(ti+m3i)*gain1;

        m3r = tr;
        m3i = ti;

        tr = yr + p1r*m4r + p1i*m4i;
        ti = yi + (p1r*m4i - p1i*m4r);

        yr = 0.5*(tr+m4r);

        m4r = tr;
        m4i = ti;
        domemoff[i] = yr;
        freq = freq + freqdiff;
    }

    DOWNSAMPLE2STD;
    unit->m_freq = freq;
    unit->m_m1r = m1r;
    unit->m_m1i = m1i;
    unit->m_m2r = m2r;
    unit->m_m2i = m2i;
    unit->m_m3r = m3r;
    unit->m_m3i = m3i;
    unit->m_m4r = m4r;
    unit->m_m4i = m4i;
}


void PhasorRLPF4OS4_next_a (PhasorRLPF4OS4 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float scale = M_PI*SAMPLEDUR;
    float *freqIn = ZIN(1);
    float res = 1.f-ZIN0(2);
    float clip1 = ZIN0(3);
    float blend1 = ZIN0(4);
    float clip2 = ZIN0(5);
    float blend2 = ZIN0(6);
    float clip3 = ZIN0(7);
    float blend3 = ZIN0(8);
    float clip4 = ZIN0(9);
    float blend4 = ZIN0(10);
    float noise = ZIN0(11);
    double m1r = unit->m_m1r;
    double m1i = unit->m_m1i;
    double m2r = unit->m_m2r;
    double m2i = unit->m_m2i;
    double m3r = unit->m_m3r;
    double m3i = unit->m_m3i;
    double m4r = unit->m_m4r;
    double m4i = unit->m_m4i;
    int nsamps  = inNumSamples*2;

    RGET

    UPSAMPLE2;
    float freq;

    for(int i=0; i<nsamps; ++i) {
        if (!(i&1)) {
            freq = ZXP(freqIn) * scale;
        }
        float w = 2.f*tanf(freq);
        // butterworth prototype and freq scale
        double p0r = (-0.38268343) * w;
        double p0i =  0.92387953 * w;
        double p1r = -0.92387953 * w;
        double p1i =  0.38268343 * w;
        p0r = res * p0r;
        float gain0 = sqrtf(p0r*p0r+p0i*p0i);
        float gain1 = sqrtf(p1r*p1r+p1i*p1i);
        //gain *= gain;
        
        double denomr = (2.0-p0r);
        double denomi = -p0i;
        double d = (denomr*denomr+denomi*denomi);
        //double g = sqrtf(d);
        p0r = (2.0+p0r);
        double tmp = (p0r*denomr + p0i*denomi) / d;
        p0i = (p0i*denomr - p0r*denomi) / d;
        p0r = tmp;
        gain0 /= sqrtf(d);
        
        denomr = (2.0-p1r);
        denomi = -p1i;
        d = (denomr*denomr+denomi*denomi);
        //g = sqrtf(d);
        p1r = (2.0+p1r);
        tmp = (p1r*denomr + p1i*denomi) / d;
        p1i = (p1i*denomr - p1r*denomi) / d;
        p1r = tmp;
        gain1 /= sqrtf(d);

        float x = domemoff[i] * gain0;
        double tr, ti;
        double yr, yi;

        tr = x + p0r*m1r - p0i*m1i + noise*frand2(s1,s2,s3);
        ti = (p0r*m1i + p0i*m1r) + noise*frand2(s1,s2,s3);

        tr += blend1 * (tanh_a(clip1*tr)/clip1 - tr);
        ti += blend2 * (tanh_a(clip2*ti)/clip2 - ti);

        yr = 0.5*(tr+m1r)*gain0;
        yi = 0.5*(ti+m1i)*gain0;

        m1r = tr;
        m1i = ti;

        tr = yr + p0r*m2r + p0i*m2i + noise*frand2(s1,s2,s3);
        ti = yi + (p0r*m2i - p0i*m2r) + noise*frand2(s1,s2,s3);

        tr += blend3 * (tanh_a(clip3*tr)/clip3 - tr);
        ti += blend4 * (tanh_a(clip4*ti)/clip4 - ti);

        yr = 0.5*(tr+m2r)*gain1;
        yi = 0.5*(ti+m2i)*gain1;

        yr += blend3 * (tanh_a(clip3*yr)/clip3 - yr);

        m2r = tr;
        m2i = ti;

        tr = yr + p1r*m3r - p1i*m3i;
        ti = yi + (p1r*m3i + p1i*m3r);

        yr = 0.5*(tr+m3r)*gain1;
        yi = 0.5*(ti+m3i)*gain1;

        m3r = tr;
        m3i = ti;

        tr = yr + p1r*m4r + p1i*m4i;
        ti = yi + (p1r*m4i - p1i*m4r);

        yr = 0.5*(tr+m4r);

        m4r = tr;
        m4i = ti;
       
        domemoff[i] = yr;
    }

    DOWNSAMPLE2STD;

    unit->m_m1r = m1r;
    unit->m_m1i = m1i;
    unit->m_m2r = m2r;
    unit->m_m2i = m2i;
    unit->m_m3r = m3r;
    unit->m_m3i = m3i;
    unit->m_m4r = m4r;
    unit->m_m4i = m4i;
    RPUT
}


void
load(InterfaceTable *inTable)
{
  ft = inTable;
  DefineSimpleUnit(PhasorRLPF2);
  DefineSimpleUnit(PhasorRLPF4);
  DefineSimpleUnit(PhasorRLPF4OS4);
}
