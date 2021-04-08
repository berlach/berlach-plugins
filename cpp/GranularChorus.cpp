/* GranularChorus.cpp
   Copyright (C) by Bjoern Erlach

   Chorus mainly useful for monophonic inputs with known frequency.
*/

#include "SC_PlugIn.h"
#include <stdio.h>

static InterfaceTable *ft;


struct GranularChorus : public Unit
{
    double m_delayphs[8192];
    float  m_delay[8192];
    double m_oldinpinc;
    double m_oldoutpinc;
    double m_outphs;
    double m_inphs;
    int m_idx;
    double m_outidxf;
    double m_outidxf2;
    float m_fade;
    int m_trigready;
};


 extern "C"
 {
     void load(InterfaceTable *inTable);
     int api_version(void);

     void GranularChorus_Ctor(GranularChorus *unit);
     void GranularChorus_next(GranularChorus *unit, int inNumSamples);
}


 int api_version(void) 
{ 
    return sc_api_version; 
}


void GranularChorus_Ctor (GranularChorus *unit)
{
    float freq = ZIN0(1);
    float outfreq = ZIN0(2);

    SETCALC(GranularChorus_next);

    unit->m_idx = 0;
    unit->m_oldinpinc = freq/SAMPLERATE;
    unit->m_oldoutpinc = outfreq/SAMPLERATE;
    unit->m_outidxf = 0.0;
    unit->m_outidxf2 = SAMPLERATE*3.f/freq;
    unit->m_inphs = 0.0;
    unit->m_fade = 0.f;
    unit->m_trigready = 1;
    memset(unit->m_delay, 0, sizeof(float) * 8192);
    memset(unit->m_delayphs, 0, sizeof(double) * 8192);
    ZOUT0(0) = 0.f;
}


void GranularChorus_next (GranularChorus *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float freq = ZIN0(1);
    float outfreq = ZIN0(2);
    float *delay = unit->m_delay;
    double *delayphs = unit->m_delayphs;
    int idx = unit->m_idx;
    double inpinc = unit->m_oldinpinc;
    double newinpinc = freq/(double)SAMPLERATE;
    double inpincinc = (newinpinc-inpinc) / inNumSamples;
    double outpinc = unit->m_oldoutpinc;
    double newoutpinc = outfreq/(double)SAMPLERATE;
    double outpincinc = (newoutpinc-outpinc) / inNumSamples;
    double inphs = unit->m_inphs;
    double ratio;
    double outidxf = unit->m_outidxf;
    double outidxf2 = unit->m_outidxf2;
    // fundamental period out input signal in samples
    double period = 1.0 / (0.5 * (newinpinc + inpinc));
    float fade = unit->m_fade;

    unit->m_oldinpinc = newinpinc;
    unit->m_oldoutpinc = newoutpinc;
    ratio = outpinc/inpinc;    

    if (idx > 8192) {
        idx -= 8192;
        outidxf -= 8192.0;
        outidxf2 -= 8192.0;
    }

    for (int k=0; k<inNumSamples; ++k) {
        float x = ZXP(in);
        float ff;
        delay[idx&8191] = x;
        delayphs[idx&8191] = inphs;
        
        int outidx = (int) outidxf;
        float frac = outidxf - (float) outidx;
        int outidx2 = (int) outidxf2;
        float frac2 = outidxf2 - (float) outidx2;


        if (fade>=2.f) {
            // start grain 1
            double outphs = (1.f-frac) * delayphs[outidx&8191] + frac * delayphs[(outidx+1)&8191];
            int dist = idx - outidx;
            int ideal = (int)(period*2);
            if (dist > ideal) {
                outidxf += period;  
                dist = idx - (int) outidxf;
                if (dist > ideal) {
                    outidxf += period;  
                }
            } else {
                outidxf -= period;  
                dist = idx - (int)outidxf;
                if (dist < ideal) {
                    outidxf -= period;  
                }
            }
            outidx = (int) outidxf;
            double phs1 = delayphs[(outidx)&8191];
            double phs2 = delayphs[(outidx+1)&8191];
            if ((phs1 > outphs) || (phs2 < outphs) || ((phs2-phs1) < 0.00003)) {
                frac = outidxf - (float) outidx;
            } else {
                frac = (outphs - phs1) / (phs2-phs1);
            }
            outidxf = ((double) outidx) + frac;
            unit->m_trigready = 1;
            ff = fade = 0.f;
        } else if (fade>=1.f) { 
            if (unit->m_trigready) {
            // start grain 2
                unit->m_trigready = 0;
                double outphs = (1.f-frac2) * delayphs[outidx2&8191] + frac2 * delayphs[(outidx2+1)&8191];
                int dist = idx - outidx2;
                int ideal = (int)(period*3);
                if (dist > ideal) {
                    outidxf2 += period;  
                    dist = idx - (int) outidxf2;
                    if (dist > ideal) {
                        outidxf2 += period;  
                    }
                } else {
                    outidxf2 -= period;  
                    dist = idx - (int)outidxf2;
                    if (dist < ideal) {
                        outidxf2 -= period;  
                    }
                }
                outidx2 = (int) outidxf2;
                double phs1 = delayphs[(outidx2)&8191];
                double phs2 = delayphs[(outidx2+1)&8191];
                if ((phs1 >= outphs) || (phs2 < outphs) || ((phs2-phs1) < 0.00003)) {
                    frac2 = outidxf2 - (float) outidx2;
                } else {
                    frac2 = (outphs - phs1) / (phs2-phs1);
                }
                outidxf2 = ((double) outidx2) + frac2;
            }
            ff = 2.f-fade; 
        } else { 
            ff = fade; 
        }
        
        float x10 = delay[(outidx)&8191];
        float x11 = delay[(outidx+1)&8191];
        float outsmp = frac * x11 + (1.f-frac) * x10;

        float x20 = delay[(outidx2)&8191];
        float x21 = delay[(outidx2+1)&8191];
        float outsmp2 = frac2 * x21 + (1.f-frac2) * x20;

        ZXP(out) = ff * outsmp + (1.f - ff) * outsmp2;

        fade += 0.0016;
        outpinc += outpincinc;
        inpinc += inpincinc;
        inphs += inpinc;
        if (inphs>=1.0) { inphs -= 1.0; }
        ratio = outpinc/inpinc;
        idx++;
        outidxf += ratio;
        outidxf2 += ratio;
    }

    unit->m_fade = fade;
    unit->m_outidxf = outidxf;
    unit->m_outidxf2 = outidxf2;
    unit->m_idx = idx;
    unit->m_inphs = inphs;
}


void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(GranularChorus);
}
