/* ReverbUGens.cpp
   Copyright (C) by Bjoern Erlach

   Collection of standard building blocks for reverbs.

   Contains variations of Schroeder Allpass blocks
   different feedback-delay-networks
   velvet noise convolver
   tapped delay lines   
   nested allpass filters
   + some less common lossy allpass filters

   Units are designed to allow to reuse delay line outputs wherever possible
   to build dense reverbs efficiently.
*/

#include "SC_PlugIn.h"
#include <boost/align/is_aligned.hpp>
#include <stdio.h>

static InterfaceTable *ft;


struct AllpassLPN : public Unit
{
    static const int minDelaySamples = 1;
    float *m_dlybuf;
    float m_y1;
    float m_feedbk, m_decaytime;
    float m_dsamp, m_fdelaylen;
    float m_delaytime, m_maxdelaytime;
    long m_iwrphase, m_idelaylen, m_mask;
    long m_numoutput;
};


struct SimpleFDN2 : public Unit
{
    float m_mem1[4096*8];
    float m_mem2[4096*8];
    int m_mask;
    int m_idx;
};


struct SimpleFDN4 : public Unit
{
    float m_mem1[4096*8];
    float m_mem2[4096*8];
    float m_mem3[4096*8];
    float m_mem4[4096*8];
    int m_mask;
    int m_idx;
};


struct JotFDN3 : public Unit
{
    float m_mem1[4096*8];
    float m_mem2[4096*8];
    float m_mem3[4096*8];
    int m_mask;
    int m_idx;
};


struct JotFDN3LPAP : public Unit
{
    float m_mem1[4096*8];
    float m_mem2[4096*8];
    float m_mem3[4096*8];
    float m_mem4[4096*8];
    float m_mem5[4096*8];
    int m_mask;
    int m_idx;
    float m_y1;
};


struct Velvet : public Unit
{
    float *m_mem1;
    int m_postaps[256];
    int m_negtaps[256];
    int m_nnegtaps;
    int m_npostaps;
    int m_mask;
    int m_idx;
};


struct Velvet1 : public Unit
{
    float *m_mem1;
    int m_postaps[256];
    int m_npostaps;
    int m_mask;
    int m_idx;
};


struct VelvetTaps : public Unit
{
    SndBuf *m_buf;    
    int m_postaps[256];
    int m_negtaps[256];
    int m_nnegtaps;
    int m_npostaps;
    int m_mask;
    int m_idx;
};


struct BufNestedLBCF : public Unit
{
    SndBuf *m_buf;
    float *m_mem2;
    float m_apy1[16];
    float m_lpy1;
    float m_firx1;
    float m_firx2;
    int m_mask2;
    int m_idx;
};


struct BufLBCF : public Unit
{
    SndBuf *m_buf;
    float m_apy1[16];
    float m_lpy1;
    int m_mask;
    int m_idx;
};


struct BufLAPF : public Unit
{
    SndBuf *m_buf;
    float m_apy1[16];
    float m_lpy1;
    int m_mask;
    int m_idx;
};


struct BufOneTap : public Unit
{
    float m_fbufnum;
    SndBuf *m_buf;
    int m_idx;
    int m_mask;    
};


struct BufOneTapC : public Unit
{
    float m_fbufnum;
    SndBuf *m_buf;
    float m_lastdel;
    int m_idx;
    int m_mask;    
};


struct BufWrite : public Unit
{
    float m_fbufnum;
    SndBuf *m_buf;
    int m_idx;
    int m_mask;    
};


struct LBCF : public Unit
{
    float *m_mem1;
    float m_apy1[16];
    float m_lpy1;
    float m_firx1;
    float m_firx2;
    int m_mask;
    int m_idx;
};


struct NestedLBCF : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float m_apy1[16];
    float m_lpy1;
    float m_firx1;
    float m_firx2;
    int m_mask1;
    int m_mask2;
    int m_idx;
};


struct AllpassLatc : public Unit
{
    float *m_mem1[16];
    int m_nsect;
    int *m_mask[16];    
    //float m_apy1[16];
    //float m_lpy1;
    int m_idx;
};


struct AllpassLatcLP : public Unit
{
    float *m_mem1[16];
    int m_nsect;
    int *m_mask[16];
    float m_y1[16];
    //float m_apy1[16];
    //float m_lpy1;
    int m_idx;
};


struct Hadamard8Delay : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float *m_mem3;
    float *m_mem4;
    float *m_mem5;
    float *m_mem6;
    float *m_mem7;
    int m_mask;
    int m_idx;
};


struct Hadamard4Delay : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float *m_mem3;
    int m_mask;
    int m_idx;
};


struct Ortho3Delay : public Unit
{
    float *m_mem1;
    float *m_mem2;
    int m_mask;
    int m_idx;
    float m_m0, m_m1, m_m2, m_m3, m_m4, m_m5, m_m6, m_m7, m_m8;
};


struct Nested2AllpassN : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float *m_mem3;
    float m_y1;
    int m_mask1;
    int m_mask2;
    int m_mask3;
    int m_idx;
};


struct Allpass2N : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float m_freq;
    float m_rq;
    int m_mask;
    int m_idx;
    float m_b1;
    float m_a0;
};


struct NestedAllpass2N : public Unit
{
    float *m_mem1;
    float *m_mem2;
    float *m_mem3;
    float m_freq;
    float m_rq;
    int m_mask;
    int m_mask3;
    int m_idx;
    float m_b1;
    float m_a0;
};


struct FbTapper : Unit {
    float m_lpy1;
    float *m_mem1;
    int m_idx;
    int m_mask;
};



extern "C"
{
    void load(InterfaceTable *inTable);
    int api_version(void);

    void FbTapper_Ctor(FbTapper *unit);
    void FbTapper_Dtor(FbTapper *unit);
    void FbTapper_next(FbTapper *unit, int inNumSamples);

    void AllpassLPN_Ctor(AllpassLPN *unit);
    void AllpassLPN_next(AllpassLPN *unit, int inNumSamples);
    void AllpassLPN_next_z(AllpassLPN *unit, int inNumSamples);
    void AllpassLPN_next_a(AllpassLPN *unit, int inNumSamples);
    void AllpassLPN_next_a_z(AllpassLPN *unit, int inNumSamples);

    void Allpass2N_Ctor(Allpass2N *unit);
    void Allpass2N_Dtor(Allpass2N *unit);
    void Allpass2N_next(Allpass2N *unit, int inNumSamples);

    void NestedAllpass2N_Ctor(NestedAllpass2N *unit);
    void NestedAllpass2N_Dtor(NestedAllpass2N *unit);
    void NestedAllpass2N_next(NestedAllpass2N *unit, int inNumSamples);

    void Nested2AllpassN_Ctor(Nested2AllpassN *unit);
    void Nested2AllpassN_Dtor(Nested2AllpassN *unit);
    void Nested2AllpassN_next(Nested2AllpassN *unit, int inNumSamples);

    void LBCF_Ctor(LBCF *unit);
    void LBCF_next(LBCF *unit, int inNumSamples);
    void LBCF_next_fir(LBCF *unit, int inNumSamples);

    void NestedLBCF_Ctor(NestedLBCF *unit);
    void NestedLBCF_next(NestedLBCF *unit, int inNumSamples);
    void NestedLBCF_next_fir(NestedLBCF *unit, int inNumSamples);

    void BufNestedLBCF_Ctor(BufNestedLBCF *unit);
    void BufNestedLBCF_next(BufNestedLBCF *unit, int inNumSamples);
    void BufNestedLBCF_next_fir(BufNestedLBCF *unit, int inNumSamples);

    void BufOneTap_Ctor(BufOneTap *unit);
    void BufOneTap_next(BufOneTap *unit, int inNumSamples);

    void BufOneTapC_Ctor(BufOneTapC *unit);
    void BufOneTapC_next(BufOneTapC *unit, int inNumSamples);

    void Hadamard8Delay_Ctor(Hadamard8Delay *unit);
    void Hadamard8Delay_next(Hadamard8Delay *unit, int inNumSamples);

    void Hadamard4Delay_Ctor(Hadamard4Delay *unit);
    void Hadamard4Delay_next(Hadamard4Delay *unit, int inNumSamples);

    void Ortho3Delay_Ctor(Ortho3Delay *unit);
    void Ortho3Delay_next(Ortho3Delay *unit, int inNumSamples);

    void BufWrite_Ctor(BufWrite *unit);
    void BufWrite_next(BufWrite *unit, int inNumSamples);

    void AllpassLatc_Ctor(AllpassLatc *unit);
    void AllpassLatc_next(AllpassLatc *unit, int inNumSamples);
    void AllpassLatc_Dtor(AllpassLatc *unit);

    void AllpassLatcLP_Ctor(AllpassLatcLP *unit);
    void AllpassLatcLP_next(AllpassLatcLP *unit, int inNumSamples);
    void AllpassLatcLP_Dtor(AllpassLatcLP *unit);
    
    void BufLBCF_Ctor(BufLBCF *unit);
    void BufLBCF_next(BufLBCF *unit, int inNumSamples);

    void BufLAPF_Ctor(BufLAPF *unit);
    void BufLAPF_next(BufLAPF *unit, int inNumSamples);

    void SimpleFDN2_Ctor(SimpleFDN2 *unit);
    void SimpleFDN2_next(SimpleFDN2 *unit, int inNumSamples);

    void SimpleFDN4_Ctor(SimpleFDN4 *unit);
    void SimpleFDN4_next(SimpleFDN4 *unit, int inNumSamples);

    void JotFDN3_Ctor(JotFDN3 *unit);
    void JotFDN3_next(JotFDN3 *unit, int inNumSamples);

    void JotFDN3LPAP_Ctor(JotFDN3LPAP *unit);
    void JotFDN3LPAP_next(JotFDN3LPAP *unit, int inNumSamples);

    void Velvet_Ctor(Velvet *unit);
    void Velvet_Dtor(Velvet *unit);
    void Velvet_next(Velvet *unit, int inNumSamples);

    void Velvet1_Ctor(Velvet1 *unit);
    void Velvet1_Dtor(Velvet1 *unit);
    void Velvet1_next(Velvet1 *unit, int inNumSamples);

    void VelvetTaps_Ctor(VelvetTaps *unit);
    void VelvetTaps_next(VelvetTaps *unit, int inNumSamples);
}

int api_version(void) 
{ 
    return sc_api_version; 
}


#define CTOR_GET_BUF                                    \
    float fbufnum  = ZIN0(0);                           \
    fbufnum = sc_max(0.f, fbufnum);                     \
    uint32 bufnum = (int)fbufnum;                       \
    World *world = unit->mWorld;                        \
    SndBuf *buf;                                        \
    if (bufnum >= world->mNumSndBufs) {                 \
        int localBufNum = bufnum - world->mNumSndBufs;  \
        Graph *parent = unit->mParent;                  \
        if(localBufNum <= parent->localBufNum) {        \
            buf = parent->mLocalSndBufs + localBufNum;  \
        } else {                                        \
            bufnum = 0;                                 \
            buf = world->mSndBufs + bufnum;             \
        }                                               \
    } else {                                            \
        buf = world->mSndBufs + bufnum;                 \
    }


#define CHECK_BUF                               \
    if (!bufData) {                             \
        unit->mDone = true;                     \
        ClearUnitOutputs(unit, inNumSamples);   \
        return;                                 \
    }


template <typename Unit>
static float CalcDelay(Unit *unit, float delaytime)
{
    float minDelay = Unit::minDelaySamples;
    float next_dsamp = delaytime * (float)SAMPLERATE;
    return sc_clip(next_dsamp, minDelay, unit->m_fdelaylen);
}


static bool AllpassLPN_AllocDelayLine(AllpassLPN *unit, const char * className)
{
    long delaybufsize = (long)ceil(unit->m_maxdelaytime * SAMPLERATE + 1.f);
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    unit->m_fdelaylen = unit->m_idelaylen = delaybufsize;

    if (unit->m_dlybuf)
        RTFree(unit->mWorld, unit->m_dlybuf);
    unit->m_dlybuf = (float*)RTAlloc(unit->mWorld, delaybufsize * sizeof(float));

#if 0 // for debugging we may want to fill the buffer with nans
    std::fill_n(unit->m_dlybuf, delaybufsize, std::numeric_limits<float>::signaling_NaN());
#endif
    memset(unit->m_dlybuf, 0, sizeof(float) * delaybufsize);
        
    if (unit->m_dlybuf == NULL) {
        SETCALC(ft->fClearUnitOutputs);
        ClearUnitOutputs(unit, 1);

        if(unit->mWorld->mVerbosity > -2)
            Print("Failed to allocate memory for %s ugen.\n", className);
    }

    unit->m_mask = delaybufsize - 1;
    return (unit->m_dlybuf != NULL);
}


template <typename Unit>
static bool DelayUnit_Reset(Unit *unit, const char * className)
{
    unit->m_maxdelaytime = ZIN0(1);
    unit->m_delaytime = ZIN0(2);
    unit->m_dlybuf = 0;

    if (!AllpassLPN_AllocDelayLine(unit, className))
        return false;

    unit->m_dsamp = CalcDelay(unit, unit->m_delaytime);

    unit->m_numoutput = 0;
    unit->m_iwrphase = 0;
    return true;
}

template <typename Unit>
static bool FeedbackDelay_Reset(Unit *unit, const char * className)
{
    unit->m_decaytime = ZIN0(3);

    bool allocationSucessful = DelayUnit_Reset(unit, className);
    if (!allocationSucessful)
        return false;

    //unit->m_feedbk = sc_CalcFeedback(unit->m_delaytime, unit->m_decaytime);
    unit->m_feedbk = unit->m_decaytime;
    return true;
}



void AllpassLPN_Ctor(AllpassLPN *unit)
{
    bool allocationSucessful = FeedbackDelay_Reset(unit, "AllpassLPN");
    if (!allocationSucessful)
        return;

    if(INRATE(2) == calc_FullRate)
        SETCALC(AllpassLPN_next);
    else
        SETCALC(AllpassLPN_next);
    unit->m_y1 = 0.f;
    ZOUT0(0) = 0.f;
}


void AllpassLPN_Dtor(AllpassLPN *unit)
{
    RTFree(unit->mWorld, unit->m_dlybuf);
}


void AllpassLPN_next(AllpassLPN *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    const float *in = ZIN(0);
    float delaytime = ZIN0(2);
    float decaytime = ZIN0(3);
    float absorb = ZIN0(4);
    float lp = ZIN0(5);
    float lp1 = (1.f-lp) * absorb;
    float y1 = unit->m_y1;

    float *dlybuf = unit->m_dlybuf;
    long iwrphase = unit->m_iwrphase;
    float dsamp = unit->m_dsamp;
    float feedbk = unit->m_feedbk;
    long mask = unit->m_mask;

    if (delaytime == unit->m_delaytime) {
        long irdphase = iwrphase - (long)dsamp;
        float* dlybuf1 = dlybuf - ZOFF;
        float* dlyrd   = dlybuf1 + (irdphase & mask);
        float* dlywr   = dlybuf1 + (iwrphase & mask);
        float* dlyN    = dlybuf1 + unit->m_idelaylen;
        if (decaytime == unit->m_decaytime) {
            long remain = inNumSamples;
            while (remain) {
                long rdspace = dlyN - dlyrd;
                long wrspace = dlyN - dlywr;
                long nsmps = sc_min(rdspace, wrspace);
                nsmps = sc_min(remain, nsmps);
                remain -= nsmps;
                LOOP(nsmps,
                     float inp = ZXP(in);
                     float value = ZXP(dlyrd);
                     float dwr = value;                                     
                     y1 = lp1 * dwr + lp * y1;
                     ZXP(dlywr) = inp - feedbk * y1;
                     ZXP(out) = y1 + feedbk * inp;
                     );
                if (dlyrd == dlyN) dlyrd = dlybuf1;
                if (dlywr == dlyN) dlywr = dlybuf1;
            }
        } else {
            //float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
            float next_feedbk = decaytime;
            float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
            long remain = inNumSamples;
            while (remain) {
                long rdspace = dlyN - dlyrd;
                long wrspace = dlyN - dlywr;
                long nsmps = sc_min(rdspace, wrspace);
                nsmps = sc_min(remain, nsmps);
                remain -= nsmps;

                LOOP(nsmps,
                     float inp = ZXP(in);
                     float value = ZXP(dlyrd);
                     float dwr = value;                                     
                     y1 = lp1 * dwr + lp * y1;
                     ZXP(dlywr) = inp - feedbk * y1;
                     ZXP(out) = y1 + feedbk * inp;
                     feedbk += feedbk_slope;
                     );
                if (dlyrd == dlyN) dlyrd = dlybuf1;
                if (dlywr == dlyN) dlywr = dlybuf1;
            }
            unit->m_feedbk = feedbk;
            unit->m_decaytime = decaytime;
        }
        iwrphase += inNumSamples;
    } else {
        float next_dsamp = CalcDelay(unit, delaytime);
        float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);
        //float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
        float next_feedbk = decaytime;
        float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);
        LOOP1(inNumSamples,
              dsamp += dsamp_slope;
              feedbk += feedbk_slope;
              float inp = ZXP(in);
              long irdphase = iwrphase - (long)dsamp;
              float value = dlybuf[irdphase & mask];
              float dwr = value;
              y1 = lp1 * dwr + lp * y1;
              dlybuf[iwrphase & mask] = -feedbk * y1 + inp;
              ZXP(out) = y1 + feedbk * inp;
              ++iwrphase;
              );
        unit->m_feedbk = feedbk;
        unit->m_dsamp = dsamp;
        unit->m_delaytime = delaytime;
        unit->m_decaytime = decaytime;
    }
    unit->m_y1 = y1;
    unit->m_iwrphase = iwrphase;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////

void SimpleFDN2_Ctor(SimpleFDN2 *unit)
{
    SETCALC(SimpleFDN2_next);
    unit->m_idx = 0;
    unit->m_mask = 4096*8-1;
    memset(unit->m_mem1, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem2, 0, sizeof(float) * 4096*8);
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    ZOUT0(0) = 0.0f;
}


void SimpleFDN2_next(SimpleFDN2 *unit, int inNumSamples)
{
    float *in1 = ZIN(0);
    float *in2 = ZIN(1);
    float *out1 = ZOUT(0);
    float *out2 = ZOUT(1);
    float sr = (float) SAMPLERATE;
    int d1 = (int) (ZIN0(2) * sr);
    int d2 = (int) (ZIN0(3) * sr);
    float a1 = ZIN0(4);
    float a2 = ZIN0(5);
    float a3 = ZIN0(6);
    float a4 = ZIN0(7);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;

    LOOP(inNumSamples,
         float x1;
         float x2;
         float tmp;
         float xin1;
         float xin2;
         xin1 = ZXP(in1);
         xin2 = ZXP(in2); 
         x1 = mem1[(idx-d1)&mask];
         x2 = mem2[(idx-d2)&mask];
         tmp = a1 * x1 + a2 * x2 + xin1;
         x2 = a3 * x1 + a4 * x2 + xin2;
         mem1[idx&mask] = tmp;
         mem2[idx&mask] = x2;  
         ZXP(out1) = x1;
         ZXP(out2) = x2;
         idx++;
         );
    

    unit->m_idx = idx;
}


//////////////////////////////////////////////////////////////////////////


void SimpleFDN4_Ctor(SimpleFDN4 *unit)
{
    SETCALC(SimpleFDN4_next);
    unit->m_idx = 0;
    unit->m_mask = 4096*8-1;
    memset(unit->m_mem1, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem2, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem3, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem4, 0, sizeof(float) * 4096*8);
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    ZOUT0(0) = 0.0f;
}


void SimpleFDN4_next(SimpleFDN4 *unit, int inNumSamples)
{
    float *in1 = ZIN(0);
    float *in2 = ZIN(1);
    float *in3 = ZIN(2);
    float *in4 = ZIN(3);
    float *out1 = ZOUT(0);
    float *out2 = ZOUT(1);
    float *out3 = ZOUT(2);
    float *out4 = ZOUT(3);
    float sr = (float) SAMPLERATE;
    int d1 = (int) (ZIN0(4) * sr);
    int d2 = (int) (ZIN0(5) * sr);
    int d3 = (int) (ZIN0(6) * sr);
    int d4 = (int) (ZIN0(7) * sr);
    float g = ZIN0(8);
    int type = (int) ZIN0(9);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    float *mem4 = unit->m_mem4;
    int idx = unit->m_idx;

    switch(type) {
    case 0:
        // stautner puckette fdn
        LOOP(inNumSamples,
             float x1;
             float x2;
             float x3;
             float x4;
             float xin1;
             float xin2;
             float xin3;
             float xin4;
             float y1;
             float y2;
             float y3;
             float y4;
             xin1 = ZXP(in1);
             xin2 = ZXP(in2); 
             xin3 = ZXP(in3); 
             xin4 = ZXP(in4); 
             x1 = g * mem1[(idx-d1)&mask];
             x2 = g * mem2[(idx-d2)&mask];
             x3 = g * mem3[(idx-d3)&mask];
             x4 = g * mem4[(idx-d4)&mask];
             y1 = x2+x3 + xin1;
             y2 = -x1-x4 + xin2;
             y3 = x1-x4 + xin3;
             y4 = x2-x3 + xin4;
          
             mem1[idx&mask] = y1;
             mem2[idx&mask] = y2;  
             mem3[idx&mask] = y3;  
             mem4[idx&mask] = y4;  
             ZXP(out1) = y1;
             ZXP(out2) = y2;
             ZXP(out3) = y3;
             ZXP(out4) = y4;

             idx++;
             );
        break;
    case 1:
        // stautner puckette fdn
        LOOP(inNumSamples,
             float x1;
             float x2;
             float x3;
             float x4;
             float xin1;
             float xin2;
             float xin3;
             float xin4;
             float y1;
             float y2;
             float y3;
             float y4;
             xin1 = ZXP(in1);
             xin2 = ZXP(in2); 
             xin3 = ZXP(in3); 
             xin4 = ZXP(in4); 
             x1 = g * mem1[(idx-d1)&mask];
             x2 = g * mem2[(idx-d2)&mask];
             x3 = g * mem3[(idx-d3)&mask];
             x4 = g * mem4[(idx-d4)&mask];
             float x1x2 = x1+x2;
             float x1mx2 = x1-x2;
             y1 = x1x2+x3+x4+xin1;
             y2 = x1mx2+x3-x4+xin2;
             y3 = x1x2-x3-x4+xin3;
             y4 = x1mx2-x3+x4+xin4;
             mem1[idx&mask] = y1;
             mem2[idx&mask] = y2;  
             mem3[idx&mask] = y3;  
             mem4[idx&mask] = y4;  
             ZXP(out1) = y1;
             ZXP(out2) = y2;
             ZXP(out3) = y3;
             ZXP(out4) = y4;
             idx++;
             );
        break;
    }
    

    unit->m_idx = idx;
}





//////////////////////////////////////////////////////////////////////////



void JotFDN3_Ctor(JotFDN3 *unit)
{
    SETCALC(JotFDN3_next);
    unit->m_idx = 0;
    unit->m_mask = 4096*8-1;
    memset(unit->m_mem1, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem2, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem3, 0, sizeof(float) * 4096*8);
    ZOUT0(0) = 0.0f;
}


void JotFDN3_next(JotFDN3 *unit, int inNumSamples)
{
    float *in1 = ZIN(0);
    float *in2 = ZIN(1);
    float *in3 = ZIN(2);
    float *out1 = ZOUT(0);
    float *out2 = ZOUT(1);
    float *out3 = ZOUT(2);
    float sr = (float) SAMPLERATE;
    int d1 = (int) (ZIN0(3) * sr);
    int d2 = (int) (ZIN0(4) * sr);
    int d3 = (int) (ZIN0(5) * sr);
    float l1 = ZIN0(6);
    float l2 = ZIN0(7);
    float l3 = ZIN0(8);
    float a = ZIN0(9);
    float b = ZIN0(10);
    float c = ZIN0(11);
    float g = ZIN0(12);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    int idx = unit->m_idx;

    LOOP(inNumSamples,
         float x1;
         float x2;
         float x3;
         float xin1;
         float xin2;
         float xin3;
         float y1;
         float y2;
         float y3;
         xin1 = ZXP(in1);
         xin2 = ZXP(in2); 
         xin3 = ZXP(in3); 
         x1 = g * mem1[(idx-d1)&mask];
         x2 = g * mem2[(idx-d2)&mask];
         x3 = g * mem3[(idx-d3)&mask];
         y1 = l1 * x1 + xin1;
         y2 = a*x1 + l2*x2 + xin2;
         y3 = b*x1 + c*x2 + l3*x3 + xin3;
         mem1[idx&mask] = y1;
         mem2[idx&mask] = y2;  
         mem3[idx&mask] = y3;  
         ZXP(out1) = y1;
         ZXP(out2) = y2;
         ZXP(out3) = y3;
         idx++;
         );
    
    unit->m_idx = idx;
}

//////////////////////////////////////////////////////////////////////////



void JotFDN3LPAP_Ctor(JotFDN3LPAP *unit)
{
    SETCALC(JotFDN3LPAP_next);
    unit->m_idx = 0;
    unit->m_y1 = 0.f;
    unit->m_mask = 4096*8-1;
    memset(unit->m_mem1, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem2, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem3, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem4, 0, sizeof(float) * 4096*8);
    memset(unit->m_mem5, 0, sizeof(float) * 4096*8);
    ZOUT0(0) = 0.0f;
}


void JotFDN3LPAP_next(JotFDN3LPAP *unit, int inNumSamples)
{
    float *in1 = ZIN(0);
    float *in2 = ZIN(1);
    float *in3 = ZIN(2);
    float *out1 = ZOUT(0);
    float *out2 = ZOUT(1);
    float *out3 = ZOUT(2);
    float sr = (float) SAMPLERATE;
    int d1 = (int) (ZIN0(3) * sr);
    int d2 = (int) (ZIN0(4) * sr);
    int d3 = (int) (ZIN0(5) * sr);
    int d4 = (int) (ZIN0(6) * sr); // allpass
    int d5 = (int) (ZIN0(7) * sr); // allpass
    float l1 = ZIN0(8);
    float l2 = ZIN0(9);
    float l3 = ZIN0(10);
    float a = ZIN0(11);
    float b = ZIN0(12);
    float c = ZIN0(13);
    float g = ZIN0(14);
    float apgain = ZIN0(15);
    float apgain2 = ZIN0(16);
    float lp = ZIN0(17);

    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    float *mem4 = unit->m_mem4;
    float *mem5 = unit->m_mem5;
    float lpy1 = unit->m_y1;
    int idx = unit->m_idx;

    LOOP(inNumSamples,
         float x1;
         float x2;
         float x3;
         float x4;
         float x5;
         float xin1;
         float xin2;
         float xin3;
         float y1;
         float y2;
         float y3;
         float ap;
         float ap2;
         xin1 = ZXP(in1);
         xin2 = ZXP(in2); 
         xin3 = ZXP(in3); 
         x1 = g * mem1[(idx-d1)&mask];
         x2 = g * mem2[(idx-d2)&mask];
         x3 = g * mem3[(idx-d3)&mask];
         x4 = mem4[(idx-d4)&mask];
         x5 = mem5[(idx-d5)&mask];
         ap = x3 - apgain * x4;
         ap2 = x1 - apgain2 * x5;
         mem4[idx&mask] = ap;
         mem5[idx&mask] = ap2;
         x1 = x5 + apgain2 * ap2;
         lpy1 = lp*lpy1 + (1.f-lp) * x1;
         x1 = lpy1;
         x3 = x4 + apgain * ap;
         y1 = l1 * x1 + xin1;
         y2 = a*x1 + l2*x2 + xin2;
         y3 = b*x1 + c*x2 + l3*x3 + xin3;
         mem1[idx&mask] = y1;
         mem2[idx&mask] = y2;  
         mem3[idx&mask] = y3;  
         ZXP(out1) = y1;
         ZXP(out2) = y2;
         ZXP(out3) = y3;
         idx++;
         );

    unit->m_y1 = lpy1;
    unit->m_idx = idx;
}

///////////////////////////////////////////////////////////////////////////////


void LBCF_Ctor(LBCF *unit)
{
    int interp = (int) ZIN0(8);
    unit->m_idx = 0;
    float maxdelay = ZIN0(1);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask = delaybufsize-1;     
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    unit->m_lpy1 = 0.f;
    memset(unit->m_apy1, 0, sizeof(float)*16);
    unit->m_firx2 = 0.f;
    unit->m_firx1 = 0.f;
    switch (interp) {
    default:
        if ((ZIN0(9) == 0.f) && (ZIN0(10) == 0.f)) {
            SETCALC(LBCF_next);
        } else {
            SETCALC(LBCF_next_fir);
        }
    }
    ZOUT0(0) = 0.0f;
}

void LBCF_Dtor(LBCF *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
}


void LBCF_next_fir(LBCF *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del = (int) (ZIN0(2) * sr) - 1;
    int del2 = (int) (ZIN0(3) * sr);
    float fb = ZIN0(4);
    float lp = ZIN0(5); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(6);    
    float disp = ZIN0(7);
    float feedforward = ZIN0(8);
    float fir1 = ZIN0(10);
    float fir3 = ZIN0(11);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float firx1 = unit->m_firx1;
    float firx2 = unit->m_firx2;
    int idx = unit->m_idx;
    float fir2 = 1.f-(fir1+fir3);

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 lpy1 = lp1*fbx + lp*lpy1;
                 mem1[idx&mask] = x + fb * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 ZXP(out) = mem1[(idx-del2)&mask];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 lpy1 = lp1*fbx + lp*lpy1;
                 fbx = x + fb * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //fbx = x + fb * lpy1;
                 mem1[idx&mask] = fbx;
                 ZXP(out) = lpy1 + feedforward * fbx;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }
                 lpy1 = lp1*fbx + lp*lpy1; 
                 mem1[idx&mask] = x + fb * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //mem1[idx&mask] = x + fb * lpy1;
                 ZXP(out) = mem1[(idx-del2)&mask];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }
                 lpy1 = lp1*fbx + lp*lpy1;
                 fbx = x + fb * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //fbx = x + fb * lpy1;
                 mem1[idx&mask] = fbx;
                 ZXP(out) = lpy1 + feedforward * fbx;
                 idx++;
                 );
        }
    }

    unit->m_firx1 = firx1;
    unit->m_firx2 = firx2;
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}


void LBCF_next(LBCF *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del = (int) (ZIN0(2) * sr);
    int del2 = (int) (ZIN0(3) * sr);
    float fb = ZIN0(4);
    float lp = ZIN0(5); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(6);    
    float disp = ZIN0(7);
    float feedforward = ZIN0(8);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    int idx = unit->m_idx;

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 lpy1 = lp1*fbx + lp*lpy1; 
                 mem1[idx&mask] = x + fb * lpy1;
                 ZXP(out) = mem1[(idx-del2)&mask];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 lpy1 = lp1*fbx + lp*lpy1;
                 fbx = x + fb * lpy1;
                 mem1[idx&mask] = fbx;
                 ZXP(out) = lpy1 + feedforward * fbx;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }
                 lpy1 = lp1*fbx + lp*lpy1; 
                 mem1[idx&mask] = x + fb * lpy1;
                 ZXP(out) = mem1[(idx-del2)&mask];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 x = ZXP(in);
                 fbx = mem1[(idx-del)&mask];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }
                 lpy1 = lp1*fbx + lp*lpy1;
                 fbx = x + fb * lpy1;
                 mem1[idx&mask] = fbx;
                 ZXP(out) = lpy1 + feedforward * fbx;
                 idx++;
                 );
        }
    }
     
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////


void Velvet_Ctor(Velvet *unit)
{
    SETCALC(Velvet_next);
    unit->m_idx = 0;
    float sr = SAMPLERATE;
    int num_params = (unit->mNumInputs-2)/2;
    int maxdelay = 0;
    int nnegtaps = 0;
    int npostaps = 0;
    for (int k=1; k<num_params; k++) {
        int delay = (int) (ZIN0(k*2) * sr);
        if (delay > maxdelay) {
            maxdelay = delay;
        }
        if (ZIN0(k*2+1) < 0.f) {
            nnegtaps ++;
            unit->m_negtaps[nnegtaps] = delay;
        } else {
            npostaps ++;
            unit->m_postaps[npostaps] = delay;
        }
    }
    unit->m_nnegtaps = nnegtaps;
    unit->m_npostaps = npostaps;
    int delaybufsize = maxdelay;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask = delaybufsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    ZOUT0(0) = 0.0f;
}


void Velvet_Dtor(Velvet *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
}


void Velvet_next(Velvet *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float gain = ZIN0(1);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    int idx = unit->m_idx;
    int *postaps = unit->m_postaps;
    int *negtaps = unit->m_negtaps;
    int npostaps = unit->m_npostaps;
    int nnegtaps = unit->m_nnegtaps;
    int npos4 = 4*(npostaps/4);
    int nneg4 = 4*(nnegtaps/4);

    LOOP(inNumSamples,
         float sum;
         float x;
         int n;
         x = ZXP(in);
         mem1[idx&mask] = x * gain;
         sum = 0.f;
         for (n=0; n<npos4; n+=4) {
             sum += mem1[(idx-postaps[n])&mask];
             sum += mem1[(idx-postaps[n+1])&mask];
             sum += mem1[(idx-postaps[n+2])&mask];
             sum += mem1[(idx-postaps[n+3])&mask];
         }
         for (; n<npostaps; n++) {
             sum += mem1[(idx-postaps[n])&mask];
         }
         for (n=0; n<nneg4; n+=4) {
             sum -= mem1[(idx-negtaps[n])&mask]; 
             sum -= mem1[(idx-negtaps[n+1])&mask];
             sum -= mem1[(idx-negtaps[n+2])&mask];
             sum -= mem1[(idx-negtaps[n+3])&mask];
         }
         for (; n<nnegtaps; n++) {
             sum -= mem1[(idx-negtaps[n])&mask]; 
         }
         ZXP(out) = sum;
         idx++;
         );
    
    unit->m_idx = idx;
}


/////////////////////////////////////////////////////////////


void Velvet1_Ctor(Velvet1 *unit)
{
    SETCALC(Velvet1_next);
    unit->m_idx = 0;
    float sr = SAMPLERATE;
    int num_params = (unit->mNumInputs-2);
    int maxdelay = 0;
    int npostaps = 0;
    for (int k=2; k<num_params; k++) {
        int delay = (int) (ZIN0(k) * sr);
        if (delay > maxdelay) {
            maxdelay = delay;
        }
        unit->m_postaps[npostaps] = delay;
        npostaps++;
    }
    unit->m_npostaps = npostaps;
    int delaybufsize = maxdelay;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask = delaybufsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    ZOUT0(0) = 0.0f;
}


void Velvet1_Dtor(Velvet1 *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
}


void Velvet1_next(Velvet1 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float gain = ZIN0(1);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    int idx = unit->m_idx;
    int *postaps = unit->m_postaps;
    int npostaps = unit->m_npostaps;
    int npos4 = 4*(npostaps/4);
    //fprintf(stderr, "%d, %f\n", npostaps, gain);

    LOOP(inNumSamples,
         float sum;
         float x;
         int n;
         x = ZXP(in);
         mem1[idx&mask] = x * gain;
         sum = 0.f;
         for (n=0; n<npos4; n+=4) {
             sum += mem1[(idx-postaps[n])&mask];
             sum += mem1[(idx-postaps[n+1])&mask];
             sum += mem1[(idx-postaps[n+2])&mask];
             sum += mem1[(idx-postaps[n+3])&mask];
         }
         for (; n<npostaps; n++) {
             sum += mem1[(idx-postaps[n])&mask];
         }
         ZXP(out) = sum;
         idx++;
         );
    
    unit->m_idx = idx;
}


/////////////////////////////////////////////////////////////


void VelvetTaps_Ctor(VelvetTaps *unit)
{
    SETCALC(VelvetTaps_next);
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_idx = 0;
    unit->m_mask = buf->mask;

    float sr = SAMPLERATE;
    int num_params = (unit->mNumInputs-3)/2;
    int maxdelay = 0;
    int nnegtaps = 0;
    int npostaps = 0;
    for (int k=0; k<num_params; k++) {
        int delay = (int) (ZIN0(k*2+3) * sr);
        if (delay > maxdelay) {
            maxdelay = delay;
        }
        if (ZIN0(k*2+4) < 0.f) {
            unit->m_negtaps[nnegtaps] = delay;
            fprintf(stderr, "pos: %d\n", delay);
            nnegtaps ++;
        } else {
            unit->m_postaps[npostaps] = delay;
            fprintf(stderr, "neg: %d\n", delay);
            npostaps ++;
        }
    }
    unit->m_nnegtaps = nnegtaps;
    unit->m_npostaps = npostaps;
    ZOUT0(0) = 0.0f;
}



void VelvetTaps_next(VelvetTaps *unit, int inNumSamples)
{
    float *in = ZIN(1);
    float *out = ZOUT(0);
    float gain = ZIN0(2);
    int mask = unit->m_buf->mask;
    float *mem1 = unit->m_buf->data;
    int idx = unit->m_idx;
    int *postaps = unit->m_postaps;
    int *negtaps = unit->m_negtaps;
    int npostaps = unit->m_npostaps;
    int nnegtaps = unit->m_nnegtaps;
    int npos4 = 4*(npostaps/4);
    int nneg4 = 4*(nnegtaps/4);

    LOOP(inNumSamples,
         float sum;
         float x;
         int n;
         x = ZXP(in);
         mem1[idx&mask] = x * gain;
         sum = 0.f;
         for (n=0; n<npos4; n+=4) {
             sum += mem1[(idx-postaps[n])&mask];
             sum += mem1[(idx-postaps[n+1])&mask];
             sum += mem1[(idx-postaps[n+2])&mask];
             sum += mem1[(idx-postaps[n+3])&mask];
         }
         for (; n<npostaps; n++) {
             sum += mem1[(idx-postaps[n])&mask];
         }
         for (n=0; n<nneg4; n+=4) {
             sum -= mem1[(idx-negtaps[n])&mask]; 
             sum -= mem1[(idx-negtaps[n+1])&mask];
             sum -= mem1[(idx-negtaps[n+2])&mask];
             sum -= mem1[(idx-negtaps[n+3])&mask];
         }
         for (; n<nnegtaps; n++) {
             sum -= mem1[(idx-negtaps[n])&mask]; 
         }
         ZXP(out) = sum;
         idx++;
         );
    
    unit->m_idx = idx;
}


///////////////////////////////////////////////////////////////////////////


void BufLBCF_Ctor(BufLBCF *unit)
{
    int interp = (int) ZIN0(7);
    unit->m_idx = 0;
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_mask = buf->mask;
    unit->m_lpy1 = 0.f;
    memset(unit->m_apy1, 0, sizeof(float)*16);
    switch (interp) {
    default:
        SETCALC(BufLBCF_next);
    }
    ZOUT0(0) = 0.0f;
}


void BufLBCF_next(BufLBCF *unit, int inNumSamples)
{
    float *in = ZIN(1);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del = (int) (ZIN0(2) * sr);
    float fb = ZIN0(3);
    float lp = ZIN0(4); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(5);    
    float disp = ZIN0(6);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask = unit->m_mask;
    float *mem1 = unit->m_buf->data;
    int idx = unit->m_idx;

    if (nfilt==0) {
        LOOP(inNumSamples,
             float x;
             float fbx;
             x = ZXP(in);
             fbx = mem1[(idx-del)&mask];
             lpy1 = lp1*fbx + lp*lpy1;
             fbx = x + fb * lpy1;
             mem1[idx&mask] = fbx;
             ZXP(out) = fbx;
             idx++;
             );
    } else {
        LOOP(inNumSamples,
             float x;
             float fbx;
             x = ZXP(in);
             fbx = mem1[(idx-del)&mask];
             for (int n=0; n<nfilt; ++n) {
                 float tmp = apy1[n];
                 float y0 = fbx - disp*tmp;
                 fbx = disp * y0 + tmp;
                 apy1[n] = y0;
             }
             lpy1 = lp1*fbx + lp*lpy1;
             fbx = x + fb * lpy1;
             mem1[idx&mask] = fbx;
             ZXP(out) = fbx;
             idx++;
             );
    }
     
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////


void BufLAPF_Ctor(BufLAPF *unit)
{
    int interp = (int) ZIN0(7);
    unit->m_idx = 0;
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_mask = buf->mask;
    unit->m_lpy1 = 0.f;
    memset(unit->m_apy1, 0, sizeof(float)*16);
    switch (interp) {
    default:
        SETCALC(BufLAPF_next);
    }
    ZOUT0(0) = 0.0f;
}


void BufLAPF_next(BufLAPF *unit, int inNumSamples)
{
    float *in = ZIN(1);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del = (int) (ZIN0(2) * sr);
    float fb = ZIN0(3);
    float lp = ZIN0(4); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(5);    
    float disp = ZIN0(6);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask = unit->m_mask;
    float *mem1 = unit->m_buf->data;
    int idx = unit->m_idx;

    if (nfilt==0) {
        LOOP(inNumSamples,
             float x;
             float fbx;
             x = ZXP(in);
             fbx = mem1[(idx-del)&mask];
             lpy1 = lp1*fbx + lp*lpy1;
             fbx = x - fb * lpy1;
             mem1[idx&mask] = fb*fbx + lpy1;
             ZXP(out) = fbx;
             idx++;
             );
    } else {
        LOOP(inNumSamples,
             float x;
             float fbx;
             x = ZXP(in);
             fbx = mem1[(idx-del)&mask];
             for (int n=0; n<nfilt; ++n) {
                 float tmp = apy1[n];
                 float y0 = fbx - disp*tmp;
                 fbx = disp * y0 + tmp;
                 apy1[n] = y0;
             }
             lpy1 = lp1*fbx + lp*lpy1;
             fbx = x - fb * lpy1;
             mem1[idx&mask] = fb*fbx+lpy1;
             ZXP(out) = fbx;
             idx++;
             );
    }
     
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}


///////////////////////////////////////////////////////////////////////////////
// tapped allpass chain

void AllpassLatc_Ctor(AllpassLatc *unit)
{
    int nsect = (unit->mNumInputs-1)/2;
    unit->m_idx = 0;
    unit->m_nsect = nsect;

    for (int k=0; k<nsect; k++) {
        unit->m_mem1[k] = (float*) RTAlloc(unit->mWorld, sizeof(float) * 8192);
        memset(unit->m_mem1[k], 0, sizeof(float) * 8192);
        ZOUT0(k) = 0.0f;
    }
    
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    //unit->m_lpy1 = 0.f;
    // memset(unit->m_apy1, 0, sizeof(float)*16);

    SETCALC(AllpassLatc_next);
}


void AllpassLatc_Dtor(AllpassLatc *unit)
{
    int nsect = unit->m_nsect;
    for (int n=0; n<nsect; n++) {
        RTFree(unit->mWorld, unit->m_mem1[n]);
    }
}



void AllpassLatc_next(AllpassLatc *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *outs[16];
    int nsect = unit->m_nsect;
    float sr = (float) SAMPLERATE;
    float k[16];
    int del[16];
    float **mem1 = unit->m_mem1;
    int idx = unit->m_idx;
    int nsect1 = nsect-1;

    //fprintf(stderr, "n: %d\n", nsect);
    for (int n=0; n<nsect; n++) {
        k[n] = ZIN0(n+1+nsect);
        del[n] = (int) (ZIN0(n+1) * sr);
        outs[n] = ZOUT(n);
    }
    outs[nsect] = ZOUT(nsect);
         
    LOOP(inNumSamples, 
         float sum = ZXP(in);
         float stage;
         float mem = mem1[0][(idx-del[0])&8191];
         sum -= k[0] * mem;
         stage = mem + k[0] * sum;
         ZXP(outs[0]) = stage;

         for (int i=1; i<(nsect); i++) {	  
             mem = mem1[i][(idx-del[i])&8191];
             sum -= k[i] * mem;
             stage = mem + k[i] * sum;
             mem1[i-1][idx&8191] = stage;
             ZXP(outs[i]) = stage;
         }

         mem1[nsect1][idx&8191] = sum;
         //stage = mem + k[nsect1] * sum;
         ZXP(outs[nsect]) = stage ;
         idx ++;
         )

        unit->m_idx = idx;
}


///////////////////////////////////////////////////////////////////////////////
// tapped allpass chain

void AllpassLatcLP_Ctor(AllpassLatcLP *unit)
{
    //int interp = (int) ZIN0(7);
    int nsect = (unit->mNumInputs-1)/3;
    unit->m_idx = 0;
    unit->m_nsect = nsect;
    for (int k=0; k<nsect; k++) {
        unit->m_mem1[k] = (float*) RTAlloc(unit->mWorld, sizeof(float) * 8192);
        memset(unit->m_mem1[k], 0, sizeof(float) * 8192);
        unit->m_y1[k] = 0.f;
        ZOUT0(k) = 0.0f;
    }
    SETCALC(AllpassLatcLP_next);
}


void AllpassLatcLP_Dtor(AllpassLatcLP *unit)
{
    int nsect = unit->m_nsect;
    for (int n=0; n<nsect; n++) {
        RTFree(unit->mWorld, unit->m_mem1[n]);
    }
}



void AllpassLatcLP_next(AllpassLatcLP *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *outs[16];
    int nsect = unit->m_nsect;
    float sr = (float) SAMPLERATE;
    float k[16];
    int del[16];
    float y1[16];
    float lp[16];
    float **mem1 = unit->m_mem1;
    int idx = unit->m_idx;
    int nsect1 = nsect-1;
     

    //fprintf(stderr, "n: %d\n", nsect);
    for (int n=0; n<nsect; n++) {
        k[n] = ZIN0(n+1+nsect);
        del[n] = (int) (ZIN0(n+1) * sr);
        lp[n] = ZIN0(n+1+nsect+nsect);
        //fprintf(stderr, "d: %d k: %f\n", del[n], k[n]);
        outs[n] = ZOUT(n);
        y1[n] = unit->m_y1[n];
    }
    outs[nsect] = ZOUT(nsect);
         
    LOOP(inNumSamples, 
         float sum = ZXP(in);
         float stage;
         float mem = mem1[0][(idx-del[0])&8191];
         mem = y1[0] = lp[0]*y1[0] + (1.f-lp[0])*mem;
         sum -= k[0] * mem;
         stage = mem + k[0] * sum;
         ZXP(outs[0]) = stage;
          
         for (int i=1; i<(nsect); i++) {	  
             mem = mem1[i][(idx-del[i])&8191];
             mem = y1[i] = lp[i]*y1[i] + (1.f-lp[i])*mem;
             sum -= k[i] * mem;
             stage = mem + k[i] * sum;
             mem1[i-1][idx&8191] = stage;
             ZXP(outs[i]) = stage;
         }
         mem1[nsect1][idx&8191] = sum;
         //stage = mem + k[nsect1] * sum;
         ZXP(outs[nsect]) = stage ;
         idx ++;
         )

    unit->m_idx = idx;

    for (int n=0; n<nsect; n++) {
        unit->m_y1[n] = y1[n];
    }
}

/////////////////////////////////////////////////////////////



void BufOneTap_Ctor(BufOneTap *unit)
{
    SETCALC(BufOneTap_next);
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_idx = 0;
    unit->m_mask = buf->mask;
    //memset(unit->m_mem, 0, sizeof(float) * 4096*8);
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    int numOutputs = unit->mNumOutputs;
    for (uint32 channel=0; channel<numOutputs; ++channel) {
        ZOUT0(channel) = 0.0f;
    }
}


void BufOneTap_next(BufOneTap *unit, int inNumSamples)
{
    float sr = (float) SAMPLERATE;
    int del = (int) (ZIN0(1) * sr);
    float gain = ZIN0(2);
    int mask = unit->m_mask;
    float *mem = unit->m_buf->data;
    int idx = unit->m_idx;
    int numOutputs = unit->mNumOutputs;
    int idxdel = idx-del;
    float *out[numOutputs];

    //fprintf(stderr, "mask: %d nout: %d\n", mask, numOutputs);
     
    for (uint32 channel=0; channel<numOutputs; ++channel) {
        out[channel] = ZOUT(channel);
    }
     
    LOOP(inNumSamples,
         //ZXP(out) = gain * mem[(idx-del)&mask];
         const float* table1 = mem + ((idxdel*numOutputs)&mask);
         for (uint32 channel=0; channel<numOutputs; ++channel) {
             ZXP(out[channel]) = gain * *table1;
             table1++;
         }
         idxdel++;
         );
    

    unit->m_idx = idx + inNumSamples;
}

/////////////////////////////////////////////////////////////



void BufOneTapC_Ctor(BufOneTapC *unit)
{
    SETCALC(BufOneTapC_next);
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_idx = 0;
    unit->m_lastdel = (ZIN0(1) * SAMPLERATE);
    unit->m_mask = buf->mask;
    //memset(unit->m_mem, 0, sizeof(float) * 4096*8);
    //int num_params = unit->mNumInputs-1;
    //int num_taps = num_params/2;
    int numOutputs = unit->mNumOutputs;
    for (uint32 channel=0; channel<numOutputs; ++channel) {
        ZOUT0(channel) = 0.0f;
    }
}


void BufOneTapC_next(BufOneTapC *unit, int inNumSamples)
{
    float sr = (float) SAMPLERATE;
    float newdel = (ZIN0(1) * sr);
    float del = unit->m_lastdel;
    float delinc = (newdel-del)/inNumSamples;
    float gain = ZIN0(2);
    int mask = unit->m_mask;
    float *mem = unit->m_buf->data;
    int idx = unit->m_idx;
    int numOutputs = unit->mNumOutputs;
    float *out[numOutputs];
    unit->m_lastdel = newdel;
    //fprintf(stderr, "mask: %d nout: %d\n", mask, numOutputs);
     
    for (uint32 channel=0; channel<numOutputs; ++channel) {
        out[channel] = ZOUT(channel);
    }
     
    LOOP(inNumSamples,
         int idel = (int)del;
         float frac = del-(float)idel;
         int idxdel = idx-idel;
         del += delinc;
         for (uint32 channel=0; channel<numOutputs; ++channel) {
             float x1 = mem[((idxdel-2)*numOutputs+channel)&mask];
             float x2 = mem[((idxdel-1)*numOutputs+channel)&mask];
             float x3 = mem[(idxdel*numOutputs+channel)&mask];
             float x4 = mem[((idxdel+1)*numOutputs+channel)&mask];
             ZXP(out[channel]) = gain * cubicinterp(frac, x1, x2, x3, x4);              
         }
         idx++;
         );
    
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void BufWrite_Ctor(BufWrite *unit)
{
    SETCALC(BufWrite_next);
    CTOR_GET_BUF
        unit->m_buf = buf;
    unit->m_idx = 0;
    unit->m_mask = buf->mask;
}


void BufWrite_next(BufWrite *unit, int inNumSamples)
{
    int numChannels = unit->mNumInputs-1;
    int mask = unit->m_mask;
    float *mem = unit->m_buf->data;
    int idx = unit->m_idx;
    float *in[numChannels];

    //fprintf(stderr, "mask: %d nout: %d\n", mask, numOutputs);
    for (uint32 channel=0; channel<numChannels; ++channel) {
        in[channel] = ZIN(channel+1);
    }
     
    LOOP(inNumSamples,
         float* table0 = mem + ((idx*numChannels)&mask);
         for (uint32 channel=0; channel<numChannels; ++channel) {
             table0[channel] = ZXP(in[channel]);
         }
         idx++;
         );
    
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////////////////////////////////


void Hadamard8Delay_Ctor(Hadamard8Delay *unit)
{
    SETCALC(Hadamard8Delay_next);
    int memsize = (int) (ZIN0(8) * SAMPLERATE);    
    unit->m_idx = 0;
    memsize = NEXTPOWEROFTWO(memsize);
    unit->m_mask = memsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem1, 0, sizeof(float) * memsize);
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem2, 0, sizeof(float) * memsize);
    unit->m_mem3 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem3, 0, sizeof(float) * memsize);
    unit->m_mem4 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem4, 0, sizeof(float) * memsize);
    unit->m_mem5 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem5, 0, sizeof(float) * memsize);
    unit->m_mem6 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem6, 0, sizeof(float) * memsize);
    unit->m_mem7 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem7, 0, sizeof(float) * memsize);
    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
    ZOUT0(2) = 0.f;
    ZOUT0(3) = 0.f;
    ZOUT0(4) = 0.f;
    ZOUT0(5) = 0.f;
    ZOUT0(6) = 0.f;
    ZOUT0(7) = 0.f;
}


void Hadamard8Delay_Dtor(Hadamard8Delay *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
    RTFree(unit->mWorld, unit->m_mem3);
    RTFree(unit->mWorld, unit->m_mem4);
    RTFree(unit->mWorld, unit->m_mem5);
    RTFree(unit->mWorld, unit->m_mem6);
    RTFree(unit->mWorld, unit->m_mem7);
}


void Hadamard8Delay_next(Hadamard8Delay *unit, int inNumSamples)
{
    float *out0 = ZOUT(0);
    float *out1 = ZOUT(1);
    float *out2 = ZOUT(2);
    float *out3 = ZOUT(3);
    float *out4 = ZOUT(4);
    float *out5 = ZOUT(5);
    float *out6 = ZOUT(6);
    float *out7 = ZOUT(7);
    float *in0 = ZIN(0);
    float *in1 = ZIN(1);
    float *in2 = ZIN(2);
    float *in3 = ZIN(3);
    float *in4 = ZIN(4);
    float *in5 = ZIN(5);
    float *in6 = ZIN(6);
    float *in7 = ZIN(7);
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    float *mem4 = unit->m_mem4;
    float *mem5 = unit->m_mem5;
    float *mem6 = unit->m_mem6;
    float *mem7 = unit->m_mem7;
    float sr = SAMPLERATE;
    int mask = unit->m_mask;
    int idx = unit->m_idx;
    int del1 = ZIN0(9) * sr;
    int del2 = ZIN0(10) * sr;
    int del3 = ZIN0(11) * sr;
    int del4 = ZIN0(12) * sr;
    int del5 = ZIN0(13) * sr;
    int del6 = ZIN0(14) * sr;
    int del7 = ZIN0(15) * sr;
    float gain = ZIN0(16) * 0.35355339059327f;

    LOOP(inNumSamples, 
         float t;
         float x0;
         float x1;
         float x2;
         float x3;
         float x4;
         float x5;
         float x6;
         float x7;
         x0 = ZXP(in0) * gain;
         x1 = mem1[(idx-del1)&mask] * gain;
         mem1[idx&mask] = ZXP(in1);
         x2 = mem2[(idx-del2)&mask] * gain;
         mem2[idx&mask] = ZXP(in2);
         x3 = mem3[(idx-del3)&mask] * gain;
         mem3[idx&mask] = ZXP(in3);
         x4 = mem4[(idx-del4)&mask] * gain;
         mem4[idx&mask] = ZXP(in4);
         x5 = mem5[(idx-del5)&mask] * gain;
         mem5[idx&mask] = ZXP(in5);
         x6 = mem6[(idx-del6)&mask] * gain;
         mem6[idx&mask] = ZXP(in6);
         x7 = mem7[(idx-del7)&mask] * gain;
         mem7[idx&mask] = ZXP(in7);

         t = x0 - x1; x0 += x1;  x1 = t;
         t = x2 - x3; x2 += x3;  x3 = t;
         t = x4 - x5; x4 += x5;  x5 = t;
         t = x6 - x7; x6 += x7;  x7 = t;
         t = x0 - x2; x0 += x2;  x2 = t;
         t = x1 - x3; x1 += x3;  x3 = t;
         t = x4 - x6; x4 += x6;  x6 = t;
         t = x5 - x7; x5 += x7;  x7 = t;
         t = x0 - x4; x0 += x4;  x4 = t;
         t = x1 - x5; x1 += x5;  x5 = t;
         t = x2 - x6; x2 += x6;  x6 = t;
         t = x3 - x7; x3 += x7;  x7 = t;

         ZXP(out0) = x0;
         ZXP(out1) = x1;
         ZXP(out2) = x2;
         ZXP(out3) = x3;
         ZXP(out4) = x4;
         ZXP(out5) = x5;
         ZXP(out6) = x6;
         ZXP(out7) = x7;
         idx++;
         )

        unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////////////////////



void Hadamard4Delay_Ctor(Hadamard4Delay *unit)
{
    SETCALC(Hadamard4Delay_next);
    int memsize = (int) (ZIN0(4) * SAMPLERATE);    
    unit->m_idx = 0;
    memsize = NEXTPOWEROFTWO(memsize);
    unit->m_mask = memsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem1, 0, sizeof(float) * memsize);
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem2, 0, sizeof(float) * memsize);
    unit->m_mem3 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem3, 0, sizeof(float) * memsize);
    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
    ZOUT0(2) = 0.f;
    ZOUT0(3) = 0.f;
}


void Hadamard4Delay_Dtor(Hadamard4Delay *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
    RTFree(unit->mWorld, unit->m_mem3);
}


void Hadamard4Delay_next(Hadamard4Delay *unit, int inNumSamples)
{
    float *out0 = ZOUT(0);
    float *out1 = ZOUT(1);
    float *out2 = ZOUT(2);
    float *out3 = ZOUT(3);
    float *in0 = ZIN(0);
    float *in1 = ZIN(1);
    float *in2 = ZIN(2);
    float *in3 = ZIN(3);
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    float sr = SAMPLERATE;
    int mask = unit->m_mask;
    int idx = unit->m_idx;
    int del1 = ZIN0(5) * sr;
    int del2 = ZIN0(6) * sr;
    int del3 = ZIN0(7) * sr;
    float gain = ZIN0(8) * 0.5f;

    LOOP(inNumSamples, 
         float t0;
         float t1;
         float t2;
         float t3;
         float x0;
         float x1;
         float x2;
         float x3;
         x0 = ZXP(in0) * gain;
         x1 = mem1[(idx-del1)&mask] * gain;
         mem1[idx&mask] = ZXP(in1);
         x2 = mem2[(idx-del2)&mask] * gain;
         mem2[idx&mask] = ZXP(in2);
         x3 = mem3[(idx-del3)&mask] * gain;
         mem3[idx&mask] = ZXP(in3);

         t0 = x0+x1;
         t1 = x2+x3;
         t2 = x0-x1;
         t3 = x2-x3;
         x0 = t0+t1;  // x0+x1+x2+x3
         x1 = t0-t1;  // x0+x1-x2-x3
         x2 = t2+t3;  // x0-x1+x2-x3
         x3 = t2-t3; 

         ZXP(out0) = x0;
         ZXP(out1) = x1;
         ZXP(out2) = x2;
         ZXP(out3) = x3;
         idx++;
         )

        unit->m_idx = idx;
}


/////////////////////////////////////////////////////////////////////////////



void Ortho3Delay_Ctor(Ortho3Delay *unit)
{
    SETCALC(Ortho3Delay_next);
    int memsize = (int) (ZIN0(3) * SAMPLERATE);    
    unit->m_idx = 0;
    memsize = NEXTPOWEROFTWO(memsize);
    unit->m_mask = memsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem1, 0, sizeof(float) * memsize);
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem2, 0, sizeof(float) * memsize);

    float a=ZIN0(7);
    float b=ZIN0(8);
    float c=ZIN0(9);

    unit->m_m0 = cos(a)*cos(b)-cos(c)*sin(a)*sin(b);
    unit->m_m1 = -cos(a)*sin(b)-cos(c)*sin(a)*cos(b);
    unit->m_m2 = sin(a)*sin(c);
    unit->m_m3 = sin(a)*cos(b)+cos(c)*cos(a)*sin(b);
    unit->m_m4 = -sin(a)*sin(b)+cos(c)*cos(a)*cos(b);
    unit->m_m5 = -cos(a)*sin(c);
    unit->m_m6 = sin(b)*sin(c);
    unit->m_m7 = cos(b)*sin(c);
    unit->m_m8 = cos(c);

    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
    ZOUT0(2) = 0.f;
}


void Ortho3Delay_Dtor(Ortho3Delay *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
}


void Ortho3Delay_next(Ortho3Delay *unit, int inNumSamples)
{
    float *out0 = ZOUT(0);
    float *out1 = ZOUT(1);
    float *out2 = ZOUT(2);
    float *in0 = ZIN(0);
    float *in1 = ZIN(1);
    float *in2 = ZIN(2);
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float sr = SAMPLERATE;
    int mask = unit->m_mask;
    int idx = unit->m_idx;
    int del1 = ZIN0(4) * sr;
    int del2 = ZIN0(5) * sr;
    float gain = ZIN0(6) * 0.5f;
    float m0 = unit->m_m0;
    float m1 = unit->m_m1;
    float m2 = unit->m_m2;
    float m3 = unit->m_m3;
    float m4 = unit->m_m4;
    float m5 = unit->m_m5;
    float m6 = unit->m_m6;
    float m7 = unit->m_m7;
    float m8 = unit->m_m8;

    LOOP(inNumSamples, 
         float t0;
         float t1;
         float t2;
         float x0;
         float x1;
         float x2;
         x0 = ZXP(in0) * gain;
         x1 = mem1[(idx-del1)&mask] * gain;
         mem1[idx&mask] = ZXP(in1);
         x2 = mem2[(idx-del2)&mask] * gain;
         mem2[idx&mask] = ZXP(in2);

         //t0 = 0.7071067811865477*(x1-x0);
         //t1 = -0.6532814824381881*(x0+x1)+0.3826834323650898*x2;
         //t2 = 0.2705980500730985*(x0+x1)+0.9238795325112867*x2;
         t0 = m0*x0 + m1*x1 + m2*x2;
         t1 = m3*x0 + m4*x1 + m5*x2;
         t2 = m6*x0 + m7*x1 + m8*x2;

         ZXP(out0) = t0;
         ZXP(out1) = t1;
         ZXP(out2) = t2;
         idx++;
         )

        unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////////////////////



void Nested2AllpassN_Ctor(Nested2AllpassN *unit)
{
    SETCALC(Nested2AllpassN_next);
    unit->m_idx = 0;
    int memsize = (int) (ZIN0(1) * SAMPLERATE);    
    memsize = NEXTPOWEROFTWO(memsize);
    //fprintf(stderr, "%d\n", memsize); 
    unit->m_mask1 = memsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem1, 0, sizeof(float) * memsize);

    memsize = (int) (ZIN0(4) * SAMPLERATE);    
    memsize = NEXTPOWEROFTWO(memsize);
    //fprintf(stderr, "%d\n", memsize); 
    unit->m_mask2 = memsize-1;
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem2, 0, sizeof(float) * memsize);

    memsize = (int) (ZIN0(7) * SAMPLERATE);    
    memsize = NEXTPOWEROFTWO(memsize);
    //fprintf(stderr, "%d\n", memsize); 
    unit->m_mask3 = memsize-1;
    unit->m_mem3 = (float*) RTAlloc(unit->mWorld, sizeof(float) * memsize);
    memset(unit->m_mem3, 0, sizeof(float) * memsize);

    unit->m_y1 = 0.f;
    ZOUT0(0) = 0.f;
}


void Nested2AllpassN_Dtor(Nested2AllpassN *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
    RTFree(unit->mWorld, unit->m_mem3);
}


void Nested2AllpassN_next(Nested2AllpassN *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float fb1 = ZIN0(3);
    float fb2 = ZIN0(6);
    float fb3 = ZIN0(9);
    float lp = ZIN0(10);
    float lp1 = (1.f-lp);
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    float y1 = unit->m_y1;
    float sr = SAMPLERATE;
    int mask1 = unit->m_mask1;
    int mask2 = unit->m_mask2;
    int mask3 = unit->m_mask3;
    int idx = unit->m_idx;
    int del1 = ZIN0(2) * sr;
    int del2 = ZIN0(5) * sr;
    int del3 = ZIN0(8) * sr;

    //fprintf(stderr, "%d %f %d %f %d %f\n", del1, fb1, del2, fb2, del3, fb3); 

    LOOP(inNumSamples, 
         float x;
         float mem;
         float tmp;
         float tmp2;
         tmp2 = mem3[(idx-del3)&mask3];

         mem = mem1[(idx-del1)&mask1];
         tmp = tmp2 - fb1 * mem;
         mem1[idx&mask1] = tmp;
         tmp2 = mem + fb1 * tmp; 
         mem = mem2[(idx-del2)&mask2];
         tmp = tmp2 - fb2 * mem;
         mem2[idx&mask2] = tmp;
         tmp2 = mem + fb2 * tmp; 

         x = ZXP(in);
         tmp = x - fb3 * tmp2;
         y1 = lp*y1 + lp1*tmp;
         mem3[idx&mask3] = y1;

         ZXP(out) = fb3 * tmp + tmp2;
         idx++;
         )

        unit->m_y1 = y1;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////


void NestedLBCF_Ctor(NestedLBCF *unit)
{
    int interp = (int) ZIN0(12);
    unit->m_idx = 0;
    float maxdelay = ZIN0(1);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask1 = delaybufsize-1;
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);

    maxdelay = ZIN0(4);
    delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask2 = delaybufsize-1;
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem2, 0, sizeof(float) * delaybufsize);

    unit->m_lpy1 = 0.f;
    memset(unit->m_apy1, 0, sizeof(float)*16);
    unit->m_firx1 = 0.f;
    unit->m_firx2 = 0.f;
    switch (interp) {
    default:
        if ((ZIN0(13)==0.f) && (ZIN0(14)==0.f)) {
            SETCALC(NestedLBCF_next);
        } else {
            SETCALC(NestedLBCF_next_fir);
        }
    }
    ZOUT0(0) = 0.0f;
}

void NestedLBCF_Dtor(NestedLBCF *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
}


void NestedLBCF_next(NestedLBCF *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del1 = (int) (ZIN0(2) * sr);
    int del2 = (int) (ZIN0(5) * sr);
    int del3 = (int) (ZIN0(7) * sr);
    float fb1 = ZIN0(3);
    float fb2 = ZIN0(6);
    float lp = ZIN0(8); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(9);    
    float disp = ZIN0(10);
    float feedforward = ZIN0(11);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask1 = unit->m_mask1;
    int mask2 = unit->m_mask2;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    }
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}


void NestedLBCF_next_fir(NestedLBCF *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del1 = (int) (ZIN0(2) * sr)-1;
    int del2 = (int) (ZIN0(5) * sr);
    int del3 = (int) (ZIN0(7) * sr);
    float fb1 = ZIN0(3);
    float fb2 = ZIN0(6);
    float lp = ZIN0(8); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(9);    
    float disp = ZIN0(10);
    float feedforward = ZIN0(11);
    float fir1 = ZIN0(13);
    float fir3 = ZIN0(14);
    float firx1 = unit->m_firx1;
    float firx2 = unit->m_firx2;
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask1 = unit->m_mask1;
    int mask2 = unit->m_mask2;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;
    float fir2 = 1.f-(fir1+fir3);

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;

                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    }

    unit->m_firx1 = firx1;
    unit->m_firx2 = firx2;          
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////


void BufNestedLBCF_Ctor(BufNestedLBCF *unit)
{
    int interp = (int) ZIN0(12);
    unit->m_idx = 0;
    CTOR_GET_BUF
        unit->m_buf = buf;
    int maxdelay = ZIN0(4);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask2 = delaybufsize-1;
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem2, 0, sizeof(float) * delaybufsize);

    unit->m_lpy1 = 0.f;
    memset(unit->m_apy1, 0, sizeof(float)*16);
    unit->m_firx1 = 0.f;
    unit->m_firx2 = 0.f;
    switch (interp) {
    default:
        if ((ZIN0(13)==0.f) && (ZIN0(14)==0.f)) {
            SETCALC(BufNestedLBCF_next);
        } else {
            SETCALC(BufNestedLBCF_next_fir);
        }
    }
    ZOUT0(0) = 0.0f;
}



void BufNestedLBCF_next(BufNestedLBCF *unit, int inNumSamples)
{
    float *in = ZIN(1);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del1 = (int) (ZIN0(2) * sr);
    int del2 = (int) (ZIN0(5) * sr);
    int del3 = (int) (ZIN0(7) * sr);
    float fb1 = ZIN0(3);
    float fb2 = ZIN0(6);
    float lp = ZIN0(8); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(9);    
    float disp = ZIN0(10);
    float feedforward = ZIN0(11);
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask1 = unit->m_buf->mask;
    int mask2 = unit->m_mask2;
    float *mem1 = unit->m_buf->data;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    }
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}


void BufNestedLBCF_next_fir(BufNestedLBCF *unit, int inNumSamples)
{
    float *in = ZIN(1);
    float *out = ZOUT(0);
    float sr = (float) SAMPLERATE;
    int del1 = (int) (ZIN0(2) * sr)-1;
    int del2 = (int) (ZIN0(5) * sr);
    int del3 = (int) (ZIN0(7) * sr);
    float fb1 = ZIN0(3);
    float fb2 = ZIN0(6);
    float lp = ZIN0(8); 
    float lp1 = (1.f-lp); 
    int nfilt = (int) ZIN0(9);    
    float disp = ZIN0(10);
    float feedforward = ZIN0(11);
    float fir1 = ZIN0(13);
    float fir3 = ZIN0(14);
    float firx1 = unit->m_firx1;
    float firx2 = unit->m_firx2;
    float lpy1 = unit->m_lpy1;
    float *apy1 = unit->m_apy1;
    int mask1 = unit->m_buf->mask;
    int mask2 = unit->m_mask2;
    float *mem1 = unit->m_buf->data;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;
    float fir2 = 1.f-(fir1+fir3);

    if (nfilt==0) {
        if (feedforward==0.f) {
            // feedforward comes without additional delay
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    } else {
        if (feedforward==0.f) {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 mem1[idx&mask1] = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //mem1[idx&mask1] = x + fb1*lpy1;
                 ZXP(out) = mem1[(idx-del3)&mask1];
                 idx++;
                 );
        } else {
            LOOP(inNumSamples,
                 float x;
                 float fbx;
                 float fbx2;
                 float tmp;
                 fbx = mem1[(idx-del1)&mask1];
                 for (int n=0; n<nfilt; ++n) {
                     float tmp = apy1[n];
                     float y0 = fbx - disp*tmp;
                     fbx = disp * y0 + tmp;
                     apy1[n] = y0;
                 }                  
                 // inner filter
                 fbx2 = mem2[(idx-del2)&mask2];
                 tmp = fbx - fb2*fbx2;
                 mem2[idx&mask2] = tmp;
                 tmp = fb2*tmp + fbx2;
                 lpy1 = lp1*tmp + lp*lpy1; 

                 x = ZXP(in);
                 tmp = x + fb1 * (fir1*firx1 + fir2*firx2 + fir3*lpy1);
                 firx2 = lpy1;
                 firx1 = firx2;
                 //tmp = x + fb1*lpy1;
                 mem1[idx&mask1] = tmp;
                 ZXP(out) = lpy1 + feedforward * tmp;
                 idx++;
                 );
        }
    }

    unit->m_firx1 = firx1;
    unit->m_firx2 = firx2;          
    unit->m_lpy1 = lpy1;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// tapped allpass chain

void Allpass2N_Ctor(Allpass2N *unit)
{
    unit->m_idx = 0;

    float maxdelay = ZIN0(1);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask = delaybufsize-1;     
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem2, 0, sizeof(float) * delaybufsize);
    unit->m_freq = -10022;
    unit->m_rq = -10022;
    ZOUT0(0) = 0.0f;
    SETCALC(Allpass2N_next);
}


void Allpass2N_Dtor(Allpass2N *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
}


void Allpass2N_next(Allpass2N *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    int del = (int) (ZIN0(2) * SAMPLERATE);
    float nextfreq = ZIN0(3);
    float nextrq = ZIN0(4);
    int mask = unit->m_mask;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    int idx = unit->m_idx;
    float a0, a1, a2, b1, b2;

    if ((unit->m_freq != nextfreq) || (unit->m_rq != nextrq)) {
        float w0 = twopi * (double)nextfreq * SAMPLEDUR;
        float alpha = sin(w0) * (double)nextrq * 0.5;
        float b0rz = 1. / (1. + alpha);
        b1 = 2. * b0rz * cos(w0);
        a0 = (1. - alpha) * b0rz;
        a1 = -b1;
        a2 = 1.;
        b2 = -a0;
        //fprintf(stderr, "b1 %f a0 %f\n", b1, a0);
        unit->m_b1 = b1;
        unit->m_a0 = a0;
        unit->m_freq = nextfreq;
        unit->m_rq = nextrq;
    } else {
        b1 = unit->m_b1;
        a0 = unit->m_a0;
        a1 = -b1;
        a2 = 1.;
        b2 = -a0;
    }
         
    LOOP(inNumSamples, 
         float x = ZXP(in);
         float y0;
         float y1;
         float y2;

         y1 = mem1[(idx-del) & mask];
         y2 = mem2[(idx-del) & mask];

         y0 = x + b1 * y1 + b2 * y2;
         ZXP(out) = a0 * y0 + a1 * y1 + a2 * y2;

         mem2[idx&mask] = y1;
         mem1[idx&mask] = y0;
         idx++;
         )

        unit->m_idx = idx;
}

///////////////////////////////////////////////////////////////////////////////
// tapped allpass chain

void NestedAllpass2N_Ctor(NestedAllpass2N *unit)
{
    unit->m_idx = 0;

    float maxdelay = ZIN0(4);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two

    unit->m_mask = delaybufsize-1;     
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    unit->m_mem2 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem2, 0, sizeof(float) * delaybufsize);

    maxdelay = ZIN0(1);
    delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    unit->m_mask3 = delaybufsize-1;
    unit->m_mem3 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem3, 0, sizeof(float) * delaybufsize);
    unit->m_freq = -10022;
    unit->m_rq = -10022;
    ZOUT0(0) = 0.0f;
    SETCALC(NestedAllpass2N_next);
}


void NestedAllpass2N_Dtor(NestedAllpass2N *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
    RTFree(unit->mWorld, unit->m_mem2);
    RTFree(unit->mWorld, unit->m_mem3);
}


void NestedAllpass2N_next(NestedAllpass2N *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    int del = (int) (ZIN0(5) * SAMPLERATE);
    int del3 = (int) (ZIN0(2) * SAMPLERATE);
    float coef = ZIN0(3);
    float nextfreq = ZIN0(6);
    float nextrq = ZIN0(7);
    int mask = unit->m_mask;
    int mask3 = unit->m_mask3;
    float *mem1 = unit->m_mem1;
    float *mem2 = unit->m_mem2;
    float *mem3 = unit->m_mem3;
    int idx = unit->m_idx;
    float a0, a1, a2, b1, b2;

    if ((unit->m_freq != nextfreq) || (unit->m_rq != nextrq)) {
        float w0 = twopi * (double)nextfreq * SAMPLEDUR;
        float alpha = sin(w0) * (double)nextrq * 0.5;
        float b0rz = 1. / (1. + alpha);
        b1 = 2. * b0rz * cos(w0);
        a0 = (1. - alpha) * b0rz;
        a1 = -b1;
        a2 = 1.;
        b2 = -a0;

        unit->m_b1 = b1;
        unit->m_a0 = a0;
        unit->m_freq = nextfreq;
        unit->m_rq = nextrq;
    } else {
        b1 = unit->m_b1;
        a0 = unit->m_a0;
        a1 = -b1;
        a2 = 1.;
        b2 = -a0;
    }
         
    LOOP(inNumSamples, 
         float xin = ZXP(in);
         float x;
         float y0;
         float y1;
         float y2;
         float del1 = mem3[(idx-del3)&mask3];

         y1 = mem1[(idx-del) & mask];
         y2 = mem2[(idx-del) & mask];

         y0 = del1 + b1 * y1 + b2 * y2;
         x = (a0 * y0 + a1 * y1 + a2 * y2);

         mem2[idx&mask] = y1;
         mem1[idx&mask] = y0;

         ZXP(out) = coef * xin + x;

         mem3[idx&mask3] = -coef * x + xin;

         idx++;
         )

        unit->m_idx = idx;
}


///////////////////////////////////////////////////////////////////////////////


void FbTapper_Ctor(FbTapper *unit)
{
    unit->m_idx = 0;
    float maxdelay = ZIN0(1);
    int delaybufsize = maxdelay * SAMPLERATE;
    delaybufsize = delaybufsize + BUFLENGTH;
    delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
    //unit->m_delaybufsize = delaybufsize;
    unit->m_mask = delaybufsize-1;     
    unit->m_mem1 = (float*) RTAlloc(unit->mWorld, sizeof(float) * delaybufsize);
    memset(unit->m_mem1, 0, sizeof(float) * delaybufsize);
    unit->m_lpy1 = 0.f;
    SETCALC(FbTapper_next);
    ZOUT0(0) = 0.0f;
}

void FbTapper_Dtor(FbTapper *unit)
{
    RTFree(unit->mWorld, unit->m_mem1);
}


void FbTapper_next(FbTapper *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float fb = ZIN0(2);
    float z = ZIN0(3);
    float b1 = ZIN0(4);
    float b1abs = fabsf(b1);
    float zabs = fabsf(z);
    float a0 = (1.f-b1abs) * (1.f-zabs) * fb;
    float a1 = (1.f-b1abs) * zabs * fb;
    float lpy = unit->m_lpy1;
    float sr = (float) SAMPLERATE;
    float *mem1 = unit->m_mem1;
    int mask = unit->m_mask;
    int idx0 = unit->m_idx;
    int fbdels[4];
    float fbs[4];
    int outdels[8];
    float *outputs[8];
    float outamps[8];
    int idx;

    int numintaps = (int) ZIN0(5);
    int numouttaps = (int) ZIN0(6);

    for (int k=0; k<numintaps; k++) {
        fbdels[k] = (int) (ZIN0(7+k) * sr);
        fbs[k] = ZIN0(7+k+numintaps);
    }

    for (int k=0; k<numouttaps; k++) {
        outdels[k] = (int) (ZIN0(7+k+numintaps*2) * sr);
        outamps[k] = ZIN0(7+k+numintaps*2+numouttaps);
        outputs[k] = ZOUT(k);
    }

    idx = idx0;
    switch (numintaps) {
    case 1:
        LOOP(inNumSamples,
             float x;
             float tmp;
             float fbx;
             // read feedback taps
             fbx = fbs[0] * mem1[(idx-fbdels[0])&mask];
             // filter the feedback signal
             tmp = fbx + b1*lpy;
             fbx = a0 * tmp + a1 * lpy;
             lpy = tmp;
             // add the input signal
             fbx += ZXP(in);
             mem1[idx&mask] = fbx;
             idx++;
             );
        break;
    case 2:
        LOOP(inNumSamples,
             float x;
             float tmp;
             float fbx;
             // read feedback taps
             fbx = fbs[0] * mem1[(idx-fbdels[0])&mask];
             fbx += fbs[1] * mem1[(idx-fbdels[1])&mask];
             // filter the feedback signal
             tmp = fbx + b1*lpy;
             fbx = a0 * tmp + a1 * lpy;
             lpy = tmp;
             // add the input signal
             fbx += ZXP(in);
             mem1[idx&mask] = fbx;
             idx++;
             );
        break;
    case 3:
        LOOP(inNumSamples,
             float x;
             float tmp;
             float fbx;
             // read feedback taps
             fbx = fbs[0] * mem1[(idx-fbdels[0])&mask];
             fbx += fbs[1] * mem1[(idx-fbdels[1])&mask];
             fbx += fbs[2] * mem1[(idx-fbdels[2])&mask];
             // filter the feedback signal
             tmp = fbx + b1*lpy;
             fbx = a0 * tmp + a1 * lpy;
             lpy = tmp;
             // add the input signal
             fbx += ZXP(in);
             mem1[idx&mask] = fbx;
             idx++;
             );
        break;
    case 4:
        LOOP(inNumSamples,
             float x;
             float tmp;
             float fbx;
             // read feedback taps
             fbx = fbs[0] * mem1[(idx-fbdels[0])&mask];
             fbx += fbs[1] * mem1[(idx-fbdels[1])&mask];
             fbx += fbs[2] * mem1[(idx-fbdels[2])&mask];
             fbx += fbs[3] * mem1[(idx-fbdels[3])&mask];
             // filter the feedback signal
             tmp = fbx + b1*lpy;
             fbx = a0 * tmp + a1 * lpy;
             lpy = tmp;
             // add the input signal
             fbx += ZXP(in);
             mem1[idx&mask] = fbx;
             idx++;
             );
        break;
    }

    idx = idx0;
    switch (numouttaps) {
    case 1:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             idx++;
             );
        break;
    case 2:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             ZXP(outputs[1]) = outamps[1] * mem1[(idx-outdels[1])&mask];
             idx++;
             );
        break;
    case 3:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             ZXP(outputs[1]) = outamps[1] * mem1[(idx-outdels[1])&mask];
             ZXP(outputs[2]) = outamps[2] * mem1[(idx-outdels[2])&mask];
             idx++;
             );
        break;
    case 4:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             ZXP(outputs[1]) = outamps[1] * mem1[(idx-outdels[1])&mask];
             ZXP(outputs[2]) = outamps[2] * mem1[(idx-outdels[2])&mask];
             ZXP(outputs[3]) = outamps[3] * mem1[(idx-outdels[3])&mask];
             idx++;
             );
        break;
    case 5:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             ZXP(outputs[1]) = outamps[1] * mem1[(idx-outdels[1])&mask];
             ZXP(outputs[2]) = outamps[2] * mem1[(idx-outdels[2])&mask];
             ZXP(outputs[3]) = outamps[3] * mem1[(idx-outdels[3])&mask];
             ZXP(outputs[4]) = outamps[4] * mem1[(idx-outdels[4])&mask];
             idx++;
             );
        break;
    case 6:
        LOOP(inNumSamples,
             // tap for outputs
             ZXP(outputs[0]) = outamps[0] * mem1[(idx-outdels[0])&mask];
             ZXP(outputs[1]) = outamps[1] * mem1[(idx-outdels[1])&mask];
             ZXP(outputs[2]) = outamps[2] * mem1[(idx-outdels[2])&mask];
             ZXP(outputs[3]) = outamps[3] * mem1[(idx-outdels[3])&mask];
             ZXP(outputs[4]) = outamps[4] * mem1[(idx-outdels[4])&mask];
             ZXP(outputs[5]) = outamps[5] * mem1[(idx-outdels[5])&mask];
             idx++;
             );
        break;
    }
     
    unit->m_lpy1 = lpy;
    unit->m_idx = idx;
}

/////////////////////////////////////////////////////////////



void
load(InterfaceTable *inTable)
{
    ft = inTable;
#define DefineDelayUnit(name)                                           \
    (*ft->fDefineUnit)(#name, sizeof(name), (UnitCtorFunc)&name##_Ctor, \
                       (UnitDtorFunc)&DelayUnit_Dtor, 0);

    DefineDtorCantAliasUnit(AllpassLPN);
    DefineSimpleUnit(SimpleFDN2);
    DefineSimpleUnit(SimpleFDN4);
    DefineDtorUnit(Hadamard8Delay);
    DefineDtorUnit(Hadamard4Delay);
    DefineDtorUnit(Velvet1);
    DefineDtorUnit(Nested2AllpassN);
    DefineDtorUnit(Ortho3Delay);
    DefineSimpleUnit(JotFDN3);
    DefineSimpleUnit(JotFDN3LPAP);
    DefineDtorUnit(LBCF);
    DefineDtorUnit(FbTapper);
    DefineDtorUnit(NestedLBCF);
    DefineSimpleUnit(BufNestedLBCF);
    DefineSimpleUnit(VelvetTaps);
    DefineSimpleUnit(BufLBCF);
    DefineSimpleUnit(BufLAPF);
    DefineDtorUnit(Velvet);
    DefineDtorUnit(AllpassLatc);
    DefineDtorUnit(AllpassLatcLP);
    DefineSimpleUnit(BufOneTap);
    DefineSimpleUnit(BufOneTapC);
    DefineSimpleUnit(BufWrite);
    DefineSimpleUnit(Allpass2N);
    DefineSimpleUnit(NestedAllpass2N);
}



