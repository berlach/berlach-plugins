/* APFChain.cpp
   Copyright by Bjoern Erlach

   Chains of first order allpass filters.
*/

#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>

#include <math.h>

static InterfaceTable *ft;

 
struct APFChain : public Unit
{
  float mem[256];
  int nfilters;
};


struct APFChainMod : public Unit
{
  float m_mem[512];
  int nfilters;
};


extern "C"  {

void load(InterfaceTable *inTable);

int api_version(void);

void APFChain_Ctor(APFChain *unit);
void APFChain_next(APFChain *unit, int inNumSamples);
void APFChain_next_a(APFChain *unit, int inNumSamples);

void APFChainMod_Ctor(APFChainMod *unit);
void APFChainMod_next(APFChainMod *unit, int inNumSamples);

}



int api_version(void) 
{ 
    return sc_api_version; 
}

static inline float lagrange3_read (float del, float *dline, int idx, int mask)
{
    int idel = (int) del;
    float D = idel - (float) del;
    float x0 = dline[(idx-idel-1)&mask]; 
    float x1 = dline[(idx-idel) & mask];
    float x2 = dline[(idx-idel+1) & mask];
    float x3 = dline[(idx-idel+2) & mask];
    float dm1 = D - 1.0f;
    float dm2 = D - 2.0f;
    float dp1 = D + 1.0f;
    float t1 = dm1*D * 0.166666666666666666667f;
    float t2 = dm2*dp1 * 0.5f;
    return (dp1*x3 -dm2*x0) * t1 + (dm1*x1 - D*x2) * t2;
}

static inline float lagrange3 (float D, float x0, float x1, float x2, float x3)
{
    float dm1 = D - 1.0f;
    float dm2 = D - 2.0f;
    float dp1 = D + 1.0f;
    float t1 = dm1*D * 0.166666666666666666667f;
    float t2 = dm2*dp1 * 0.5f;
    return (dp1*x3 -dm2*x0) * t1 + (dm1*x1 - D*x2) * t2;
}


void
APFChain_Ctor (APFChain *unit)
{

    if (INRATE(1) == calc_FullRate) {
        SETCALC(APFChain_next_a);
    } else {
        SETCALC(APFChain_next);
    }

    memset(unit->mem, 0, sizeof(float) * 256);
    unit->nfilters = (int) ZIN0(2);  
    ZOUT0(0) = 0.f;
}


void
APFChain_next (APFChain *unit, int inNumSamples)
{
    float *out = OUT(0);
    float *in = IN(0);
    float *outp = out;
    //float *inp = in;
    float k = ZIN0(1);
    int nfilters = unit->nfilters;
    float *y1 = unit->mem;

    for (int n=0; n<nfilters; n++) {
        float y1v = y1[n];
        for (int _j=0; _j<inNumSamples; _j++) {
            float x;
            x = *in++;
            float y0 = x + k * y1v;
            x = -k * y0 + y1v;
            y1v = y0;
            *out++ = x;
        }
        y1[n] = y1v;
        in = outp;
        out = outp;
    }
  
    for (int n=0; n<nfilters; n++){
        y1[n] = zapgremlins(y1[n]);
    }
}


void
APFChain_next_a (APFChain *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *kin = ZIN(1);
    int nfilters = unit->nfilters;
    float *y1 = unit->mem;

    LOOP(inNumSamples,
         float x;
         float k;
         x = ZXP(in);
         k = ZXP(kin);

         for (int n=0; n<nfilters; n++) {
             float y0 = x - k * y1[n];
             x = k * y0 + y1[n];
             y1[n] = y0;
         }

         ZXP(out) = x;
         );

    for (int n=0; n<nfilters; n++){
        y1[n] = zapgremlins(y1[n]);
    }

}


void
APFChainMod_Ctor (APFChainMod *unit)
{
    SETCALC(APFChainMod_next);
    memset(unit->m_mem, 0, sizeof(float) * 512);
    unit->nfilters = (int) ZIN0(2);  
    ZOUT0(0)  = 0.f;
}


void
APFChainMod_next (APFChainMod *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float lam = ZIN0(1);
    float *delin = ZIN(3);
    int nfilters = unit->nfilters;
    float *m = unit->m_mem;

    for (int k=0; k<inNumSamples; ++k) {
        float x = ZXP(in);
        float x1 = lam * (m[1] - x) + m[0];
        m[0] = x;
        for (int k = 1; k<nfilters; ++k) {
            x = lam * (m[k+1] - x1) + m[k];
            m[k] = x1;
            x1 = x;
        }
        m[nfilters] = x1;
        float del = ZXP(delin);
        int di = del;
        float frac = del - (float) di;
        x = lagrange3(frac, m[di], m[di+1], m[di+2], m[di+3]);
        ZXP(out) = x;
    }
}

void
APFChainMod_next_a (APFChainMod *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *lamin = ZIN(1);
    float *delin = ZIN(3);
    int nfilters = unit->nfilters;
    float *m = unit->m_mem;

    for (int k=0; k<inNumSamples; ++k) {
        float lam = ZXP(lamin);
        float x = ZXP(in);
        float x1 = lam * (m[1] - x) + m[0];
        m[0] = x;
        for (int k = 1; k<nfilters; ++k) {
            x = lam * (m[k+1] - x1) + m[k];
            m[k] = x1;
            x1 = x;
        }
        m[nfilters] = x1;
        float del = ZXP(delin);
        int di = del;
        float frac = del - (float) di;
        x = lagrange3(frac, m[di], m[di+1], m[di+2], m[di+3]);
        ZXP(out) = x;
    }

}

////////////////////////////////////////////////////////

void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleCantAliasUnit(APFChain);
    DefineSimpleUnit(APFChainMod);
}



