/* WarpIIRLattice.cpp
   Copyright (C) by Bjoern Erlach

   Warped IIR lattice filters of arbitrary order (up to order 64)
*/
   
#include <SC_PlugIn.h>
#include <stdio.h>

static InterfaceTable *ft;

struct WarpIIRLattice : public Unit
{
    float m_k[64];
    float m_taps[64];
    float mLastVal;
};


extern "C"
{
    void load(InterfaceTable *inTable);
    
    int api_version(void);

    void WarpIIRLattice_next(WarpIIRLattice *unit, int inNumSamples);
    void WarpIIRLattice_next_a(WarpIIRLattice *unit, int inNumSamples);
    void WarpIIRLattice_Ctor(WarpIIRLattice* unit);
}


int api_version(void) 
{ 
    return sc_api_version; 
}


///////////////////////////////////////////////////////////////////////////////


void WarpIIRLattice_Ctor(WarpIIRLattice* unit)
{
    int numInputs = (unit->mNumInputs-2);
    unit->mLastVal = 0.0;
    memset(unit->m_k, 0, sizeof(float) * 64);
    memset(unit->m_taps, 0, sizeof(float) * 64);
    int allfullrate = 1;
    int ord = (numInputs-2);
    for (int j=0; j<ord; j++) {
        if (INRATE(j+2) != calc_FullRate) {
            allfullrate = 0;
            break;
        }
    }
    if (allfullrate) {
        SETCALC(WarpIIRLattice_next_a);
        fprintf(stderr, "fullrate\n");
    } else {
        SETCALC(WarpIIRLattice_next);
    }
    ZOUT0(0) = 0.f;
}


void WarpIIRLattice_next(WarpIIRLattice* unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float lambda = ZIN0(1);
    int ord = (unit->mNumInputs-2);
    float *k = unit->m_k;
    float *r = unit->m_taps;
    float vd[64];

    for (int n=0; n<ord; n++) {
        vd[n] = (ZIN0(n+2)-k[n])/inNumSamples;
    }

    LOOP(inNumSamples,
         float y;
         float input;
         float O0;
         float f;
         float b;
         float v;
         float tmp;
         float X;
         float s;

         input = ZXP(in);

         // compute the system output with 0.0 input 
         v = 1.f-lambda * lambda;
         f = 0.f;
         O0 = 0.f;
         b = v*r[0];
         O0 = k[0] * b;
         tmp = -k[0] * b;
         for (int n=1; n<ord; n++) {
             b = v*r[n] - lambda * (b - k[n-1] * f);
             O0 = O0 + k[n] * b;
             f = tmp;
             tmp = f - k[n] * b;
         }
       
         // gain measure
         b = 1.0;    f = 1.0;    X = 0.0;    s = -lambda;
         for (int n=0; n<ord; n++) {
             X = X + k[n] * s;
             b = s - k[n] * f;
             f = f - k[n] * s;
             s = -lambda * b;
         }
       
         // calculate the output
         y = (input + O0)/(1.f-X);
       
         b = y;    f = y;
         for (int n=0; n<ord; n++) {
             tmp = b + r[n] * lambda;
             b = -lambda * tmp + r[n];
             r[n] = zapgremlins(tmp);
             tmp = f - k[n] * b;
             b = b - k[n] * f;
             f = tmp;
         }
       
         ZXP(out) = y;

         for(int n=0; n<ord; n++) {
             k[n] += vd[n];
         }
         )
}



void WarpIIRLattice_next_a(WarpIIRLattice* unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float lambda = ZIN0(1);
    int ord = (unit->mNumInputs-2);
    float *k = unit->m_k;
    float *r = unit->m_taps;
    float *coeffinputs[64];

    for (int n=0; n<ord; n++) {
        coeffinputs[n] = ZIN(n+2);
    }

    LOOP(inNumSamples,
         float y;
         float input;
         float O0;
         float f;
         float b;
         float v;
         float tmp;
         float X;
         float s;

         for (int n=0; n<ord; n++) {
             k[n] = ZXP(coeffinputs[n]);
         }

         input = ZXP(in);

         // compute the system output with 0.0 input 
         v = 1.f-lambda * lambda;
         f = 0.f;
         O0 = 0.f;
         b = v*r[0];
         O0 = k[0] * b;
         tmp = -k[0] * b;
         for (int n=1; n<ord; n++) {
             b = v*r[n] - lambda * (b - k[n-1] * f);
             O0 = O0 + k[n] * b;
             f = tmp;
             tmp = f - k[n] * b;
         }
       
         // gain measure
         b = 1.0;    f = 1.0;    X = 0.0;    s = -lambda;
         for (int n=0; n<ord; n++) {
             X = X + k[n] * s;
             b = s - k[n] * f;
             f = f - k[n] * s;
             s = -lambda * b;
         }
       
         // calculate the output
         y = (input + O0)/(1.f-X);
       
         b = y;    f = y;
         for (int n=0; n<ord; n++) {
             tmp = b + r[n] * lambda;
             b = -lambda * tmp + r[n];
             r[n] = zapgremlins(tmp);
             tmp = f - k[n] * b;
             b = b - k[n] * f;
             f = tmp;
         }
       
         ZXP(out) = y;

         )
}

//////////////////////////////////////////////////////////////////////////////

void load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(WarpIIRLattice);
}
