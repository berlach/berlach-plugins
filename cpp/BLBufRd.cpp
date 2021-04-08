/* BLBufRd.cpp
   Copyright (C) 2005 by Bjoern Erlach. 

   Bandlimited buffer reader.  Can playback audio at much higher rate with
   strongly reduced aliasing.
*/

#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>

#include <math.h>
#include <stdio.h>
#include "sinctab.h"

static InterfaceTable *ft;

struct BLBufRd : public Unit
{
	SndBuf *m_buf;
	float m_fbufnum;
};


#define SIMPLE_GET_BUF                                  \
    float fbufnum  = ZIN0(0);                           \
    if (fbufnum != unit->m_fbufnum) {                   \
        uint32 bufnum = (int)fbufnum;                   \
        World *world = unit->mWorld;                    \
        if (bufnum >= world->mNumSndBufs) bufnum = 0;   \
        unit->m_fbufnum = fbufnum;                      \
        unit->m_buf = world->mSndBufs + bufnum;         \
    }                                                   \
    SndBuf *buf = unit->m_buf;                          \

#define CHECK_BUF                               \
    if (!bufData) {                             \
        unit->mDone = true;                     \
        ClearUnitOutputs(unit, inNumSamples);   \
        return;                                 \
    }


extern "C"  {

    void load(InterfaceTable *inTable);

    void BLBufRd_Ctor(BLBufRd *unit);
    void BLBufRd_next(BLBufRd *unit, int inNumSamples);

}

void
BLBufRd_Ctor (BLBufRd *unit)
{
    unit->m_fbufnum = -1.f;
    SETCALC(BLBufRd_next);
    ZOUT0(0) = 0.f;
}


#define SINC(iphs, fphs) ((1.f-fphs) * sinctab[iphs] + fphs * sinctab[iphs+1])
#define SINCN(iphs) (sinctab[iphs])

void
BLBufRd_next (BLBufRd *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *phsin = ZIN(1);
    float ratio = ZIN0(2);

    SIMPLE_GET_BUF;


    float *bufData = buf->data;
    int bufFrames = buf->frames;

    CHECK_BUF;
  
    LOOP(inNumSamples, 
         float sum = 0.f;
         float ratio256 = 256.f*(1.f/ratio);
         float phs = ZXP(phsin);

         phs = (phs >= (float) bufFrames ? (float) (bufFrames-2) : phs);
         fprintf(stderr, "%f\n", phs);
         int pos = (int) phs;
         int nsamps = (int) (ratio * 2.f);
         float frac = phs - (float) pos;
         int n;

         n = sc_min(nsamps, pos);
         for (int k=0; k<n; k++) {
             float distance = frac+k;
             float fphs = ratio256*distance;
             int iphs = lrintf(fphs);
             sum += bufData[pos - k] * SINCN(iphs);
         }
       
         frac = 1.f - frac;
         pos ++;

         n = sc_min(nsamps, bufFrames-pos);       
         for (int k=0; k<n; k++) {
             float distance = frac+k;
             float fphs = ratio256*distance;
             int iphs = lrintf(fphs);
             sum += bufData[pos + k] * SINCN(iphs);
         }

         ZXP(out) = sum/ratio;
         );
}
////////////////////////////////////////////////////////////////////////////


void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(BLBufRd);
}



