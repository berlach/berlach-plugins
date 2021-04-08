/* WTO.cpp
   Copyright (C) by Bjoern Erlach. 

   Wavetable oscillators with anti-aliased hard-sync.  
   The wavetables can be generated by a python script (wto.py).
*/


#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>
#include <math.h>

static InterfaceTable *ft;

struct WTTab { int size; int offset; };

struct WTO : public Unit
{
    double m_phase;
    float m_basesize;
    int m_numtabs;
    float *m_bufdata;
    WTTab tabs[32];
    float m_trigfade;
    float m_trigval;
};

struct WTOM : public Unit
{
    double m_phase;
    float m_basesize;
    int m_numtabs;
    float *m_bufdata[8];
    WTTab m_tabs[8][32];
    int m_numBuffers;
    float m_pos;
    float m_trigfade;
    float m_trigval;
};


extern "C"  {

    void load(InterfaceTable *inTable);
    int api_version(void);

    void WTO_Ctor(WTO *unit);
    void WTO_next_aa (WTO *unit, int inNumSamples);
    void WTO_next_aaa (WTO *unit, int inNumSamples);
    void WTO_next_ka(WTO *unit, int inNumSamples);
    void WTO_next_kk(WTO *unit, int inNumSamples);

    void WTOM_Ctor(WTOM *unit);
    void WTOM_next_aa (WTOM *unit, int inNumSamples);
    void WTOM_next_aaa (WTOM *unit, int inNumSamples);
}


int api_version(void) 
{ 
    return sc_api_version; 
}



#define WTO_GET_BUF_UNLOCKED						\
    SndBuf *bufs;							\
    if (bufnum < 0)							\
	bufnum = 0;							\
    									\
    if (bufnum+1 >= world->mNumSndBufs) {				\
	int localBufNum = bufnum - world->mNumSndBufs;			\
	Graph *parent = unit->mParent;					\
	if(localBufNum <= parent->localBufNum) {			\
	    bufs = parent->mLocalSndBufs + localBufNum;			\
	} else {							\
	    bufnum = 0;							\
	    bufs = world->mSndBufs + bufnum;				\
	}								\
    } else {								\
	if (bufnum >= world->mNumSndBufs)				\
	    bufnum = 0;							\
	bufs = world->mSndBufs + sc_max(0, bufnum);			\
}


static inline float crossfade (float x1, float x2, float fade)
{
  return x1+(x2-x1)*fade;
}

// Erik de Castro Lopo suggests that c casts are slow and
// should be replaced with other methods.. I measured on my machine and
// c-casts performed best...

#define FLOAT2INT(x) ((int)(x))
#define INT2FLOAT(x) ((float)(x))

#include "wto_ctab.h"

static inline float 
calc_tabnum (float freq, float basesize, float sr)
{
    float fidx = basesize * (freq/sr);
    int idx = (int)fidx;
    float frac = fidx-(float)idx;
    return (1.0f-frac) * ctab[idx] + frac * ctab[idx+1];
}


static inline float 
fastlog2 (float x)
{
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;
    
    return y - 124.22551499f
             - 1.498030302f * mx.f 
             - 1.72587999f / (0.3520887068f + mx.f);
}


static inline float
fastlog (float x)
{
    return 0.69314718f * fastlog2 (x);
}



static inline float bspline_fade (float fade, float x1, float x2)
{
    // range is 1 .. 2
    float mul, D;
    if (fade < 1.f) {
	D = fade + 1.f;
	mul = (D*(D*((D-4.f)*D+6.f)-4.f)+1.f)/24.f;
    } else if (fade < 2.f) {
	D = fade;
	mul = (D*(D*((16.f-3.f*D)*D-24.f)+16.f)-4.f)/24.f;
    } else if (fade < 3.f) {
	D = fade-1.f;
	mul = (D*(D*(D*(3.f*D-20)+42.f)-20.f)+7.f)/24.f;
    } else {
	D = fade-2.f;
	mul = (D*(D*((8.f-D)*D-24.f)+32.f)+8.f)/24.f;
    }
    return (1.f-mul) * x1 + mul * x2;
}


void
WTOM_Ctor (WTOM *unit)
{
  int numBuffers = unit->mNumInputs - 5;

  if (INRATE(4) == calc_FullRate) {
      SETCALC(WTOM_next_aaa);
  } else {
      SETCALC(WTOM_next_aa);
  }

  unit->m_phase = 0.f;
  unit->m_pos = ZIN0(3);

  for (int i=0; i<numBuffers; ++i) {
      int bufnum = (int) ZIN0(5+i);
      SndBuf *buf = unit->mWorld->mSndBufs + bufnum;
      float *bufdata = buf->data;
      if (!bufdata) {
          SETCALC(ClearUnitOutputs);
      } else {
          // check if this is a wavetable buffer 
          unit->m_basesize = bufdata[1];
          int numtabs = (int) bufdata[2];
          int acc = 0;
          unit->m_numtabs = numtabs;
          for (int k=0; k<numtabs; k++) {
              unit->m_tabs[i][k].size = (int) bufdata[k+3];
              unit->m_tabs[i][k].offset = 3 + numtabs + acc;
              acc += (int) unit->m_tabs[i][k].size;
          }
          unit->m_bufdata[i] = bufdata;
      }
  }

  unit->m_trigfade = 5.f;
  unit->m_numBuffers = numBuffers;
  ZOUT0(0) = 0.f;
}


void
WTOM_next_aa (WTOM *unit, int inNumSamples)
{
  float *out = ZOUT(0);
  float freq;
  int numBuffers = unit->m_numBuffers;
  float *freqin = ZIN(0);
  float *phasein = ZIN(1);
  float **bufdata = unit->m_bufdata;
  float basesize = unit->m_basesize;
  float sr =  SAMPLERATE;
  float sd =  SAMPLEDUR;
  WTTab *tabs[8];
  float shift = ZIN0(2);
  float newpos = ZIN0(3);
  float pos = unit->m_pos;
  float posinc = (newpos-pos)/inNumSamples;
  float sync = ZIN0(4);
  double rphs = unit->m_phase;
  int numtabs = unit->m_numtabs;
  float tabnum;

  tabs[0] = unit->m_tabs[0];
  tabs[1] = unit->m_tabs[1];
  tabs[2] = unit->m_tabs[2];
  tabs[3] = unit->m_tabs[3];
  tabs[4] = unit->m_tabs[4];
  tabs[5] = unit->m_tabs[5];
  tabs[6] = unit->m_tabs[6];
  tabs[7] = unit->m_tabs[7];


  LOOP(inNumSamples,
       float x1; float x2; float x3; float x4;
       float o1; float o2;
       int idx; float fracphs;
       double phs; double pinc;
       int itabnum; 
       float fade;
       int tabsize1;
       int tabsize2;
       float wto1;
       freq = ZXP(freqin);
       pinc = (freq*sd);
       rphs += pinc;
       phs = rphs + ZXP(phasein);

       int w1 = (int) pos;
       float frac = pos - (float) w1;
       int w2 = (w1 + 1) % numBuffers;
       w1 %= numBuffers;
       
       tabnum = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[w1][itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
           wto1 = o1;
       } else {
	   itabnum = ((int) tabnum);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[w1][itabnum].size;
	   tabsize2 = tabs[w1][itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata[w1] + tabs[w1][itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum - INT2FLOAT(itabnum);
           wto1 = crossfade(o1, o2, fade);
       }

       tabnum = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[w2][itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
           ZXP(out) = crossfade(wto1, o1, frac);
       } else {
	   itabnum = ((int) tabnum);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[w2][itabnum].size;
	   tabsize2 = tabs[w2][itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata[w2] + tabs[w2][itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum - INT2FLOAT(itabnum);
           ZXP(out) = crossfade(wto1, crossfade(o1, o2, fade), frac);
       }

       pos += posinc;
       )
      unit->m_phase = fmod(rphs, 1.0);
  unit->m_pos = pos;
}



void
WTOM_next_aaa (WTOM *unit, int inNumSamples)
{
  float *out = ZOUT(0);
  float freq;
  int numBuffers = unit->m_numBuffers;
  float *freqin = ZIN(0);
  float *phasein = ZIN(1);
  float **bufdata = unit->m_bufdata;
  float basesize = unit->m_basesize;
  float sr =  SAMPLERATE;
  float sd =  SAMPLEDUR;
  WTTab *tabs[8];
  float shift = ZIN0(2);
  float newpos = ZIN0(3);
  float pos = unit->m_pos;
  float posinc = (newpos-pos)/inNumSamples;
  float *sync = ZIN(4);
  double rphs = unit->m_phase;
  int numtabs = unit->m_numtabs;
  float trigfade = unit->m_trigfade;
  float trigval = unit->m_trigval;


  tabs[0] = unit->m_tabs[0];
  tabs[1] = unit->m_tabs[1];
  tabs[2] = unit->m_tabs[2];
  tabs[3] = unit->m_tabs[3];
  tabs[4] = unit->m_tabs[4];
  tabs[5] = unit->m_tabs[5];
  tabs[6] = unit->m_tabs[6];
  tabs[7] = unit->m_tabs[7];


  LOOP(inNumSamples,
       float x1; float x2; float x3; float x4;
       float o1; float o2;
       float tabnum;
       int idx; float fracphs;
       double phs; double pinc;
       int itabnum; 
       float fade;
       float tabnum2;
       int tabsize1;
       int tabsize2;
       float wto1;
       float trig;
       float y;
       float yt;
       double startphs = ZXP(phasein);
       trig = ZXP(sync);
       freq = ZXP(freqin);
       pinc = (freq*sd);
       phs = rphs + startphs;

       int w1 = (int) pos;
       float frac = pos - (float) w1;
       int w2 = (w1 + 1) % numBuffers;
       w1 %= numBuffers;
       
       tabnum = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[w1][itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
           wto1 = o1;
       } else {
	   itabnum = ((int) tabnum);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[w1][itabnum].size;
	   tabsize2 = tabs[w1][itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata[w1] + tabs[w1][itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum - INT2FLOAT(itabnum);
           wto1 = crossfade(o1, o2, fade);
       }

       tabnum2 = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum2 < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[w2][itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
           y = crossfade(wto1, o1, frac);
       } else {
	   itabnum = ((int) tabnum2);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[w2][itabnum].size;
	   tabsize2 = tabs[w2][itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata[w2] + tabs[w2][itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum2 - INT2FLOAT(itabnum);
           y = crossfade(wto1, crossfade(o1, o2, fade), frac);
        }

	if (trig >= 0.1f) {	    
            float frac = (trig-0.1f);
	    trigfade = (1.f-frac);
	    trigval = trigfade * pinc + startphs - pinc;
	} 

       if (trigfade < 4.f) {
	    trigval += pinc;
            if (tabnum < 0.0) { 
                itabnum = 0;
                tabsize1 = tabs[w1][itabnum].size;
                // translate phase to fractional sample index
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                wto1 = o1;
            } else {
                itabnum = ((int) tabnum);
                if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
                tabsize1 = tabs[w1][itabnum].size;
                tabsize2 = tabs[w1][itabnum+1].size;
                
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata[w1] + tabs[w1][itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                idx = FLOAT2INT(trigval * tabsize2);
                fracphs = (trigval*tabsize2) - (float)(idx);
                idx += tabsize2;
                bufoffset = (bufdata[w1] + tabs[w1][itabnum+1].offset);
                x1 = bufoffset[(idx) % tabsize2];
                x2 = bufoffset[(idx+1) % tabsize2];
                x3 = bufoffset[(idx+2) % tabsize2];
                x4 = bufoffset[(idx+3) % tabsize2];
                o2 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                fade = tabnum - INT2FLOAT(itabnum);
                wto1 = crossfade(o1, o2, fade);
            }
            
            if (tabnum2 < 0.0) { 
                itabnum = 0;
                tabsize1 = tabs[w2][itabnum].size;
                // translate phase to fractional sample index
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                yt = crossfade(wto1, o1, frac);
            } else {
                itabnum = ((int) tabnum2);
                if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
                tabsize1 = tabs[w2][itabnum].size;
                tabsize2 = tabs[w2][itabnum+1].size;
                
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata[w2] + tabs[w2][itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                idx = FLOAT2INT(trigval * tabsize2);
                fracphs = (trigval*tabsize2) - (float)(idx);
                idx += tabsize2;
                bufoffset = (bufdata[w2] + tabs[w2][itabnum+1].offset);
                x1 = bufoffset[(idx) % tabsize2];
                x2 = bufoffset[(idx+1) % tabsize2];
                x3 = bufoffset[(idx+2) % tabsize2];
                x4 = bufoffset[(idx+3) % tabsize2];
                o2 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                fade = tabnum2 - INT2FLOAT(itabnum);
                yt = crossfade(wto1, crossfade(o1, o2, fade), frac);
            }

	    y = bspline_fade(trigfade, y, yt);
	    trigfade += 1.f;
            if (trigfade >= 4.f) {
                rphs = trigval;
            }
	}

       rphs += pinc;
       ZXP(out) = y;
       pos += posinc;
       )

  unit->m_phase = fmod(rphs, 1.0);
  unit->m_pos = pos;
  unit->m_trigfade = trigfade;
  unit->m_trigval = trigval;
}


///////////////////////////////////////////////////////////////////////////////////

void
WTO_Ctor (WTO *unit)
{
  if (INRATE(4) == calc_FullRate) {
      SETCALC(WTO_next_aaa);
  } else {
      SETCALC(WTO_next_aa);
  }
  unit->m_phase = 0.f;
  int bufnum = (int) ZIN0(0);
  SndBuf *buf = unit->mWorld->mSndBufs + bufnum;
  float *bufdata = buf->data;
  if (!bufdata) {
      SETCALC(ClearUnitOutputs);
  } else {
      // check if this is a wavetable buffer 
      unit->m_basesize = bufdata[1];
      int numtabs = (int) bufdata[2];
      int acc = 0;
      unit->m_numtabs = numtabs;
      for (int k=0; k<numtabs; k++) {
	  unit->tabs[k].size = (int) bufdata[k+3];
	  unit->tabs[k].offset = 3 + numtabs + acc;
	  acc += (int) unit->tabs[k].size;
      }
  }

  unit->m_bufdata = bufdata;
  unit->m_trigfade = 5.f;
  ZOUT0(0) = 0.f;
}


void
WTO_next_aa (WTO *unit, int inNumSamples)
{
  float *out = ZOUT(0);
  float freq;
  float *freqin = ZIN(1);
  float *phasein = ZIN(2);
  float *bufdata = unit->m_bufdata;
  float basesize = unit->m_basesize;
  float sr =  SAMPLERATE;
  float sd =  SAMPLEDUR;
  WTTab *tabs = unit->tabs;
  float shift = ZIN0(3);
  double rphs = unit->m_phase;
  int numtabs = unit->m_numtabs;
  float tabnum;

  LOOP(inNumSamples,
       float x1; float x2; float x3; float x4;
       float o1; float o2;
       int idx; float fracphs;
       double phs; double pinc;
       int itabnum; 
       float fade;
       int tabsize1;
       int tabsize2;
       freq = ZXP(freqin);
       pinc = (freq*sd);
       rphs += pinc;
       phs = rphs + ZXP(phasein);
       tabnum = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata + tabs[itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
	   ZXP(out) = o1;
       } else {
	   itabnum = ((int) tabnum);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[itabnum].size;
	   tabsize2 = tabs[itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata + tabs[itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata + tabs[itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum - INT2FLOAT(itabnum);
	   ZXP(out) = crossfade(o1, o2, fade);
       }
       )

      unit->m_phase = fmod(rphs, 1.0);
}



void
WTO_next_aaa (WTO *unit, int inNumSamples)
{
  float *out = ZOUT(0);
  float freq;
  //int bufnum = (int) ZIN0(0);
  float *freqin = ZIN(1);
  float *phasein = ZIN(2);
  float *bufdata = unit->m_bufdata;
  float basesize = unit->m_basesize;
  float sr =  SAMPLERATE;
  float sd =  SAMPLEDUR;
  WTTab *tabs = unit->tabs;
  float *trigin = ZIN(4);
  float shift = ZIN0(3);
  double rphs = unit->m_phase;
  int numtabs = unit->m_numtabs;
  float tabnum;
  float trigfade = unit->m_trigfade;
  float trigval = unit->m_trigval;

  LOOP(inNumSamples,
       float x1; float x2; float x3; float x4;
       float o1; float o2;
       int idx; float fracphs;
       double phs; double pinc;
       int itabnum; 
       float fade;
       float y;
       float trig;
       int tabsize1;
       int tabsize2;
       float startphs = ZXP(phasein);
       freq = ZXP(freqin);
       trig = ZXP(trigin);
       pinc = (freq*sd);
       phs = rphs + startphs;
       tabnum = calc_tabnum(freq, basesize, sr) - shift;

       if (tabnum < 0.0) { 
	   itabnum = 0;
	   tabsize1 = tabs[itabnum].size;
	   // translate phase to fractional sample index
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata + tabs[itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);
	   y = o1;
       } else {
	   itabnum = ((int) tabnum);
	   if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
	   tabsize1 = tabs[itabnum].size;
	   tabsize2 = tabs[itabnum+1].size;
       
	   idx = FLOAT2INT(phs * tabsize1);
	   fracphs = (phs*tabsize1) - (float)(idx);
	   idx += tabsize1;
	   float *bufoffset = (bufdata + tabs[itabnum].offset);
	   x1 = bufoffset[(idx) % tabsize1];
	   x2 = bufoffset[(idx+1) % tabsize1];
	   x3 = bufoffset[(idx+2) % tabsize1];
	   x4 = bufoffset[(idx+3) % tabsize1];
	   o1 = cubicinterp(fracphs, x1, x2, x3, x4);

	   idx = FLOAT2INT(phs * tabsize2);
	   fracphs = (phs*tabsize2) - (float)(idx);
	   idx += tabsize2;
	   bufoffset = (bufdata + tabs[itabnum+1].offset);
	   x1 = bufoffset[(idx) % tabsize2];
	   x2 = bufoffset[(idx+1) % tabsize2];
	   x3 = bufoffset[(idx+2) % tabsize2];
	   x4 = bufoffset[(idx+3) % tabsize2];
	   o2 = cubicinterp(fracphs, x1, x2, x3, x4);

	   fade = tabnum - INT2FLOAT(itabnum);
	   y = crossfade(o1, o2, fade);
       }

	if (trig >= 0.1f) {	    
            float frac = (trig-0.1f);
	    trigfade = (1.f-frac);
	    trigval = trigfade * pinc + startphs - pinc;
	} 

       if (trigfade < 4.f) {
           float yt;
	    trigval += pinc;
            if (tabnum < 0.0) { 
                itabnum = 0;
                tabsize1 = tabs[itabnum].size;
                // translate phase to fractional sample index
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata + tabs[itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                yt = o1;
            } else {
                itabnum = ((int) tabnum);
                if (itabnum >= (numtabs-2)) itabnum = numtabs-2;
                tabsize1 = tabs[itabnum].size;
                tabsize2 = tabs[itabnum+1].size;
                
                idx = FLOAT2INT(trigval * tabsize1);
                fracphs = (trigval*tabsize1) - (float)(idx);
                idx += tabsize1;
                float *bufoffset = (bufdata + tabs[itabnum].offset);
                x1 = bufoffset[(idx) % tabsize1];
                x2 = bufoffset[(idx+1) % tabsize1];
                x3 = bufoffset[(idx+2) % tabsize1];
                x4 = bufoffset[(idx+3) % tabsize1];
                o1 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                idx = FLOAT2INT(trigval * tabsize2);
                fracphs = (trigval*tabsize2) - (float)(idx);
                idx += tabsize2;
                bufoffset = (bufdata + tabs[itabnum+1].offset);
                x1 = bufoffset[(idx) % tabsize2];
                x2 = bufoffset[(idx+1) % tabsize2];
                x3 = bufoffset[(idx+2) % tabsize2];
                x4 = bufoffset[(idx+3) % tabsize2];
                o2 = cubicinterp(fracphs, x1, x2, x3, x4);
                
                fade = tabnum - INT2FLOAT(itabnum);
                yt = crossfade(o1, o2, fade);
            }


	    y = bspline_fade(trigfade, y, yt);
	    trigfade += 1.f;
            if (trigfade >= 4.f) {
                rphs = trigval;
            }
	}

       ZXP(out) = y;
       rphs += pinc;       
       )

  unit->m_phase = fmod(rphs, 1.0);
  unit->m_trigfade = trigfade;
  unit->m_trigval = trigval;
}


void
load(InterfaceTable *inTable)
{
  ft = inTable;
  DefineSimpleUnit(WTO);
  DefineSimpleUnit(WTOM);
}