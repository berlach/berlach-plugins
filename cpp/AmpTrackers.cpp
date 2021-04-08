/* AmpAllpass.cpp
   Copyright (C) 2018 by Bjoern Erlach. 

   Several amplitude/trackers to build flexible Compressors.
*/

#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>

static InterfaceTable *ft;


struct AmplitudeHold : public Unit
{
    float m_amp;
    float m_smooth;
    int m_sampshold;
};


struct AmpAllpass : public Unit
{
    float m_m11;
    float m_m12;
    float m_m21;
    float m_m22;
    float m_m31;
    float m_m32;
};

struct SmoothRMS : public Unit
{
    float m_sum;
    float m_mem[512];
    float m_m1;
    float m_att;
    float m_rel;
    int m_idx;
};


struct UpDownCompressor : public Unit
{
    float m_gain;
};


struct UpCompressor : public Unit
{
    float m_gain;
};


struct DownCompressor : public Unit
{
    float m_gain;
};



extern "C"  {

    void load(InterfaceTable *inTable);
    int api_version(void);

    void AmplitudeHold_Ctor(AmplitudeHold *unit);
    void AmplitudeHold_next(AmplitudeHold *unit, int inNumSamples);

    void AmpAllpass_Ctor(AmpAllpass *unit);
    void AmpAllpass_next(AmpAllpass *unit, int inNumSamples);

    void SmoothRMS_Ctor(SmoothRMS *unit);
    void SmoothRMS_next(SmoothRMS *unit, int inNumSamples);

    void UpDownCompressor_Ctor(UpDownCompressor *unit);
    void UpDownCompressor_next(UpDownCompressor *unit, int inNumSamples);

    void UpCompressor_Ctor(UpCompressor *unit);
    void UpCompressor_next(UpCompressor *unit, int inNumSamples);

    void DownCompressor_Ctor(DownCompressor *unit);
    void DownCompressor_next(DownCompressor *unit, int inNumSamples);
}


int api_version(void) 
{ 
    return sc_api_version; 
}


////////////////////////////////////////////////////////////////



union ficast {
	int i;
	float f;
};


// power of 2 approximation
static inline float pow2_appx(float x)
{
  union ficast fi;
  float frac;
  int ipart = ((int)(x));
  if (x<0.f) {
	  frac = 1.f+x-(float)ipart;
	  ipart--;
  } else {
	  frac = x-(float)ipart;
  }
  float mant = frac * 
	  (frac * (0.0791252185430753179895546623f * 
		   frac + 0.2249463110926884779061651898f) 
	   + 0.6959284703642362179820679557f);
  fi.i = (((127+ipart)<<23)|(int)(8388607.f*mant));
  return fi.f;
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


#define db2lin(x) pow2_appx((x)*3.321928094887363f/20.f)
#define lin2db(x) (((float)20.f/log2(10.0))*fastlog2(x))


static inline float softclip_gain(float gain, float tame, float tamer)
{
    gain *= tamer;
    gain = 1.f - (1.f/(1.f+gain*(gain*(gain+1.f)+1.f)));
    gain *= tame;
    return gain;
}


static inline float upward_gain_computer(float xdb, float thres, float rratio, float knee)
{
    float gaindif = (xdb-thres);
    float gaindif2 = 2.f*gaindif;
    float ydb;
    // gain computer (upward)
    if (gaindif2 < -knee) {
        ydb = thres + gaindif * rratio;
    } else if (fabsf(gaindif2) <= knee) {
        float tmp = (thres - xdb + knee*0.5f);
        ydb = xdb + (1.f-rratio) * tmp * tmp / (2.f*knee);
    } else {
        ydb = xdb;
    }
    return  ydb - xdb;
}

static inline float downward_gain_computer(float xdb, float thres, float rratio, float knee)
{
    float gaindif = (xdb-thres);
    float gaindif2 = 2.f*gaindif;
    float ydb;
    // gain computer (downward)
    if (gaindif2 < -knee) {
        ydb = xdb;
    } else if (fabsf(gaindif2) <= knee) {
        float tmp = (xdb - thres + knee*0.5f);
        ydb = xdb + (rratio-1.f) * tmp * tmp / (2.f*knee);
    } else {
        ydb = thres + gaindif * rratio;
    }
    return  ydb - xdb;
}


static inline float updown_gain_computer(float xdb,
					 float upthres, float uprratio, float upknee,
					 float dothres, float dorratio, float doknee
					 )
{
    float ydb;
    if (xdb >= (dothres-0.5f*doknee)) {
        float gaindif = (xdb-dothres);
        float gaindif2 = 2.f*gaindif;
        if (fabsf(gaindif2) <= doknee) {
            float tmp = (xdb - dothres + doknee*0.5f);
            ydb = xdb + (dorratio-1.f) * tmp * tmp / (2.f*doknee);
        } else {
            ydb = dothres + gaindif * dorratio;
        }
        return  ydb - xdb;
    } else {
        float gaindif = (xdb-upthres);
        float gaindif2 = 2.f*gaindif;
        float ydb;
        // gain computer (upward)
        if (gaindif2 < -upknee) {
            ydb = upthres + gaindif * uprratio;
        } else if (fabsf(gaindif2) <= upknee) {
            float tmp = (upthres - xdb + upknee*0.5f);
            ydb = xdb + (1.f-uprratio) * tmp * tmp / (2.f*upknee);
        } else {
            ydb = xdb;
        }
        return  ydb - xdb;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////

void
AmplitudeHold_Ctor (AmplitudeHold *unit)
{
    unit->m_amp = 0.f;
    unit->m_amp = 0.f;
    SETCALC(AmplitudeHold_next);
    OUT0(0) = 0.f;
}

void
AmplitudeHold_next (AmplitudeHold *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float risecoef = ZIN0(1);
    float fallcoef = ZIN0(2);
    int hold = (int) ZIN0(3);
    float amp = unit->m_amp;
    float smooth = unit->m_smooth;
    float sampshold = unit->m_sampshold;


    LOOP(inNumSamples,
         float x = ZXP(in);
         x = fabsf(x);
         if (x > amp) {
             sampshold = hold;
             amp = x - risecoef * (x - amp);
         } else {
             if (sampshold>0) {
                 sampshold--;
             } else {
                 amp = x - fallcoef * (x - amp);
             }
         }
         ZXP(out) = amp;
         );

    unit->m_smooth = smooth;
    unit->m_amp = amp;
    unit->m_sampshold = sampshold;
}


//////////////////////////////////////////////////////////////////////////////////////////

void
UpDownCompressor_Ctor (UpDownCompressor *unit)
{
    SETCALC(UpDownCompressor_next);
    unit->m_gain = 0.0; 
    OUT0(0) = 0.f;
}

void
UpDownCompressor_next (UpDownCompressor *unit, int inNumSamples)
{
    float *out1 = ZOUT(0);
    float *in1 = ZIN(0);
    float upthres = ZIN0(1);
    float upratio = ZIN0(2);
    float upknee = ZIN0(3);
    float dothres = ZIN0(4);
    float doratio = ZIN0(5);
    float doknee = ZIN0(6);
    float attack = ZIN0(7);
    float release = ZIN0(8);
    float makeup = ZIN0(9);
    float tame = ZIN0(10);
    float gain = unit->m_gain;
    float tamer = 1.f/tame;
    float attcoef = 1.f-expf(-1.f / (attack * SAMPLERATE));
    float relcoef = 1.f-expf(-1.f / (release * SAMPLERATE));
    float uprratio = 1.f/upratio;
    float dorratio = 1.f/doratio;

    LOOP(inNumSamples,
         float xl = ZXP(in1);
         float xmax = fabsf(xl);
         float xdb;
         float yl;
         float rawgain;
         float gainmul;
           
         xdb = lin2db(xmax);
           
         yl = updown_gain_computer(xdb, upthres, uprratio, upknee,
                                   dothres, dorratio, doknee);
           
         gain += (yl<gain?attcoef:relcoef) * (yl - gain);    

         if (gain > 0.f) {
             gain = softclip_gain(gain, tame, tamer);
         }

         rawgain = db2lin(gain);
         gainmul = rawgain * makeup;
           
         ZXP(out1) = gainmul;
         );
  

    unit->m_gain = gain;
}

/////////////////////////////////////////////////////////////////////

void
UpCompressor_Ctor (UpCompressor *unit)
{
    SETCALC(UpCompressor_next);
    unit->m_gain = 0.0; 
    OUT0(0) = 0.f;
}


void
UpCompressor_next (UpCompressor *unit, int inNumSamples)
{
    float *out1 = ZOUT(0);
    float *in1 = ZIN(0);
    float upthres = ZIN0(1);
    float upratio = ZIN0(2);
    float upknee = ZIN0(3);
    float attack = ZIN0(4);
    float release = ZIN0(5);
    float makeup = ZIN0(6);
    float tame = ZIN0(7);
    float gain = unit->m_gain;
    float tamer = 1.f/tame;
    float attcoef = 1.f-expf(-1.f / (attack * SAMPLERATE));
    float relcoef = 1.f-expf(-1.f / (release * SAMPLERATE));
    float uprratio = 1.f/upratio;

    LOOP(inNumSamples,
         float xl = ZXP(in1);
         float xmax = fabsf(xl);
         float xdb;
         float yl;
         float rawgain;
         float gainmul;
       
         xdb = lin2db(xmax);
         yl = upward_gain_computer(xdb, upthres, uprratio, upknee);
         gain += (yl<gain?attcoef:relcoef) * (yl - gain);    
         gain = softclip_gain(gain, tame, tamer);
         rawgain = db2lin(gain);
         gainmul = rawgain * makeup;
         ZXP(out1) = gainmul;           
         );  

    unit->m_gain = gain;
}

/////////////////////////////////////////////////////////////////////

void
DownCompressor_Ctor (DownCompressor *unit)
{
    SETCALC(DownCompressor_next);
    unit->m_gain = 0.0; 
    OUT0(0) = 0.f;
}



void
DownCompressor_next (DownCompressor *unit, int inNumSamples)
{
    float *out1 = ZOUT(0);
    float *in1 = ZIN(0);
    float upthres = ZIN0(1);
    float upratio = ZIN0(2);
    float upknee = ZIN0(3);
    float attack = ZIN0(4);
    float release = ZIN0(5);
    float makeup = ZIN0(6);
    float gain = unit->m_gain;
    float attcoef = 1.f-expf(-1.f / (attack * SAMPLERATE));
    float relcoef = 1.f-expf(-1.f / (release * SAMPLERATE));
    float uprratio = 1.f/upratio;

    LOOP(inNumSamples,
         float xl = ZXP(in1);
         float xmax = fabsf(xl);
         float xdb;
         float yl;
         float rawgain;
         float gainmul;
           
         xdb = lin2db(xmax);
         yl = downward_gain_computer(xdb, upthres, uprratio, upknee);
         gain += (yl<gain?attcoef:relcoef) * (yl - gain);    
         rawgain = db2lin(gain);
         gainmul = rawgain * makeup;
         ZXP(out1) = gainmul;           
         );


    unit->m_gain = gain;
}


/////////////////////////////////////////////////////////////////////

void
AmpAllpass_Ctor (AmpAllpass *unit)
{
    SETCALC(AmpAllpass_next);
    unit->m_m11 = 0.f;
    unit->m_m12 = 0.f;
    unit->m_m21 = 0.f;
    unit->m_m22 = 0.f;
    unit->m_m31 = 0.f;
    unit->m_m32 = 0.f;

    OUT0(0) = 0.f;
}

void
AmpAllpass_next (AmpAllpass *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float m11 = unit->m_m11;
    float m12 = unit->m_m12;
    float m21 = unit->m_m21;
    float m22 = unit->m_m22;
    float m31 = unit->m_m31;
    float m32 = unit->m_m32;

    LOOP(inNumSamples,
         float x;
         float ymax;
         float y;
         float tmp;
         x = ZXP(in);
         ymax = fabsf(x);
       
         tmp = 0.98f * m12 + x;
         y = 0.85f * fabsf(0.98 * tmp - m12);
         m12 = m11;
         m11 = tmp;
         if (y > ymax) {
             ymax = y;
         }

         tmp = 0.92f * m22 + x;
         y = 0.85f * fabsf(0.92 * tmp - m22);
         m22 = m21;
         m21 = tmp;
         if (y > ymax) {
             ymax = y;
         }

         tmp = 0.956f * m32 + x;
         y = 0.85f * fabsf(0.956f * tmp - m32);
         m32 = m31;
         m31 = tmp;
         if (y > ymax) {
             ymax = y;
         }

         ZXP(out) = ymax
         );

    unit->m_m11 = zapgremlins(m11);
    unit->m_m12 = zapgremlins(m12);
    unit->m_m21 = zapgremlins(m21);
    unit->m_m22 = zapgremlins(m22);
    unit->m_m31 = zapgremlins(m31);
    unit->m_m32 = zapgremlins(m32);
}

/////////////////////////////////////////////////////////////////////////////



void
SmoothRMS_Ctor (SmoothRMS *unit)
{
    SETCALC(AmpAllpass_next);
    unit->m_sum = 0.f;
    memset(unit->m_mem, 0, sizeof(float)*512);
    unit->m_m1 = 0.f;
    unit->m_idx = 0;

    OUT0(0) = 0.f;
}

void
SmoothRMS_next (SmoothRMS *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    int winsize = (int)ZIN0(1);
    float att = ZIN0(2);
    float rel = ZIN0(3);
    float m1 = unit->m_m1;
    float *mem = unit->m_mem;
    float sum = unit->m_sum;
    int idx = unit->m_idx;

    if ((att!=unit->m_att) || (rel!=unit->m_rel)) {
        unit->m_att = att;
        unit->m_rel = rel;
        float sr = SAMPLERATE;
        att = 1-exp(-1.0/(sr*att));
        rel = 1-exp(-1.0/(sr*rel));
    } else {
        att = unit->m_att;
        rel = unit->m_rel;
    }

    LOOP(inNumSamples,
         float x;
         float xx;
         x = ZXP(in);
         xx = x*x;
         mem[idx&511] = xx;
         sum += xx - mem[(idx-winsize)&511];
         m1 += (sum>m1?att:rel) * (sum-m1);
         ZXP(out) = sqrtf(m1);
         );

    unit->m_m1 = zapgremlins(m1);
    unit->m_sum = sum;
    unit->m_idx = idx;
}
    


void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(AmpAllpass);
    DefineSimpleUnit(SmoothRMS);
    DefineSimpleUnit(AmplitudeHold);
    DefineSimpleUnit(UpDownCompressor);
    DefineSimpleUnit(DownCompressor);
    DefineSimpleUnit(UpCompressor);
}
