/* BerlachFilters.cpp
   Copyright (C) by Bjoern Erlach

   MoogImproved - Based on "An improved analgo model of the Moog ladder filter"
   ICASSP-88. by Stefano D'Angelo and￼Vesa Välimäki
￼
   BMoogFF - Moog filter emulation. Variation of MoogFF with noise insertion and
             inputs for tiny modulations of the coefficients to create subtle extra
             dynamic
             
   Korg35 - Korg filter based on Will Pirkle's design with extra non-linearity added
            and internal noise source

   Korg35OS - Oversampled version of Korg35

   WDFPeak - Peaking equalizer implemented with wave digital structure

   WDFPeakClip - WDFPeak with internal clipping

   ShelfWDF - cascade of wave digital structured shelfing filters.

   LPWDF - wdf low pass filter

   HPWDF1 - wdf highpass filter

   WDFBP - wdf bandpass filter

   BRLPF - Variation of RLPF resonant lowpass filter

   Several first order lowpass implementations with different
   placing of the normalization giving different time varying behaviours:

   LPOnePole
   LP1
   LP1NL
   LP1GP
   LP1GP2
*/

#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>
#include "tanhtab.h"
#include "Oversamp.h"

static InterfaceTable *ft;

#include "OSNonlinearClippers.h"

static const float uninitializedControl = std::numeric_limits<float>::quiet_NaN();


struct MoogImproved : public Oversamp3
{
    double m_tV0;
    double m_tV1;
    double m_tV2;
    double m_tV3;
    double m_V0;
    double m_V1;
    double m_V2;
    double m_V3;
    float m_freq;
};


struct BMoogFF : public Unit
{
    float m_freq, m_k;
    double m_b0, m_a1; 
    // Resonant freq and corresponding vals; stored because we need to compare against prev vals
    double m_wcD;
    float m_mod1mem, m_mod2mem;
    double m_s1, m_s2, m_s3, m_s4; // 1st order filter states
};


struct Korg35OS : public Oversamp4
{
    float m_alpha;
    float m_alpha0;
    float m_beta2;
    float m_beta3;
    float m_s1, m_s2, m_s3;
    float m_sat, m_Fc, m_K;
};


struct Korg35 : public Unit
{
    float m_alpha;
    float m_alpha0;
    float m_beta2;
    float m_beta3;
    float m_s1, m_s2, m_s3;
    float m_sat, m_Fc, m_K;
};


struct WDFPeak : public Unit
{
    float m_mem1;
    float m_mem2;
    float m_m1;
    float m_m2;
    float m_lastfreq;
};


struct WDFBP : public Unit
{
    float m_mem1;
    float m_mem2;
    float m_m1;
    float m_m2;
    float m_lastfreq;
};


struct WDFPeakClip : public Unit
{
    float m_mem1;
    float m_mem2;
    float m_m1;
    float m_m2;
    float m_lastfreq;
};


struct LPOnePole : public Unit
{
    float m_y1;
    float m_p;
    float m_lastfreq;
};

struct LP1 : public Unit
{
    float m_y1;
    float m_p;
    float m_lastfreq;
};

struct LP1NL : public Unit
{
    float m_y1;
    float m_p;
    float m_lastfreq;
};

struct LP1GP : public Unit
{
    float m_y1;
    float m_p;
    float m_lastfreq;
};


struct LP1GP2 : public Unit
{
    float m_y1;
    float m_p;
    float m_lastfreq;
};

struct LPWDF1 : public Unit
{
    float m_mem;
    float m_m;
    float m_lastfreq;
};


struct HPWDF1 : public Unit
{
    float m_mem;
    float m_m;
    float m_lastfreq;
};

struct ShelfWDF : public Unit
{
    float m_mem[6];
    float m_m1[6];
    float m_m2[6];
};


struct BRLPF : public Unit
{
    float m_freq, m_reson;
    double m_y1, m_y2, m_a0, m_b1, m_b2;
    float m_mem[4096];
    int m_idx;
};



extern "C"  {
    void load(InterfaceTable *inTable);
    int api_version(void);

    void BRLPF_next(BRLPF *unit, int inNumSamples);
    void BRLPF_Ctor(BRLPF *unit);

    void BMoogFF_next(BMoogFF *unit, int inNumSamples);
    void BMoogFF_next_a(BMoogFF *unit, int inNumSamples);
    void BMoogFF_Ctor(BMoogFF* unit);

    void MoogImproved_Ctor (MoogImproved *unit);
    void MoogImproved_next (MoogImproved *unit, int inNumSamples);

    void LP1GP_Ctor(LP1GP *unit);
    void LP1GP_next(LP1GP *unit, int inNumSamples);
    void LP1GP_next_a(LP1GP *unit, int inNumSamples);

    void LP1NL_Ctor(LP1NL *unit);
    void LP1NL_next(LP1NL *unit, int inNumSamples);
    void LP1NL_next_a(LP1NL *unit, int inNumSamples);

    void LP1GP2_Ctor(LP1GP2 *unit);
    void LP1GP2_next(LP1GP2 *unit, int inNumSamples);
    void LP1GP2_next_a(LP1GP2 *unit, int inNumSamples);

    void LPOnePole_Ctor(LPOnePole *unit);
    void LPOnePole_next(LPOnePole *unit, int inNumSamples);
    void LPOnePole_next_a(LPOnePole *unit, int inNumSamples);

    void LP1_Ctor(LP1 *unit);
    void LP1_next(LP1 *unit, int inNumSamples);
    void LP1_next_a(LP1 *unit, int inNumSamples);

    void LPWDF1_Ctor(LPWDF1 *unit);
    void LPWDF1_next(LPWDF1 *unit, int inNumSamples);
    void LPWDF1_next_a(LPWDF1 *unit, int inNumSamples);

    void HPWDF1_Ctor(HPWDF1 *unit);
    void HPWDF1_next(HPWDF1 *unit, int inNumSamples);
    void HPWDF1_next_a(HPWDF1 *unit, int inNumSamples);

    void ShelfWDF_Ctor(ShelfWDF *unit);
    void ShelfWDF_next(ShelfWDF *unit, int inNumSamples);
    void ShelfWDF_next_a(ShelfWDF *unit, int inNumSamples);

    void Korg35_Ctor(Korg35 *unit);
    void Korg35_next(Korg35 *unit, int inNumSamples);
    // with external noise source
    void Korg35_next_a(Korg35 *unit, int inNumSamples);

    void Korg35OS_Ctor(Korg35OS *unit);
    void Korg35OS_next(Korg35OS *unit, int inNumSamples);

    void WDFBP_Ctor(WDFBP *unit);
    void WDFBP_next(WDFBP *unit, int inNumSamples);

    void WDFPeak_Ctor(WDFPeak *unit);
    void WDFPeak_next(WDFPeak *unit, int inNumSamples);
    void WDFPeak_next_a(WDFPeak *unit, int inNumSamples);

    void WDFPeakClip_Ctor(WDFPeakClip *unit);
    void WDFPeakClip_next(WDFPeakClip *unit, int inNumSamples);
    void WDFPeakClip_next_a(WDFPeakClip *unit, int inNumSamples);
}


int api_version(void) 
{ 
    return sc_api_version; 
}


////////////////////////////////////////////////////////////////////////


// Themal voltage (26 milliwats at room temperature)
#define VT 0.312


void MoogImproved_Ctor (MoogImproved *unit)
{
    OVERSAMPLE3_INIT;
    unit->m_tV0 = 0.0;
    unit->m_tV1 = 0.0;
    unit->m_tV2 = 0.0;
    unit->m_tV3 = 0.0;
    unit->m_V0 = 0.0;
    unit->m_V1 = 0.0;
    unit->m_V2 = 0.0;
    unit->m_V3 = 0.0;
    unit->m_freq = ZIN0(1);

    SETCALC(MoogImproved_next);
    ZOUT0(0) = 0;
}


void MoogImproved_next (MoogImproved *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float drive = ZIN0(3);
    float noise = ZIN0(4);
    float newfreq = ZIN0(1);
    float res = ZIN0(2);
    float c1 = ZIN0(5);
    float c2 = ZIN0(6);
    float c3 = ZIN0(7);
    float c4 = ZIN0(8);
    float c5 = ZIN0(9);
    float freq = unit->m_freq;
    int nsamps = inNumSamples * 3;
    float freqdif = (newfreq-freq)/nsamps;
    float srinv = SAMPLEDUR / 3;
    float sr = SAMPLERATE * 3;

    double tV0 = unit->m_tV0;
    double tV1 = unit->m_tV1;
    double tV2 = unit->m_tV2;
    double tV3 = unit->m_tV3;
    double V0 = unit->m_V0;
    double V1 = unit->m_V1;
    double V2 = unit->m_V2;
    double V3 = unit->m_V3;

    RGET

    UPSAMPLE3;

    for (int i=0; i<nsamps; ++i) {
        double dV0;
        double dV1;
        double dV2;
        double dV3;

        //double fc =  freq * srinv;
        //double fc2 = fc * fc;
        //double fc3 = fc * fc * fc;
        //double fcr = 1.8730 * fc3 + 0.4955 * fc2 - 0.6490 * fc + 0.9988;
        //acr = -3.9364 * fc2 + 1.8409 * fc + 0.9968;
        //double g = (1.0 - exp(-((2 * MOOG_PI) * f * fcr))) / ; 
        //double g2 = -(-0.69346 * fc3 - 0.59515 * fc2 + 3.2937 * fc - 1.0072);

        double x = (0.5 * M_PI * freq) * srinv;
        double g = 4.0 * M_PI * VT * freq * (1.0 - x) / (1.0 + x);
        //double g2 = (x * (1.0 - x) / (1.0 + x)) / (2.0*VT);

        dV0 = -g * (tanh((drive * domemoff[i] + res * V3) / (2.0 * VT)) + tV0);
        V0 += (dV0 + dV0) / (2.0 * sr);
        tV0 = tanh_a(V0 / (2.0 * VT));         

        dV1 = g * (tV0 - tV1);
        V1 += (dV1 + dV1 + noise * frand2(s1,s2,s3)) / (2.0 * sr);
        tV1 = tanh_a(V1 / (2.0 * VT));
         
        dV2 = g * (tV1 - tV2);
        V2 += (dV2 + dV2+ noise * frand2(s1,s2,s3)) / (2.0 * sr);
        tV2 = tanh_a(V2 / (2.0 * VT));
         
        dV3 = g * (tV2 - tV3);
        V3 += (dV3 + dV3) / (2.0 * sr);
        tV3 = tanh_a(V3 / (2.0 * VT));
         
        domemoff[i] = c5 * V3 + c4 * tV2 + c3 * tV1 + c2 * tV0 + c1 * domemoff[i];

        freq += freqdif;
    }

    DOWNSAMPLE3;

    RPUT

    unit->m_tV0 = zapgremlins(tV0);
    unit->m_tV1 = zapgremlins(tV1);
    unit->m_tV2 = zapgremlins(tV2);
    unit->m_tV3 = zapgremlins(tV3);
    unit->m_V0 = zapgremlins(V0);
    unit->m_V1 = zapgremlins(V1);
    unit->m_V2 = zapgremlins(V2);
    unit->m_V3 = zapgremlins(V3);
    unit->m_freq = freq;
}



/////////////////////////////////////////////////////////////////////


void BMoogFF_Ctor(BMoogFF* unit)
{
    if (INRATE(1) == calc_FullRate) {
	SETCALC(BMoogFF_next_a);
    } else {
	SETCALC(BMoogFF_next);
    }
    // initialize the unit generator state variables.
    //unit->m_freq = uninitializedControl;
    unit->m_freq = 200;
    unit->m_k    = IN0(2);
    unit->m_mod1mem = 0.f;
    unit->m_mod2mem = 0.f;
    unit->m_s1 = 0.f;
    unit->m_s2 = 0.f;
    unit->m_s3 = 0.f;
    unit->m_s4 = 0.f;

    // calculate one sample of output.
    BMoogFF_next(unit, 1);
}

void BMoogFF_next(BMoogFF *unit, int inNumSamples)
{
    float *out = ZOUT(0);

    float *in = ZIN(0);
    float k = IN0(2);
    float noise = ZIN0(4);

    k = k > 4.f? 4.f : (k<0.f ? 0.f : k);

    // Load state from the struct
    double ss1 = unit->m_s1;
    double ss2 = unit->m_s2;
    double ss3 = unit->m_s3;
    double ss4 = unit->m_s4;

    // Reset filter state if requested
    if(IN0(3)>0)
        ss1 = ss2 = ss3 = ss4 = 0.f;

    double a1 = unit->m_a1, b0 = unit->m_b0; // Filter coefficient parameters
    double o, u; // System's null response, loop input


    RGET
	// Update filter coefficients, but only if freq changes since it involves some expensive operations

	float freqIn = IN0(1);
    if(unit->m_freq != freqIn) {
        //Print("Updated freq to %g\n", freq);
        double wcD=unit->m_wcD;
        double T = SAMPLEDUR;
        wcD = 2.0 * tan ( T * M_PI * freqIn ) * SAMPLERATE;
        if(wcD<0)
            wcD = 0; // Protect against negative cutoff freq
        double TwcD = T*wcD;
        b0 = (float)(TwcD/(TwcD + 2.));
        a1 = (float)((TwcD - 2.)/(TwcD + 2.));
        unit->m_freq = freqIn;
        unit->m_b0 = b0;
        unit->m_a1 = a1;
        unit->m_wcD = wcD;
    }

    if (unit->m_k == k) {
        LOOP1(inNumSamples,
              // compute loop values
              o = ss4 + b0*(ss3 + b0*(ss2 + b0*ss1));
              double ins = ZXP(in);
              double outs = (b0*b0*b0*b0*ins + o) * sc_reciprocal(1.0 + b0*b0*b0*b0*k);
              ZXP(out) = outs;
              u = ins - k*outs;

              // update 1st order filter states
              double past = u;
              double future = b0*past + ss1;
              ss1 = b0*past - a1*future;

              past = future + noise*frand2(s1,s2,s3);
              past = 0.5 * (sc_softclip(past) + past);
              future = b0*past + ss2;
              ss2 = b0*past - a1*future;

              past = future + noise*frand2(s1,s2,s3);
              past = 0.5 * (sc_softclip(past) + past);
              //past = sc_softclip(past);
              future = b0*past + ss3;
              ss3 = b0*past - a1*future;

              ss4 = b0*future - a1*outs;
              )
            } else {
        float new_k = k;
        float old_k = unit->m_k;
        float slope_k = CALCSLOPE(new_k, old_k);
        k = old_k;

        LOOP1(inNumSamples,
              // compute loop values
              o = ss4 + b0*(ss3 + b0*(ss2 + b0*ss1));
              double ins = ZXP(in);
              double outs = (b0*b0*b0*b0*ins + o) * sc_reciprocal(1.0 + b0*b0*b0*b0*k);
              ZXP(out) = outs;
              u = ins - k*outs;

              // update 1st order filter states
              double past = u;
              double future = b0*past + ss1;
              ss1 = b0*past - a1*future;

              past = future;
              future = b0*past + ss2;
              ss2 = b0*past - a1*future;

              past = future;
              future = b0*past + ss3;
              ss3 = b0*past - a1*future;

              ss4 = b0*future - a1*outs;
              k += slope_k;
              );
        unit->m_k = new_k;
    }

    RPUT

    // Store state
    unit->m_s1 = ss1;
    unit->m_s2 = ss2;
    unit->m_s3 = ss3;
    unit->m_s4 = ss4;
}


void BMoogFF_next_a(BMoogFF *unit, int inNumSamples)
{
    float *out = ZOUT(0);

    float *in = ZIN(0);
    float k = IN0(2);
    float noise = ZIN0(4);
    float sat = ZIN0(5);
    float fuck1 = ZIN0(6);
    float fuck2 = ZIN0(7);

    float mod1mem = unit->m_mod1mem;
    float mod2mem = unit->m_mod2mem;

    k = k > 4.f? 4.f : (k<0.f ? 0.f : k);

    // Load state from the struct
    double ss1 = unit->m_s1;
    double ss2 = unit->m_s2;
    double ss3 = unit->m_s3;
    double ss4 = unit->m_s4;

    // Reset filter state if requested
    if(IN0(3)>0)
        ss1 = ss2 = ss3 = ss4 = 0.f;

    double a1 = unit->m_a1, b0 = unit->m_b0; // Filter coefficient parameters
    double o, u; // System's null response, loop input


    RGET
	// Update filter coefficients, but only if freq changes since it involves some expensive operations

	float *freqIn = ZIN(1);
	

    if (unit->m_k == k) {
        LOOP1(inNumSamples,
              //Print("Updated freq to %g\n", freq);
              float freq = ZXP(freqIn);
              double wcD=unit->m_wcD;
              double T = SAMPLEDUR;
              wcD = 2.0 * tan ( T * M_PI * freq ) * SAMPLERATE;
              if(wcD<0)
                  wcD = 0; // Protect against negative cutoff freq
              double TwcD = T*wcD;
              b0 = (float)(TwcD/(TwcD + 2.));
              a1 = (float)((TwcD - 2.)/(TwcD + 2.));

              // compute loop values
              o = ss4 + b0*(ss3 + b0*(ss2 + b0*ss1));
              double ins = ZXP(in);
              double outs = (b0*b0*b0*b0*ins + o) * sc_reciprocal(1.0 + b0*b0*b0*b0*k);
              ZXP(out) = outs;
              u = ins - k*outs;
		      
              // update 1st order filter states
              double past = u;
              // fuckery
              mod1mem = (0.1 * fuck1) + 0.9 * mod1mem;
              double b0e = b0 + mod1mem;
              double future = b0e*past + ss1;
              ss1 = b0e*past - a1*future;
		      
              past = future + noise*frand2(s1,s2,s3);
              past = sat * 0.5 * (sc_softclip(past) + past);
              future = b0*past + ss2;
              ss2 = b0*past - a1*future;
		      
              past = future + noise*frand2(s1,s2,s3);
              past = 0.5 * (sc_softclip(past) + past);
              //past = sc_softclip(past);
              mod2mem = (0.1 * fuck2) + 0.9 * mod2mem;
              b0e = b0 + mod2mem;
              future = b0e*past + ss3;
              ss3 = b0e*past - a1*future;
              ss4 = b0*future - a1*outs;
              )
            } else {
        float new_k = k;
        float old_k = unit->m_k;
        float slope_k = CALCSLOPE(new_k, old_k);
        k = old_k;

        LOOP1(inNumSamples,
              float freq = ZXP(freqIn);
              double wcD=unit->m_wcD;
              double T = SAMPLEDUR;

              wcD = 2.0 * tan ( T * M_PI * freq ) * SAMPLERATE;
              if(wcD<0)
                  wcD = 0; // Protect against negative cutoff freq
              double TwcD = T*wcD;
              b0 = (float)(TwcD/(TwcD + 2.));
              a1 = (float)((TwcD - 2.)/(TwcD + 2.));

              // compute loop values
              o = ss4 + b0*(ss3 + b0*(ss2 + b0*ss1));
              double ins = ZXP(in);
              double outs = (b0*b0*b0*b0*ins + o) * sc_reciprocal(1.0 + b0*b0*b0*b0*k);
              ZXP(out) = outs;
              u = ins - k*outs;
		      
              // update 1st order filter states
              double past = u;
              // fuckery
              mod1mem = (0.1 * fuck1) + 0.9 * mod1mem;
              double b0e = b0 + mod1mem;
              double future = b0e*past + ss1;
              ss1 = b0e*past - a1*future;
		      
              past = future + noise*frand2(s1,s2,s3);
              past = sat * 0.5 * (sc_softclip(past) + past);
              future = b0*past + ss2;
              ss2 = b0*past - a1*future;
		      
              past = future + noise*frand2(s1,s2,s3);
              past = 0.5 * (sc_softclip(past) + past);
              //past = sc_softclip(past);
              mod2mem = (0.1 * fuck2) + 0.9 * mod2mem;
              b0e = b0 + mod2mem;
              future = b0e*past + ss3;
              ss3 = b0e*past - a1*future;
              ss4 = b0*future - a1*outs;
              k += slope_k;
              );
        unit->m_k = new_k;
    }

    RPUT

    unit->m_mod1mem = mod1mem;
    unit->m_mod2mem = mod2mem;

    // Store state
    unit->m_s1 = ss1;
    unit->m_s2 = ss2;
    unit->m_s3 = ss3;
    unit->m_s4 = ss4;
}


//////////////////////////////////////////////////////////////////////////


static inline void Korg35_update (Korg35 *unit, float Fc, float sr, float K)
{
    // prewarp for BZT
    double wd = 2.0*M_PI*Fc;
    double T = 1.0/sr;
    double wa = (2/T)*tan(wd*T/2);
    double g = wa*T/2;
    // G - the feedforward coeff in the VA One Pole
    double G = g/(1.0 + g);
    unit->m_alpha = G;
    // set betas all are in the form of <something>/((1 + g)
    unit->m_beta2 = (K - K*G)/(1.0 + g);
    unit->m_beta3 = -1.0/(1.0 + g);
    // set m_dAlpha0 variable
    unit->m_alpha0 = 1.0/(1.0 - K*G + K*G*G); 
}


/*
  Filter based on Will Pirkle's Korg filter model.
  An additional non-linearity is inserted between the stages.
*/

void
Korg35_Ctor (Korg35 *unit)
{
    if (INRATE(4) == calc_FullRate) {
	SETCALC(Korg35_next_a);
    } else {
	SETCALC(Korg35_next);
    }

    float K = ZIN0(2);
    float Fc = ZIN0(1);    

    unit->m_s1 = 0.f;
    unit->m_s2 = 0.f;
    unit->m_s3 = 0.f;
    unit->m_sat = ZIN0(3);
    unit->m_Fc = Fc;
    unit->m_K = K;
    ZOUT0(0) = 0.f;
}


void Korg35_next (Korg35 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float s1 = unit->m_s1;
    float s2 = unit->m_s2;
    float s3 = unit->m_s3;
    float sr = SAMPLERATE;
    float saturation = unit->m_sat;
    float K = unit->m_K;
    float Fc = unit->m_Fc;

    float Fcinc = (ZIN0(1) - unit->m_Fc) / 16;
    float Kinc = (ZIN0(2) - unit->m_K) / 16;
    float satinc = (ZIN0(3) - unit->m_sat) / 16;
    float beta3, beta2, alpha0, alpha;
    
    for (int i=0; i<inNumSamples; ++i) {
	if ((i & 3) == 0) {	    
	    Korg35_update(unit, Fc, sr, K);
	    beta3 = unit->m_beta3;
	    beta2 = unit->m_beta2;
	    alpha0 = unit->m_alpha0;
	    alpha = unit->m_alpha;
	    Fc += Fcinc;
	    saturation += satinc;
	    K += Kinc;
	}

        float x = ZXP(in);
        float tmp = alpha * (x - s1);
        float lpf1 = s1 + tmp;
        s1 = tmp + lpf1;
        float u = alpha0 * (lpf1 + beta3*s3 + beta2*s2);
        //u = tanhf(saturation*u);
        u = lin_tanhtab(saturation*u);
        // feed u through lpf2
        tmp = alpha * (u - s2);
        float lpf2 = s2 + tmp;
        s2 = tmp + lpf2;
        float y = K * lpf2;
	 
        tmp = alpha * (y - s3);
        float fil3 = s3 + tmp;
        s3 = tmp + fil3;
        y *= 1.f/K;
	 
        ZXP(out) = y;
    }

    unit->m_s1 = s1;
    unit->m_s2 = s2;
    unit->m_s3 = s3;
    unit->m_K = K;
    unit->m_Fc = Fc;
    unit->m_sat = saturation;
}

void Korg35_next_a (Korg35 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *noisein = ZIN(4);
    float *out = ZOUT(0);
    float s1 = unit->m_s1;
    float s2 = unit->m_s2;
    float s3 = unit->m_s3;
    float sr = SAMPLERATE;
    float saturation = unit->m_sat;
    float K = unit->m_K;
    float Fc = unit->m_Fc;

    float Fcinc = (ZIN0(1) - unit->m_Fc) / 16;
    float Kinc = (ZIN0(2) - unit->m_K) / 16;
    float satinc = (ZIN0(3) - unit->m_sat) / 16;
    float beta3, beta2, alpha0, alpha;
    
    for (int i=0; i<inNumSamples; ++i) {
	if ((i & 3) == 0) {	    
	    Korg35_update(unit, Fc, sr, K);
	    beta3 = unit->m_beta3;
	    beta2 = unit->m_beta2;
	    alpha0 = unit->m_alpha0;
	    alpha = unit->m_alpha;
	    Fc += Fcinc;
	    saturation += satinc;
	    K += Kinc;
	}

        float x = ZXP(in);
        float tmp = alpha * (x - s1);
        float lpf1 = s1 + tmp;
        s1 = tmp + lpf1;
        float u = alpha0 * (lpf1 + beta3*s3 + beta2*s2);
        //u = tanhf(saturation*u);
        u += ZXP(noisein);
        u = lin_tanhtab(saturation*u);
        // feed u through lpf2
        tmp = alpha * (u - s2);
        float lpf2 = s2 + tmp;
        s2 = tmp + lpf2;
        float y = K * lpf2;
	 
        tmp = alpha * (y - s3);
        float fil3 = s3 + tmp;
        s3 = tmp + fil3;
        y *= 1.f/K;
	 
        ZXP(out) = y;
    }

    unit->m_s1 = s1;
    unit->m_s2 = s2;
    unit->m_s3 = s3;
    unit->m_K = K;
    unit->m_Fc = Fc;
    unit->m_sat = saturation;
}


////////////////////////////////////////////////////////////


static inline void Korg35OS_update (Korg35OS *unit, float Fc, float sr, float K)
{
    // prewarp for BZT
    double wd = 2.0*M_PI*Fc;
    double T = 1.0/sr;
    double wa = (2/T)*tan(wd*T/2);
    double g = wa*T/2;
    // G - the feedforward coeff in the VA One Pole
    double G = g/(1.0 + g);
    unit->m_alpha = G;
    // set betas all are in the form of <something>/((1 + g)
    unit->m_beta2 = (K - K*G)/(1.0 + g);
    unit->m_beta3 = -1.0/(1.0 + g);
    // set m_dAlpha0 variable
    unit->m_alpha0 = 1.0/(1.0 - K*G + K*G*G); 
}



void
Korg35OS_Ctor (Korg35OS *unit)
{
    OVERSAMPLE4_INIT;
    SETCALC(Korg35OS_next);
    float K = ZIN0(2);
    float Fc = ZIN0(1);    
    //Korg35_update(unit, Fc, SAMPLERATE, K);
    unit->m_s1 = 0.f;
    unit->m_s2 = 0.f;
    unit->m_s3 = 0.f;
    unit->m_sat = ZIN0(3);
    unit->m_Fc = Fc;
    unit->m_K = K;
    ZOUT0(0) = 0.f;
}


void Korg35OS_next (Korg35OS *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float s1 = unit->m_s1;
    float s2 = unit->m_s2;
    float s3 = unit->m_s3;
    float sr = SAMPLERATE * 4.f;
    float saturation = unit->m_sat;
    float K = unit->m_K;
    float Fc = unit->m_Fc;
    float cpoints = (float) inNumSamples;
    float Fcinc = (ZIN0(1) - unit->m_Fc) / cpoints;
    float Kinc = (ZIN0(2) - unit->m_K) / cpoints;
    float satinc = (ZIN0(3) - unit->m_sat) / cpoints;
    float beta3, beta2, alpha0, alpha;

    UPSAMPLE4;
    
    for (int i=0; i<(inNumSamples*4); ++i) {
	if ((i & 3) == 0) {	    
	    Korg35OS_update(unit, Fc, sr, K);
	    beta3 = unit->m_beta3;
	    beta2 = unit->m_beta2;
	    alpha0 = unit->m_alpha0;
	    alpha = unit->m_alpha;
	    Fc += Fcinc;
	    saturation += satinc;
	    K += Kinc;
	}

	float x = domemoff[i];
        float tmp = alpha * (x - s1);
        float lpf1 = s1 + tmp;
        s1 = tmp + lpf1;
        float u = alpha0 * (lpf1 + beta3*s3 + beta2*s2);
        //u = tanhf(saturation*u);
        u = lin_tanhtab(saturation*u);
        // feed u through lpf2
        tmp = alpha * (u - s2);
        float lpf2 = s2 + tmp;
        s2 = tmp + lpf2;
        float y = K * lpf2;
	 
        tmp = alpha * (y - s3);
        float fil3 = s3 + tmp;
        s3 = tmp + fil3;
        y *= 1.f/K;
	 
        domemoff[i] = y;
        //ZXP(out) = y;
    }

    DOWNSAMPLE4;

    unit->m_s1 = s1;
    unit->m_s2 = s2;
    unit->m_s3 = s3;
    unit->m_K = K;
    unit->m_Fc = Fc;
    unit->m_sat = saturation;
}


///////////////////////////////////////////////////////////////////////////

static inline float
tanapprox (float x)
{
    return x*(x*(x*(0.96369f*x-0.865157f)+0.53576f)+0.93f);
}


void
WDFBP_Ctor (WDFBP *unit)
{
    float sr = SAMPLERATE;
    float bw = ZIN0(2);
    float freq = ZIN0(1);
    float wtan = tanf(bw*M_PI/sr);
    float wtan2 = tanf(freq*M_PI/sr);

    unit->m_mem1 = 0.f;
    unit->m_mem2 = 0.f;
    unit->m_lastfreq = ZIN0(1);

    wtan2 = wtan2 * wtan2;
    unit->m_m1 = (1.f-wtan2)/(1.f+wtan2);
    unit->m_m2 = wtan/(1.f+wtan);
    SETCALC(WDFBP_next);
    ZOUT0(0) = 0.f;
}


void
WDFBP_next (WDFBP *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float bw = ZIN0(2);
    float G = ZIN0(3);
    float mem1 = unit->m_mem1;
    float mem2 = unit->m_mem2;
    float sr = SAMPLERATE;
    float wtan = tanf(bw*((float)M_PI)/sr);
    float wtan2 = tanf(freq*((float)M_PI)/sr);
    wtan2 = wtan2 * wtan2;
    float m1 = (1.f-wtan2)/(1.f+wtan2);
    float m2 = wtan/(1.f+wtan);
    float m1del = (m1 - unit->m_m1) / BUFLENGTH;
    float m2del = (m2 - unit->m_m2) / BUFLENGTH;

    LOOP(inNumSamples,
         float x;
         x = ZXP(in);
         float s1 = m1 * (mem2 - mem1);
         float b12 = mem2 + s1;
         float b22 = mem1 + s1;
         mem2 = b22;
         s1 = (x + b12);
         float s2 = -m2 * s1;
         float b11 = x + s2;
         mem1 = -(s1 + b11 + s2);
         ZXP(out) = G * b11 + s2;
         m1 += m1del;
         m2 += m2del;
         );

    unit->m_m1 = m1;
    unit->m_m2 = m2;
    unit->m_mem1 = zapgremlins(mem1);    
    unit->m_mem2 = zapgremlins(mem2);    
}


void
WDFPeak_Ctor (WDFPeak *unit)
{
    float sr = SAMPLERATE;
    float bw = ZIN0(2);
    float freq = ZIN0(1);

    unit->m_mem1 = 0.f;
    unit->m_mem2 = 0.f;
    unit->m_lastfreq = ZIN0(1);

    if (INRATE(1) == calc_FullRate) {
	float wtan = tanapprox(bw*M_PI/sr);
	float wtan2 = tanapprox(freq*M_PI/sr);
	wtan2 = wtan2 * wtan2;
	unit->m_m1 = (1.f-wtan2)/(1.f+wtan2);
	unit->m_m2 = wtan/(1.f+wtan);
	SETCALC(WDFPeak_next_a);
    } else {
	float wtan = tanf(bw*M_PI/sr);
	float wtan2 = tanf(freq*M_PI/sr);
	wtan2 = wtan2 * wtan2;
	unit->m_m1 = (1.f-wtan2)/(1.f+wtan2);
	unit->m_m2 = wtan/(1.f+wtan);
	SETCALC(WDFPeak_next);
    }

    ZOUT0(0) = 0.f;
}


void
WDFPeak_next (WDFPeak *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float bw = ZIN0(2);
    float G = ZIN0(3);
    float mem1 = unit->m_mem1;
    float mem2 = unit->m_mem2;
    float sr = SAMPLERATE;
    float wtan = tanf(bw*((float)M_PI)/sr);
    float wtan2 = tanf(freq*((float)M_PI)/sr);
    wtan2 = wtan2 * wtan2;
    float m1 = (1.f-wtan2)/(1.f+wtan2);
    float m2 = wtan/(1.f+wtan);
    float m1del = (m1 - unit->m_m1) / BUFLENGTH;
    float m2del = (m2 - unit->m_m2) / BUFLENGTH;

    LOOP(inNumSamples,
         float x;
         x = ZXP(in);
         float s1 = m1 * (mem2 - mem1);
         float b12 = mem2 + s1;
         float b22 = mem1 + s1;
         mem2 = b22;
         s1 = (x + b12);
         float s2 = -m2 * s1;
         float b11 = x + s2;
         mem1 = -(s1 + b11 + s2);
         ZXP(out) = b11 + G * s2;
         m1 += m1del;
         m2 += m2del;
         );

    unit->m_m1 = m1;
    unit->m_m2 = m2;
    unit->m_mem1 = zapgremlins(mem1);    
    unit->m_mem2 = zapgremlins(mem2);    
}

void
WDFPeak_next_a (WDFPeak *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freqin = ZIN(1);
    float *bwin = ZIN(2);
    float *Gin = ZIN(3);
    float mem1 = unit->m_mem1;
    float mem2 = unit->m_mem2;
    float sr = SAMPLERATE;

    LOOP(inNumSamples,
         float x;
         float freq = ZXP(freqin);
         float bw = ZXP(bwin);
         float G = ZXP(Gin);
         float wtan = tanapprox(bw*((float)M_PI)/sr);
         float wtan2 = tanapprox(freq*((float)M_PI)/sr);
         wtan2 = wtan2 * wtan2;
         float m1 = (1.f-wtan2)/(1.f+wtan2);
         float m2 = wtan/(1.f+wtan);
         x = ZXP(in);
         float s1 = m1 * (mem2 - mem1);
         float b12 = mem2 + s1;
         float b22 = mem1 + s1;
         mem2 = b22;
         s1 = (x + b12);
         float s2 = -m2 * s1;
         float b11 = x + s2;
         mem1 = -(s1 + b11 + s2);
         ZXP(out) = b11 + G * s2;
         );

    unit->m_mem1 = zapgremlins(mem1);    
    unit->m_mem2 = zapgremlins(mem2);    
}


////////////////////////////////////////////////////////////////////////////////


void
WDFPeakClip_Ctor (WDFPeakClip *unit)
{
    unit->m_mem1 = 0.f;
    unit->m_mem2 = 0.f;
    unit->m_lastfreq = ZIN0(1);
    float sr = SAMPLERATE;
    float bw = ZIN0(2);
    float freq = ZIN0(1);
    float wtan = tanapprox(bw*M_PI/sr);
    float wtan2 = tanapprox(freq*M_PI/sr);
    wtan2 = wtan2 * wtan2;
    unit->m_m1 = (1.f-wtan2)/(1.f+wtan2);
    unit->m_m2 = wtan/(1.f+wtan);
    if (INRATE(1) == calc_FullRate) {
	SETCALC(WDFPeakClip_next_a);
    } else {
	SETCALC(WDFPeakClip_next);
    }
    ZOUT0(0) = 0.f;
}

void
WDFPeakClip_next (WDFPeakClip *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float bw = ZIN0(2);
    float G = ZIN0(3);
    float saturation = ZIN0(4);
    float scale = ZIN0(5);
    float mem1 = unit->m_mem1;
    float mem2 = unit->m_mem2;
    float sr = SAMPLERATE;
    float wtan = tanapprox(bw*((float)M_PI)/sr);
    float wtan2 = tanapprox(freq*((float)M_PI)/sr);
    wtan2 = wtan2 * wtan2;
    float m1 = (1.f-wtan2)/(1.f+wtan2);
    float m2 = wtan/(1.f+wtan);

    LOOP(inNumSamples,
         float x;
         x = ZXP(in);
         float s1 = m1 * (mem2 - mem1);
         float b12 = mem2 + s1;
         float b22 = mem1 + s1;
         mem2 = b22;
         s1 = (x + b12);
         float s2 = -m2 * s1;
         float b11 = x + s2;
         mem1 = -(s1 + b11 + s2);
         mem1 = lin_tanhtab(saturation*mem1) * scale;
         ZXP(out) = b11 + G * s2;
         );

    unit->m_mem1 = zapgremlins(mem1);    
    unit->m_mem2 = zapgremlins(mem2);    
}


void
WDFPeakClip_next_a (WDFPeakClip *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freqin = ZIN(1);
    float *bwin = ZIN(2);
    float *Gin = ZIN(3);
    float saturation = ZIN0(4);
    float scale = ZIN0(5);
    float mem1 = unit->m_mem1;
    float mem2 = unit->m_mem2;
    float sr = SAMPLERATE;

    LOOP(inNumSamples,
         float x;
         float freq = ZXP(freqin);
         float bw = ZXP(bwin);
         float G = ZXP(Gin);
         float wtan = tanapprox(bw*((float)M_PI)/sr);
         float wtan2 = tanapprox(freq*((float)M_PI)/sr);
         wtan2 = wtan2 * wtan2;
         float m1 = (1.f-wtan2)/(1.f+wtan2);
         float m2 = wtan/(1.f+wtan);
         x = ZXP(in);
         float s1 = m1 * (mem2 - mem1);
         float b12 = mem2 + s1;
         float b22 = mem1 + s1;
         mem2 = b22;
         s1 = (x + b12);
         float s2 = -m2 * s1;
         float b11 = x + s2;
         mem1 = -(s1 + b11 + s2);
         mem1 = lin_tanhtab(saturation*mem1) * scale;
         ZXP(out) = b11 + G * s2;
         );

    unit->m_mem1 = zapgremlins(mem1);    
    unit->m_mem2 = zapgremlins(mem2);    
}


/////////////////////////////////////////////////////////////////////////////////


void
LPOnePole_Ctor (LPOnePole *unit)
{

    float p;
    unit->m_lastfreq = ZIN0(1);
    if (INRATE(1) == calc_FullRate) {
	SETCALC(LPOnePole_next_a);
	p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
    } else {
	unit->m_p = p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
	SETCALC(LPOnePole_next);
    }
    unit->m_y1 = (1.f-p)*ZIN0(0);
    ZOUT0(0) = unit->m_y1;
}

void
LPOnePole_next (LPOnePole *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float p = unit->m_p;
    float y1 = unit->m_y1;

    if (unit->m_lastfreq != freq) {
        float newp = (1.f-2.f*tanf((freq*pi/SAMPLERATE)));
        float pinc = (newp-p)/inNumSamples;
        unit->m_p=newp;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             ZXP(out) = y1 = (1.f-p)*ZXP(in)+ p*y1;
             p += pinc;
             );
    } else {
        LOOP(inNumSamples,
             ZXP(out) = y1 = (1.f-p)*ZXP(in)+ p*y1;
             );
    }

    unit->m_y1 = zapgremlins(y1);    
}


void
LPOnePole_next_a (LPOnePole *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float y1 = unit->m_y1;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float p = (1.f-2.f*tanapprox((ZXP(freq)*pioversrate)));
         ZXP(out) = y1 = (1.f-p)*ZXP(in)+ p*y1;
         );

    unit->m_y1 = zapgremlins(y1);  
}

////////////////////////////////////////////////////////////////////////////


void
LP1_Ctor (LP1 *unit)
{

    float p;
    unit->m_lastfreq = ZIN0(1);
    if (INRATE(1) == calc_FullRate) {
	SETCALC(LP1_next_a);
	p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
    } else {
	unit->m_p = p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
	SETCALC(LP1_next);
    }
    unit->m_y1 = (1.f-p)*ZIN0(0);
    ZOUT0(0) = unit->m_y1;
}

void
LP1_next (LP1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float p = unit->m_p;
    float y1 = unit->m_y1;
    //fprintf(stderr, "p:%f\n", p);
  
    if (unit->m_lastfreq != freq) {
        float newp = (1.f-2.f*tanf((freq*pi/SAMPLERATE)));
        float pinc = (newp-p)/inNumSamples;
        unit->m_p=newp;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
             ZXP(out) = tmp + y1;
             y1 = tmp;
             p += pinc;
             );
    } else {
        LOOP(inNumSamples,
             float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
             ZXP(out) = tmp + y1;
             y1 = tmp;
             );
    }
  
    unit->m_y1 = zapgremlins(y1);  
}


void
LP1_next_a (LP1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float y1 = unit->m_y1;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float p = (1.f-2.f*tanapprox((ZXP(freq)*pioversrate)));
         float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
         ZXP(out) = tmp + y1;
         y1 = tmp;
         );

    unit->m_y1 = zapgremlins(y1);  
}


//////////////////////////////////////////////////////////////////////////


void
LP1GP_Ctor (LP1GP *unit)
{

    float p;
    unit->m_lastfreq = ZIN0(1);
    if (INRATE(1) == calc_FullRate) {
	SETCALC(LP1GP_next_a);
	p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
    } else {
	unit->m_p = p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
	SETCALC(LP1GP_next);
    }
    unit->m_y1 = ZIN0(0);
    ZOUT0(0) = unit->m_y1 * (1.f-p);
}

void
LP1GP_next (LP1GP *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float p = unit->m_p;
    float y1 = unit->m_y1;

    if (unit->m_lastfreq != freq) {
        float newp = (1.f-2.f*tanf((freq*pi/SAMPLERATE)));
        float pinc = (newp-p)/inNumSamples;
        unit->m_p=newp;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             y1 = ZXP(in)+ p*y1;
             ZXP(out) = (1.f-p) * y1;
             p += pinc;
             );
    } else {
        LOOP(inNumSamples,
             y1 = ZXP(in)+ p*y1;
             ZXP(out) = (1.f-p) *  y1;
             );
    }
  
    unit->m_y1 = zapgremlins(y1);  
}


void
LP1GP_next_a (LP1GP *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float y1 = unit->m_y1;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float p = (1.f-2.f*tanapprox((ZXP(freq)*pioversrate)));
         y1 = ZXP(in)+ p*y1;
         ZXP(out) = (1.f-p) * y1;
         );

    unit->m_y1 = zapgremlins(y1);  
}



//////////////////////////////////////////////////////////////////////////


void
LP1GP2_Ctor (LP1GP2 *unit)
{

    float p;
    unit->m_lastfreq = ZIN0(1);

    if (INRATE(1) == calc_FullRate) {
	SETCALC(LP1GP2_next_a);
    } else {
	SETCALC(LP1GP2_next);
    }

    unit->m_y1 = ZIN0(0);
    ZOUT0(0) = unit->m_y1 * (1.f-p);
}

void
LP1GP2_next (LP1GP2 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float a = ZIN0(1);
    float b = ZIN0(2);
    float p = unit->m_p;
    float y1 = unit->m_y1;
    float g = (1.f-a) / (1.f+b);
    float g2 = sqrtf(g);

    LOOP(inNumSamples,
         float tmp = g2 * ZXP(in)+ a*y1;
         ZXP(out) = g2 * (tmp + b*y1);
         y1 = tmp;
         //p += pinc;
         );
  
    unit->m_y1 = zapgremlins(y1);  
}


void
LP1GP2_next_a (LP1GP2 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *ain = ZIN(1);
    float *bin = ZIN(2);
    float p = unit->m_p;
    float y1 = unit->m_y1;

    LOOP(inNumSamples,
         float a= ZXP(ain);
         float b= ZXP(bin);
         float g = (1.f-a) / (1.f+b);
         float g2 = sqrtf(g);
         float tmp = g2 * ZXP(in)+ a*y1;
         ZXP(out) = g2 * (tmp + b*y1);
         y1 = tmp;
         //p += pinc;
         );
  
    unit->m_y1 = zapgremlins(y1);  
}


////////////////////////////////////////////////////////////////////////////


static inline float wdf_coeff (float freq, float PIoverSR)
{	
    float w = freq * PIoverSR;
    float tw2 = tanf(w);
    return tw2/(1.f+tw2);
}

static inline float wdf_coeff_apx (float freq, float PIoverSR)
{	
    float w = freq * PIoverSR;
    float tw2 = tanapprox(w);
    return tw2/(1.f+tw2);
}


void
LPWDF1_Ctor (LPWDF1 *unit)
{

    float in = ZIN0(0);
    unit->m_lastfreq = ZIN0(1);
    float m = wdf_coeff(unit->m_lastfreq, pi/SAMPLERATE);

    if (INRATE(1) == calc_FullRate) {
	SETCALC(LPWDF1_next_a);
    } else {
	unit->m_m = m;
	SETCALC(LPWDF1_next);
    }

    float b3 = -m * in; 
    unit->m_mem = -2.f * (in + b3);
    ZOUT0(0) = b3;
}

void
LPWDF1_next (LPWDF1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float m = unit->m_m;
    float mem = unit->m_mem;

    if (unit->m_lastfreq != freq) {
        float newm = wdf_coeff(freq, pi/SAMPLERATE);
        float minc = (newm-m)/inNumSamples;
        unit->m_m=newm;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             float a1 = ZXP(in);
             float a2 = -mem;
             float tmp1 = a1 + a2; 
             float b3 = -m * tmp1; mem = -(tmp1 + a1 + b3 + b3);
             ZXP(out) = b3;
             m += minc;
             );
    } else {
        LOOP(inNumSamples,
             float a1 = ZXP(in);
             float a2 = -mem;
             float tmp1 = a1 + a2; 
             float b3 = -m * tmp1; mem = -(tmp1 + a1 + b3 + b3);
             ZXP(out) = b3;
             );
    }

    unit->m_mem = zapgremlins(mem);  
}


void
LPWDF1_next_a (LPWDF1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float mem = unit->m_mem;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float m = wdf_coeff_apx(ZXP(freq), pioversrate);
         float a1 = ZXP(in);
         float a2 = -mem;
         float tmp1 = a1 + a2; 
         float b3 = -m * tmp1; mem = -(tmp1 + a1 + b3 + b3);
         ZXP(out) = b3;
         );

    unit->m_mem = zapgremlins(mem);
}

//////////////////////////////////////////////////////////////////////


void
LP1NL_Ctor (LP1NL *unit)
{

    float p;
    unit->m_lastfreq = ZIN0(1);
    if (INRATE(1) == calc_FullRate) {
	SETCALC(LP1NL_next_a);
	p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
    } else {
	unit->m_p = p = (1.f-2.f*tanf((unit->m_lastfreq*pi/SAMPLERATE)));
	SETCALC(LP1NL_next);
    }
    unit->m_y1 = (1.f-p)*ZIN0(0);
    ZOUT0(0) = unit->m_y1;
}

void
LP1NL_next (LP1NL *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float sattype1 = ZIN0(2);
    float blend1 = ZIN0(3);
    float p = unit->m_p;
    float y1 = unit->m_y1;
    //fprintf(stderr, "p:%f\n", p);
  
    if (unit->m_lastfreq != freq) {
        float newp = (1.f-2.f*tanf((freq*pi/SAMPLERATE)));
        float pinc = (newp-p)/inNumSamples;
        unit->m_p=newp;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
             float dx = blend1 * clipper(tmp, sattype1) + (1.f-blend1) * tmp;
             ZXP(out) = dx + y1;
             y1 = dx;
             p += pinc;
             );
    } else {
        LOOP(inNumSamples,
             float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
             ZXP(out) = tmp + y1;
             y1 = tmp;
             );
    }
  
    unit->m_y1 = zapgremlins(y1);  
}


void
LP1NL_next_a (LP1NL *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float sattype1 = ZIN0(2);
    float blend1 = ZIN0(3);
    float y1 = unit->m_y1;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float p = (1.f-2.f*tanapprox((ZXP(freq)*pioversrate)));
         float tmp = 0.5f*(1.f-p)*ZXP(in)+ p*y1;
         float dx = blend1 * clipper(tmp, sattype1) + (1.f-blend1) * tmp;

         ZXP(out) = dx + y1;
         y1 = dx;
         );

    unit->m_y1 = zapgremlins(y1);  
}


////////////////////////////////////////////////////////////////////////////

void
HPWDF1_Ctor (HPWDF1 *unit)
{

    float in = ZIN0(0);
    unit->m_lastfreq = ZIN0(1);
    float m = wdf_coeff(unit->m_lastfreq, pi/SAMPLERATE);

    if (INRATE(1) == calc_FullRate) {
	SETCALC(HPWDF1_next_a);
    } else {
	unit->m_m = m;
	SETCALC(HPWDF1_next);
    }

    float b3 = -m * in; 
    float b1 = b3 + in;
    unit->m_mem = -(in + b1 + b3);
    ZOUT0(0) = b1;
}

void
HPWDF1_next (HPWDF1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float freq = ZIN0(1);
    float m = unit->m_m;
    float mem = unit->m_mem;

    if (unit->m_lastfreq != freq) {
        float newm = wdf_coeff(freq, pi/SAMPLERATE);
        float minc = (newm-m)/inNumSamples;
        unit->m_m=newm;
        unit->m_lastfreq=freq;
        LOOP(inNumSamples,
             float a1 = ZXP(in);
             float a2 = -mem;
             float tmp1 = a1 + a2; 
             float b3 = -m * tmp1; 
             float b1 = b3 + a1;
             mem = -(tmp1 + b1 + b3);
             ZXP(out) = b1;
             m += minc;
             );
    } else {
        LOOP(inNumSamples,
             float a1 = ZXP(in);
             float a2 = -mem;
             float tmp1 = a1 + a2; 
             float b3 = -m * tmp1;
             float b1 = b3 + a1;
             mem = -(tmp1 + b1 + b3);
             ZXP(out) = b1;
             );
    }

    unit->m_mem = zapgremlins(mem);  
}


void
HPWDF1_next_a (HPWDF1 *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float *freq = ZIN(1);
    float mem = unit->m_mem;
    float pioversrate = pi/SAMPLERATE;

    LOOP(inNumSamples,
         float m = wdf_coeff_apx(ZXP(freq), pioversrate);
         float a1 = ZXP(in);
         float a2 = -mem;
         float tmp1 = a1 + a2; 
         float b3 = -m * tmp1;
         float b1 = b3 + a1;
         mem = -(tmp1 + b1 + b3);
         ZXP(out) = b1;
         );

    unit->m_mem = zapgremlins(mem);
}

//////////////////////////////////////////////////////////////////////


void
ShelfWDF_Ctor (ShelfWDF *unit)
{
    if (INRATE(1) == calc_FullRate) {
	SETCALC(ShelfWDF_next_a);
    } else {
	SETCALC(ShelfWDF_next_a);
    }
    memset(unit->m_mem, 0, sizeof(float) * 6);
    ZOUT0(0) = 0.f;
}


void
ShelfWDF_next_a (ShelfWDF *unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    int n_lowshelfs = ZIN0(1);
    int n_filters = ZIN0(2);
    float *m_mem = unit->m_mem;
    float pioversrate = pi/SAMPLERATE;
    float mem[6];
    float *freqsIn[6];
    float *gainsIn[6];

    for (int n=0; n<n_filters; ++n) {
        mem[n] = m_mem[n];
        freqsIn[n] = ZIN(n+3);
        gainsIn[n] = ZIN(n+n_filters+3);
    }

    LOOP(inNumSamples,
         float a1 = ZXP(in);

         for (int n=0; n<n_lowshelfs; ++n) {
             float m = wdf_coeff_apx(ZXP(freqsIn[n]), pioversrate);
             float g = ZXP(gainsIn[n]);
             float b3;
             float b1;
             float a2 = -mem[n];
             float tmp1 = a1 + a2;b3 = -m * tmp1;b1 = b3 + a1;mem[n] = -(tmp1 + b1 + b3);
             a1 = b1 + g*b3;
         }

         for (int n=n_lowshelfs; n<n_filters; ++n) {
             float m = wdf_coeff_apx(ZXP(freqsIn[n]), pioversrate);
             float g = ZXP(gainsIn[n]);
             float b3;
             float b1;
             float a2 = -mem[n];
             float tmp1 = a1 + a2;b3 = -m * tmp1;b1 = b3 + a1;mem[n] = -(tmp1 + b1 + b3);
             a1 = g*b1 + b3;
         }

         ZXP(out) = a1;
         );

    for (int n=0; n<n_filters; ++n) {
        unit->m_mem[n] = zapgremlins(mem[n]);
    }
}


///////////////////////////////////////////////////////////////////////////////////


void BRLPF_Ctor(BRLPF* unit)
{
    SETCALC(BRLPF_next);
        
    memset(unit->m_mem, 0, sizeof(float) * 4096);
    unit->m_idx = 0;
    unit->m_a0 = 0.f;
    unit->m_b1 = 0.f;
    unit->m_b2 = 0.f;
    unit->m_y1 = 0.f;
    unit->m_y2 = 0.f;
    unit->m_freq = uninitializedControl;
    unit->m_reson = uninitializedControl;
    ZOUT0(0) = 0.f;
}


void BRLPF_next(BRLPF* unit, int inNumSamples)
{
    float *out = ZOUT(0);
    float *in = ZIN(0);
    float reson = ZIN0(2);
    float coeff = ZIN0(1);
    float apcoeff = ZIN0(3);
    float ffcoeff = ZIN0(4);
    float norm = ZIN0(5);
    float damp = ZIN0(6);
    int idel = 24;
    int idx = unit->m_idx;
    float *mem = unit->m_mem;
    double y0;
    double y1 = unit->m_y1;
    double y2 = unit->m_y2;
    double a0 = unit->m_a0;
    double b1 = unit->m_b1;
    double b2 = unit->m_b2;

    LOOP(inNumSamples,
         float tmp;
         float x = ZXP(in);
         float apsig = mem[(idx-idel)&4095];

         tmp = apcoeff * apsig + norm * y1;
         apsig -= ffcoeff * tmp;
         mem[idx&4095] = tmp;

         y2 += damp * coeff * apsig;
         y1 += damp * coeff * (x - y2 - reson * y1);             
         ZXP(out) = y2;
         idx++;
         );

    unit->m_idx = idx;
    unit->m_y1 = zapgremlins(y1);
    unit->m_y2 = zapgremlins(y2);
}



void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(LPOnePole);
    DefineSimpleUnit(LP1);
    DefineSimpleUnit(LP1NL);
    DefineSimpleUnit(LPWDF1);
    DefineSimpleUnit(HPWDF1);
    DefineSimpleUnit(ShelfWDF);
    DefineSimpleUnit(Korg35);
    DefineSimpleUnit(Korg35OS);
    DefineSimpleUnit(WDFPeak);
    DefineSimpleUnit(WDFBP);
    DefineSimpleUnit(WDFPeakClip);
    DefineSimpleUnit(LP1GP);
    DefineSimpleUnit(LP1GP2);
    DefineSimpleUnit(BMoogFF);
    DefineSimpleUnit(BRLPF);
    DefineSimpleUnit(MoogImproved);
}

