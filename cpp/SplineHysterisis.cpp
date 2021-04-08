/* SplineHysterisis.cpp
   Copyright (C) 2019 by Bjoern Erlach.

   Processes audio input with a shaper that exhibits hysterisis.
   It works by having two different spline based shapers and a
   signal dependent blend between the shapes.
*/
   
#include <SC_PlugIn.h>
#include <SC_InterfaceTable.h>
#include <stdio.h>
#include "Oversamp.h"

static InterfaceTable *ft;

struct SplineHysterisis : public Unit
{
    float m_x1;
    float m_blend;
    float m_upper;
    float m_lower;
};


struct SplineHysterisis4 : public Oversamp4
{
    float m_x1;
    float m_blend;
    float m_upper;
    float m_lower;
};


struct SplineHysterisis4L : public Unit
{
    float m_x1;
    float m_blend;
    float m_upper;
    float m_lower;
};


extern "C"  {
    void load(InterfaceTable *inTable);
    int api_version(void);
    void SplineHysterisis_Ctor(SplineHysterisis *unit);
    void SplineHysterisis_next(SplineHysterisis *unit, int inNumSamples);

    void SplineHysterisis4_Ctor(SplineHysterisis4 *unit);
    void SplineHysterisis4_next(SplineHysterisis4 *unit, int inNumSamples);

    void SplineHysterisis4L_Ctor(SplineHysterisis4L *unit);
    void SplineHysterisis4L_next(SplineHysterisis4L *unit, int inNumSamples);
}



int api_version(void) 
{ 
    return sc_api_version; 
}


static inline void spline_coeffs(float *norma, 
                                 float *offa, 
                                 float *c0, 
                                 float *c1, 
                                 float *c2, 
                                 float *c3, 
                                 float pa0, 
                                 float ma0, 
                                 float pa1, 
                                 float ma1)
{
    float norm = (pa1-pa0);
    float ma0adj = ma0 * norm;
    float ma1adj = ma1 * norm;
    *offa = pa0 - (ma0 * pa0);
    *c3 = (ma0adj + ma1adj + 2.f*pa0 - 2.f*pa1);
    *c2 = (-2.f*ma0adj - ma1adj - 3.f*pa0 + 3.f*pa1);
    *c1 = ma0adj;
    *c0 = pa0;
    *norma = 1.f/norm;
}


void
SplineHysterisis_Ctor(SplineHysterisis *unit) 
{
    unit->m_x1 = 0.0;
    unit->m_blend = 0.5;
    unit->m_upper = 0.0;
    unit->m_lower = 0.0;
    SETCALC(SplineHysterisis_next);
    ZOUT0(0) = 0.f;
}


void
SplineHysterisis_next(SplineHysterisis *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float cross_speed = ZIN0(1); // 0.78
    float cross_curve = ZIN0(2); // 0.29
    float shift = ZIN0(3); // 0.29
    float pa0 = ZIN0(4); 
    float ma0 = ZIN0(5);
    float pa1 = ZIN0(6);
    float ma1 = ZIN0(7);
    float pb0 = ZIN0(8);
    float mb0 = ZIN0(9);
    float pb1 = ZIN0(10);
    float mb1 = ZIN0(11);

    float ca0, ca1, ca2, ca3;
    float cb0, cb1, cb2, cb3;
    float norma, offa;
    float normb, offb;
    float x1 = unit->m_x1;
    float blend = unit->m_blend;
    float lower = unit->m_lower;
    float upper = unit->m_upper;
    float lim;

    spline_coeffs(&norma, &offa, &ca0, &ca1, &ca2, &ca3, pa0, ma0, pa1, ma1);
    spline_coeffs(&normb, &offb, &cb0, &cb1, &cb2, &cb3, pb0, mb0, pb1, mb1);

    lim = -(offa+offb)*0.5f + pa1;

    for (int i=0; i<inNumSamples; ++i) {
        float xd;
        float y;
        float x = ZXP(in);
        // xd is x-x1 (input - last input)

        xd = x - x1;
        x1 = x;

        upper = 0.9993 * upper;
        lower = 0.9993 * lower;

        if (xd > 0.f) {
            float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
            float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
            float dist = uplimit - lowlimit;
            float change = dist * cross_speed*xd/(cross_curve+blend*blend);
            blend = sc_min(change + blend, uplimit);
        } else {
            float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
            float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
            float dist = uplimit - lowlimit;
            float change = dist * cross_speed*xd/(cross_curve+(1.f-blend)*(1.f-blend));
            blend = sc_max(change + blend, lowlimit);
        }

        if (x > 0) {
            float res2, res;
            if (x > pa1) {
                y = lim + (x-pa1)*ma1;
            } else {
                if (x > pb0) {
                    float t = (x-pb0)*normb;
                    res2 = (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                } else {
                    res2 = (x*mb0 + offb);
                }
                if (x > pa0) {
                    float t = (x-pa0)*norma;
                    res = (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                } else {
                    res = (x*ma0 + offa);
                }
                y = (blend * res + (1.f-blend) * res2) - (offb + offa)*0.5f;
            }
        } else {
            float res2, res;
            if (x < -pa1) {
                y = -(lim + (-x-pa1)*ma1);
            } else {
                x = -x;
                if (x > pa0) {
                    float t = (x-pa0)*norma;
                    res2 = - (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                } else {
                    res2 = - (x*ma0 + offa);
                }
                if (x > pb0) {
                    float t = (x-pb0)*normb;
                    res = - (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                } else {
                    res = - (x*mb0 + offb);
                }
                y = (blend * res + (1.f-blend) * res2) + (offb + offa)*0.5f;
            }
        }

        if (y > upper && (y < 0.8)) {
            upper = y;
        } else if (y < lower && (y > -0.8)) {
            lower = y;
        }

        ZXP(out) = y;
    }
    
    unit->m_x1 = x1; 
    unit->m_lower = lower; 
    unit->m_upper = upper; 
    unit->m_blend = blend;
}

///////////////////////////////////////////////////////////////////////////


void
SplineHysterisis4_Ctor(SplineHysterisis4 *unit) 
{
    unit->m_x1 = 0.0;
    unit->m_blend = 0.5;
    unit->m_upper = 0.0;
    unit->m_lower = 0.0;
    OVERSAMPLE4_INIT;
    SETCALC(SplineHysterisis4_next);
    ZOUT0(0) = 0.f;
}


void
SplineHysterisis4_next(SplineHysterisis4 *unit, int inNumSamples)
{
    float *in = ZIN(0);
    float *out = ZOUT(0);
    float cross_speed = ZIN0(1); // 0.78
    float cross_curve = ZIN0(2); // 0.29
    float shift = ZIN0(3); // 0.29
    float pa0 = ZIN0(4); 
    float ma0 = ZIN0(5);
    float pa1 = ZIN0(6);
    float ma1 = ZIN0(7);
    float pb0 = ZIN0(8);
    float mb0 = ZIN0(9);
    float pb1 = ZIN0(10);
    float mb1 = ZIN0(11);

    float ca0, ca1, ca2, ca3;
    float cb0, cb1, cb2, cb3;
    float norma, offa;
    float normb, offb;
    float x1 = unit->m_x1;
    float blend = unit->m_blend;
    float lower = unit->m_lower;
    float upper = unit->m_upper;
    float nsamp = inNumSamples * 4;
    float lim;

    // pa0 = 0.55; 
    // ma0 = 1.8; 
    // pa1 = 1.2;
    // ma1 = 0.25;
    // pb0 = 0.15;
    // mb0 = 1.8;
    // pb1 = 1.2;
    // mb1 = 0.25;

    spline_coeffs(&norma, &offa, &ca0, &ca1, &ca2, &ca3, pa0, ma0, pa1, ma1);
    spline_coeffs(&normb, &offb, &cb0, &cb1, &cb2, &cb3, pb0, mb0, pb1, mb1);

    lim = -(offa+offb)*0.5 + pa1;

    UPSAMPLE4;

    for (int i=0; i<nsamp; ++i) {
        float xd;
        float y;
        float x = domemoff[i];
        // xd is x-x1 (input - last input)

        xd = x - x1;
        x1 = x;

        upper = 0.9993 * upper;
        lower = 0.9993 * lower;

        if (xd > 0.f) {
            float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
            float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
            float dist = uplimit - lowlimit;
            float change = dist * cross_speed*xd/(cross_curve+blend*blend);
            blend = sc_min(change + blend, uplimit);
        } else {
            float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
            float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
            float dist = uplimit - lowlimit;
            float change = dist * cross_speed*xd/(cross_curve+(1.f-blend)*(1.f-blend));
            blend = sc_max(change + blend, lowlimit);
        }

        if (x > 0) {
            float res2, res;
            if (x > pa1) {
                y = lim + (x-pa1)*ma1;
            } else {
                if (x > pb0) {
                    float t = (x-pb0)*normb;
                    res2 = (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                } else {
                    res2 = (x*mb0 + offb);
                }
                if (x > pa0) {
                    float t = (x-pa0)*norma;
                    res = (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                } else {
                    res = (x*ma0 + offa);
                }
                y = (blend * res + (1.f-blend) * res2) - (offb + offa)*0.5f;
            }
        } else {
            float res2, res;
            if (x < -pa1) {
                y = -(lim + (-x-pa1)*ma1);
            } else {
                x = -x;
                if (x > pa0) {
                    float t = (x-pa0)*norma;
                    res2 = - (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                } else {
                    res2 = - (x*ma0 + offa);
                }
                if (x > pb0) {
                    float t = (x-pb0)*normb;
                    res = - (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                } else {
                    res = - (x*mb0 + offb);
                }
                y = (blend * res + (1.f-blend) * res2) + (offb + offa)*0.5f;
            }
        }

        if (y > upper && (y < 0.8)) {
            upper = y;
        } else if (y < lower && (y > -0.8)) {
            lower = y;
        }

        domemoff[i] = y;
    }
    
    DOWNSAMPLE4;

    unit->m_x1 = x1; 
    unit->m_lower = lower; 
    unit->m_upper = upper; 
    unit->m_blend = blend;
}



///////////////////////////////////////////////////////////////////////////

void
SplineHysterisis4L_Ctor(SplineHysterisis4L *unit) 
{
    unit->m_x1 = 0.0;
    unit->m_blend = 0.5;
    unit->m_upper = 0.0;
    unit->m_lower = 0.0;
    SETCALC(SplineHysterisis4L_next);
    ZOUT0(0) = 0.f;
    ZOUT0(1) = 0.f;
    ZOUT0(2) = 0.f;
    ZOUT0(3) = 0.f;
}


void
SplineHysterisis4L_next(SplineHysterisis4L *unit, int inNumSamples)
{
    float *in1 = ZIN(0);
    float *in2 = ZIN(1);
    float *in3 = ZIN(2);
    float *in4 = ZIN(3);
    float *out1 = ZOUT(0);
    float *out2 = ZOUT(1);
    float *out3 = ZOUT(2);
    float *out4 = ZOUT(3);
    float cross_speed = ZIN0(4); // 0.78
    float cross_curve = ZIN0(5); // 0.29
    float shift = ZIN0(6); // 0.29
    float pa0 = ZIN0(7); 
    float ma0 = ZIN0(8);
    float pa1 = ZIN0(9);
    float ma1 = ZIN0(10);
    float pb0 = ZIN0(11);
    float mb0 = ZIN0(12);
    float pb1 = ZIN0(13);
    float mb1 = ZIN0(14);
    float *ins[4] = { in1, in2, in3, in4 };
    float *outs[4] = { out1, out2, out3, out4 };

    float ca0, ca1, ca2, ca3;
    float cb0, cb1, cb2, cb3;
    float norma, offa;
    float normb, offb;
    float x1 = unit->m_x1;
    float blend = unit->m_blend;
    float lower = unit->m_lower;
    float upper = unit->m_upper;
    float nsamp = inNumSamples;
    float lim;

    spline_coeffs(&norma, &offa, &ca0, &ca1, &ca2, &ca3, pa0, ma0, pa1, ma1);
    spline_coeffs(&normb, &offb, &cb0, &cb1, &cb2, &cb3, pb0, mb0, pb1, mb1);

    lim = -(offa+offb)*0.5 + pa1;

    for (int k=0; k<4; ++k) {
        float *in = ins[k];
        float *out = outs[k];
        for (int i=0; i<nsamp; ++i) {
            float xd;
            float y;
            float x = ZXP(in);
            // xd is x-x1 (input - last input)

            xd = x - x1;
            x1 = x;

            upper = 0.9993 * upper;
            lower = 0.9993 * lower;

            if (xd > 0.f) {
                float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
                float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
                float dist = uplimit - lowlimit;
                float change = dist * cross_speed*xd/(cross_curve+blend*blend);
                blend = sc_min(change + blend, uplimit);
            } else {
                float uplimit = sc_min(1.f-shift - lower*shift*1.25, 1.f);
                float lowlimit = sc_max(shift - upper*shift*1.25, 0.f);
                float dist = uplimit - lowlimit;
                float change = dist * cross_speed*xd/(cross_curve+(1.f-blend)*(1.f-blend));
                blend = sc_max(change + blend, lowlimit);
            }

            if (x > 0) {
                float res2, res;
                if (x > pa1) {
                    y = lim + (x-pa1)*ma1;
                } else {
                    if (x > pb0) {
                        float t = (x-pb0)*normb;
                        res2 = (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                    } else {
                        res2 = (x*mb0 + offb);
                    }
                    if (x > pa0) {
                        float t = (x-pa0)*norma;
                        res = (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                    } else {
                        res = (x*ma0 + offa);
                    }
                    y = (blend * res + (1.f-blend) * res2) - (offb + offa)*0.5f;
                }
            } else {
                float res2, res;
                if (x < -pa1) {
                    y = -(lim + (-x-pa1)*ma1);
                } else {
                    x = -x;
                    if (x > pa0) {
                        float t = (x-pa0)*norma;
                        res2 = - (ca0+ca1*t+ca2*t*t+ca3*t*t*t);
                    } else {
                        res2 = - (x*ma0 + offa);
                    }
                    if (x > pb0) {
                        float t = (x-pb0)*normb;
                        res = - (cb0+cb1*t+cb2*t*t+cb3*t*t*t);
                    } else {
                        res = - (x*mb0 + offb);
                    }
                    y = (blend * res + (1.f-blend) * res2) + (offb + offa)*0.5f;
                }
            }

            if (y > upper && (y < 0.8)) {
                upper = y;
            } else if (y < lower && (y > -0.8)) {
                lower = y;
            }

            ZXP(out) = y;
        }
    }    

    unit->m_x1 = x1; 
    unit->m_lower = lower; 
    unit->m_upper = upper; 
    unit->m_blend = blend;
}


void load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineSimpleUnit(SplineHysterisis);
    DefineSimpleUnit(SplineHysterisis4);
    DefineSimpleUnit(SplineHysterisis4L);
}
