/* BitManip4.cpp
   Copyright (C) by Bjoern Erlach

   SuperCollider plugin of an interpreter for my mini audio language BitManip.
*/

#include "SC_PlugIn.h"
#include <stdio.h>


enum {
      BM_OP_NOT = 0,
      BM_OP_AND = 1,
      BM_OP_OR = 2,
      BM_OP_XOR = 3,
      BM_OP_NOR = 4,
      BM_OP_ADD = 5,
      BM_OP_SUB = 6,
      BM_OP_MUL = 7,
      BM_OP_DIV = 8,
      BM_OP_YIELD=9,
      BM_OP_ROT=10,
      BM_OP_FOR=11,
      BM_OP_DMUL=12,
      BM_OP_JMP = 14,
      BM_OP_JMPR = 15,   // check if <0 jump to end of repeat loop and decrement
      BM_OP_ASSIGN = 16,
      BM_OP_READAT = 17,
      BM_OP_BWRAP = 18,
      BM_OP_REF = 19,
      BM_OP_CALL = 20,
      BM_OP_RETURN = 21,
      BM_OP_BUFWR = 22,
      BM_OP_SHIFT = 23,
      BM_OP_BUFREAD = 24,
      BM_OP_GOTO = 25
};


enum {
      BM_TYPE_BINOP  = 0,
      BM_TYPE_UNOP   = 1,
      BM_TYPE_INT    = 2,
      BM_TYPE_FLOAT  = 3,
      BM_TYPE_STMTS  = 4,
      BM_TYPE_LIST   = 5,
      BM_TYPE_OUTPUT = 6,
      BM_TYPE_INPUT  = 7,
      BM_TYPE_INT_INPUT = 8,
      BM_TYPE_FORLOOP = 9,
      BM_OUTPUT = 10,
      BM_INTINPUT = 11,
      BM_TYPE_VAR = 12,
      BM_TYPE_ASSIGN = 13,
      BM_TYPE_STRING = 14,
      BM_TYPE_INT_BUFFER = 15,
      BM_TYPE_FLOAT_BUFFER = 16
};


// cubic hermite interpolator 
static inline double cubic_hermite (double t, double x0, double x1, double m0, double m1)
{
    double t2 = t*t;
    return (3-2*t)*t2*x1+(t2*(2*t-3)+1)*x0+t*(t*((m1+m0)*t-m1-2*m0)+m0);
}


// kaiser windowed sinc halfband 19 point kaiser param 7.4
// needs 28 points of memory
static inline float decimate_mp4 (float *m)
{
    float stage1p1 =  m[9] * 0.5000459243617790816571755f + 
	(m[10] + m[8]) * 0.3050653325548928762600553f+
	(m[12] + m[6]) * -0.0716613428162717197578147f+
	(m[14] + m[4]) * 0.0201463119524870268306493f+
	(m[16] + m[2]) * -0.0037180247079463052059523f+
	(m[18] + m[0]) * 0.0001447608359485727503278f;
    float stage1p2 =  m[11] * 0.5000459243617790816571755f + 
	(m[12] + m[10]) * 0.3050653325548928762600553f+
	(m[14] + m[8]) * -0.0716613428162717197578147f+
	(m[16] + m[6]) * 0.0201463119524870268306493f+
	(m[18] + m[4]) * -0.0037180247079463052059523f+
	(m[20] + m[2]) * 0.0001447608359485727503278f;
    memmove(m, m+4, sizeof(float) * 17);
    float tmp = 0.405251153879167103078f * m[21] + 0.01785905009866808196283827f * (stage1p1 + m[22]);	
    float as1 = m[22] + tmp;  
    m[22] = m[21];  m[21] = stage1p1 - tmp; 
    tmp = 1.558304379470869971414f * m[23] + 0.59601645601411934460855946f * (as1 + m[24]);	
    float a1 = m[24] + tmp;  
    m[24] = m[23];  m[23] = as1 - tmp; 
    tmp = 0.709185169532354664490f * m[25] + 0.09619761143814288539832802f * (m[26] + m[27]);
    float as2 = m[26] + tmp;  
    m[26] = m[25];  m[25] = m[27] - tmp; m[27] = stage1p2;
    tmp = 1.754039233447663992393f * m[28] + 0.76186483819836081643472880f * (as2 + m[29]);	
    float a2 = m[29] + tmp;  
    m[29] = m[28];  m[28] = as2 - tmp; 
    return a1+a2;
}

static inline void write_sample_mp4 (float *m, int pp, float in)
{
    m[pp+17] = in;
}


/** Before you use this make sure that int32_t and float have the same size */
static inline float store_int_to_float (int32_t a)
{
    union { int i; float f; } ca;
    ca.i = a;
    return ca.f;
}


/** Before you use this make sure that int32_t and float have the same size */
static inline int32_t store_float_to_int (float a)
{
    union { int i; float f; } ca;
    ca.f = a;
    return ca.i;
}


#define BM_MAGIC 0x00be2013
#define INTBITS(x) (x)
#define INTSHIFT(x) (x-1)
#define BITMASK(b) ((1<<INTBITS(b))-1)
#define SCALE(x,b) (((float) BITMASK(x))/(float)-(1<<INTSHIFT(b)))
#define TOINT(x,b) ((int)(x * (float)(1<<INTSHIFT(b))))


static inline int bitrotr (int32 x, int32 n, int32 b)
{
    unsigned int mask = (1<<(INTBITS(b)-n)) -1;
    unsigned int high = (~mask) & x;
    return (x << n) + (high >> n);
}


static inline int bitrotl (int32 x, int32 y, int32 b)
{
    unsigned int mask = (1 << y) - 1;
    unsigned int low = mask & x;
    return (x >> y) + (low << (INTBITS(b)-y));
}


static inline int bitrot (int32 x, int32 n, int32 b)
{
    if (n > 0) {
        return bitrotl(x, n, b);
    } else if (n < 0) {
        return bitrotr(x, -n, b);
    }
    return x;
}



static InterfaceTable *ft;


struct NonUniConverter {
    double tin;
    double tout;
    int outidx;
    int inidx;
    int delay;
    float input_mem[256];
    double outx[512];
    double outt[512];
    int idx;
    double d1, d2, d3;
    double m1, m2;
    double x0, x1, x2, x3;
    double t0, t1, t2, t3, t4;
    double playt;
    double integrator1;
    float outmem[32];
};


struct IOShifts {
    int start;
    int end;
};


struct BitManip4 : public Unit
{
    NonUniConverter conv;
    int32 *m_bytecode;
    int32 *m_mem;
    int ip;
    uint bd;
    int m_bclen;
    double m_time;
    float *inputmem[8];
    float *int_inputmem[8];
    float *outputmem[8];
    int m_inputs[8];       // mem locations for io
    int m_int_inputs[8];
    int m_outputs[8];
    int m_num_inputs;
    int m_num_int_inputs;
    int m_num_outputs;
    IOShifts m_shifts[16];
    int m_num_shifts;
    double m_periodend;
    int m_idx;
    int m_idx2;
    double m_samptime;
    float m_prevout;
    float m_prevout2;
    int m_num_buffers;
    int m_labels[16];
    int32_t *m_buffers[16];
    int32_t m_buffer_sizes[16];
};



static inline void init_nonuni (NonUniConverter *conv)
{
    memset(conv, 0, sizeof(NonUniConverter));    
    conv->delay = 48;
    conv->d1 = 1;
    conv->d2 = 1;
    conv->d3 = 1;
    conv->outt[0] = 0;
    conv->outt[1] = 1;
    conv->outt[2] = 2;
    conv->outt[3] = 3;
    conv->outt[4] = 4;
    conv->outt[5] = 5;
    conv->outt[6] = 6;
    conv->outt[7] = 7;
    conv->outt[8] = 8;
    conv->outt[9] = 9;
    conv->outt[10] = 10;
    conv->outt[11] = 11;
    conv->outt[12] = 12;
    conv->outt[13] = 13;
    conv->outt[14] = 14;
}


// (outidx-1) points to the last sample taken

static inline void resample_nonuni_lin (NonUniConverter *conv, float *ratein, float *in,  int inNumSamples)
{
    double tin = conv->tin;
    double tout = conv->tout;
    double *outx = conv->outx;
    double *outt = conv->outt;
    int outidx = conv->outidx;
    int delay = conv->delay;
    float *input_mem = conv->input_mem;
    int inidx = conv->inidx;
    int freeAtEnd = 256 - inidx;
    if (freeAtEnd < inNumSamples) {
	memcpy(input_mem+inidx, in, sizeof(float) * freeAtEnd);
	memcpy(input_mem, in+freeAtEnd, sizeof(float) * (inNumSamples - freeAtEnd));
    } else {
	memcpy(input_mem+inidx, in, sizeof(float) * inNumSamples);
    }
    conv->inidx = (inidx + inNumSamples) & 255;

    // determine the time instances of the internal sample clock
    for (int k=0; k<inNumSamples; ++k) {
	float rate = ZXP(ratein);
	tin += 1.f;                      // this is just k but as a float - avoid casting
	while (tout < tin) {
	    outt[outidx & 511] = tout;
	    int it = (int) tout;
	    float frac = tout - sc_floor(tout);
	    double val = (1.f-frac) * input_mem[(it-delay-1)&255] + frac * input_mem[(it-delay)&255];
	    outx[outidx & 511] = val;
	    outidx++;
	    tout += rate;
	}      
    }
    conv->tin = tin;
    conv->tout = tout;
    conv->outidx = outidx;    
}


static inline void output_nonuni (NonUniConverter *conv, float *out,  int inNumSamples)
{
    double integrator1 = conv->integrator1;
    double playt = conv->playt;
    int idx = conv->idx;
    double d1 = conv->d1;
    double d2 = conv->d2;
    double d3 = conv->d3;
    double x0 = conv->x0;
    double x1 = conv->x1;
    double x2 = conv->x2;
    double x3 = conv->x3;
    double m2 = conv->m2;
    double m1 = conv->m1;
    double t0 = conv->t0;
    double t1 = conv->t1;
    double t2 = conv->t2;
    double t3 = conv->t3;
    double t4 = conv->t4;
    double *outt = conv->outt;
    double *outx = conv->outx;
    float *mem = conv->outmem;

    for (int k=0; k<inNumSamples; ++k) {
	// find the interval so that t1 < t <= t2
	while ((playt-12.f) > t2) {
	    idx++;
	    double nextt = outt[(idx+1)&511];
	    if ((nextt-t3) < 0.0000001) continue;
	    t0 = t1;  t1 = t2; t2 = t3; t3 = t4; t4 = nextt;
	    d1 = d2; d2 = d3; d3 = t4-t3;
	    x0 = x1; x1 = x2; x2 = x3; 
	    x3 = outx[idx&511]; 
	    x3 = integrator1 = 0.9999*integrator1 + x3 * d3; 
	    double dr1 = 1.0/d1;
	    double dr2 = 1.0/d2;
	    m1 = m2; m2 = 0.5*(dr2*x3+(dr1-dr2)*x2-dr1*x1);
	}
	double mm1 = d1 * m1;
	double mm2 = d1 * m2;
	float val;
	if (d1 > 0.000001) {
	    val = (float) cubic_hermite( ((playt-12.f)-t1)/d1 , x1, x2, mm1, mm2);
	} else {
	    val = 0.0;
	}
	write_sample_mp4(mem, 0, val);
	playt += 0.25;

	// find the interval so that t1 < t <= t2
	while ((playt-12.f) > t2) {
	    idx++;
	    double nextt = outt[(idx+1)&511];
	    if ((nextt-t3) < 0.0000001) continue;
	    t0 = t1;  t1 = t2; t2 = t3; t3 = t4; t4 = nextt;
	    d1 = d2; d2 = d3; d3 = t4-t3;
	    x0 = x1; x1 = x2; x2 = x3; 
	    x3 = outx[idx&511]; 
	    x3 = integrator1 = 0.9999*integrator1 + x3 * d3; 
	    double dr1 = 1.0/d1;
	    double dr2 = 1.0/d2;
	    m1 = m2; m2 = 0.5*(dr2*x3+(dr1-dr2)*x2-dr1*x1);
	}
	mm1 = d1 * m1;
	mm2 = d1 * m2;
	//fprintf(stderr, "d1: %f\n", d1);
	if (d1 > 0.000001) {
	    val = (float) cubic_hermite( ((playt-12.f)-t1)/d1 , x1, x2, mm1, mm2);
	} else {
	    val = 0.0;
	}
	write_sample_mp4(mem, 1, val);
	playt += 0.25;

	// find the interval so that t1 < t <= t2
	while ((playt-12.f) > t2) {
	    idx++;
	    double nextt = outt[(idx+1)&511];
	    if ((nextt-t3) < 0.0000001) continue;
	    t0 = t1;  t1 = t2; t2 = t3; t3 = t4; t4 = nextt;
	    d1 = d2; d2 = d3; d3 = t4-t3;
	    x0 = x1; x1 = x2; x2 = x3; 
	    x3 = outx[idx&511]; 
	    x3 = integrator1 = 0.9999*integrator1 + x3 * d3; 
	    double dr1 = 1.0/d1;
	    double dr2 = 1.0/d2;
	    m1 = m2; m2 = 0.5*(dr2*x3+(dr1-dr2)*x2-dr1*x1);
	}

	mm1 = d1 * m1;
	mm2 = d1 * m2;

	if (d1 > 0.000001) {
	    val = (float) cubic_hermite( ((playt-12.f)-t1)/d1 , x1, x2, mm1, mm2);
	} else {
	    val = 0.0;
	}
	write_sample_mp4(mem, 2, val);
	playt += 0.25;

	while ((playt-12.f) > t2) {
	    idx++;
	    double nextt = outt[(idx+1)&511];
	    if ((nextt-t3) < 0.0000001) continue;
	    t0 = t1;  t1 = t2; t2 = t3; t3 = t4; t4 = nextt;
	    d1 = d2; d2 = d3; d3 = t4-t3;
	    x0 = x1; x1 = x2; x2 = x3; 
	    x3 = outx[idx&511]; 
	    x3 = integrator1 = 0.9999*integrator1 + x3 * d3; 
	    double dr1 = 1.0/d1;
	    double dr2 = 1.0/d2;
	    //fprintf(stderr, "%f %f\n", dr1, dr2);
	    m1 = m2; m2 = 0.5*(dr2*x3+(dr1-dr2)*x2-dr1*x1);
	}

	mm1 = d1 * m1;
	mm2 = d1 * m2;

	//fprintf(stderr, "d1: %f\n", d1);
	if (d1 > 0.000001) {
	    val = (float) cubic_hermite( ((playt-12.f)-t1)/d1 , x1, x2, mm1, mm2);
	} else {
	    val = 0.0;
	}

	write_sample_mp4(mem, 3, val);
	playt += 0.25;
	val = decimate_mp4(mem);
	ZXP(out) = val;
    }

    conv->idx = idx;
    conv->d1 = d1;
    conv->d2 = d2;
    conv->d3 = d3;
    conv->x0 = x0;
    conv->x1 = x1;
    conv->x2 = x2;
    conv->x3 = x3;
    conv->m2 = m2;
    conv->m1 = m1;
    conv->t0 = t0;
    conv->t1 = t1;
    conv->t2 = t2;
    conv->t3 = t3;
    conv->t4 = t4;
    conv->playt = playt;
    conv->integrator1 = integrator1;
}


extern "C"
{
    int api_version(void);
    void load(InterfaceTable *inTable);
    void BitManip4_Dtor(BitManip4 *unit);
    void BitManip4_Ctor(BitManip4 *unit);
    void BitManip4_next(BitManip4 *unit, int inNumSamples);
}

int api_version(void) 
{ 
    return sc_api_version; 
}





void BitManip4_Dtor(BitManip4 *unit)
{
    fprintf(stderr, "dtor a\n");
    for (int n=0; n<unit->m_num_inputs; n++) {
	RTFree(unit->mWorld, unit->inputmem[n]);
    }
    for (int n=0; n<unit->m_num_int_inputs; n++) {
	RTFree(unit->mWorld, unit->int_inputmem[n]);
    }
    fprintf(stderr, "dtor b\n");
    for (int n=0; n<unit->m_num_outputs; n++) {
	RTFree(unit->mWorld, unit->outputmem[n]);
    }
    for (int n=0; n<unit->m_num_buffers; n++) {
	RTFree(unit->mWorld, unit->m_buffers[n]);
    }
    fprintf(stderr, "dtor c\n");
    RTFree(unit->mWorld, unit->m_mem);
    fprintf(stderr, "dtor d\n");
    //RTFree(unit->mWorld, unit->m_bytecode);
}

#define UINTS(x) (*(uint32_t*)x)

void BitManip4_Ctor(BitManip4 *unit)
{
    init_nonuni(&unit->conv);
    SETCALC(BitManip4_next);
    float fbufnum  = ZIN0(0); 
    int bd  = (int) ZIN0(2); 	
    uint32 bufnum = (int)fbufnum; 
    World *world = unit->mWorld; 
    if (bufnum >= world->mNumSndBufs) bufnum = 0; 
    SndBuf *buf = world->mSndBufs + bufnum ;
    float* ptr = buf->data;
    //int bufFrames = buf->frames;
    int val;
    int mem_begin;
    int bclen;
    int required_memory;
    unit->bd = bd;
    unit->m_prevout = 0.f;
    unit->m_prevout2 = 0.f;
    unit->m_samptime = 0.f;

    val = store_float_to_int(*ptr);
    if(val != BM_MAGIC) {
        fprintf(stderr, "BitManip4: not a bitmanip buffer\n");
    }
    ptr++;
    val = store_float_to_int(*ptr);
    required_memory = val;
    char *all_memory = (char*) RTAlloc(unit->mWorld, val*8+256);
    memset(all_memory, 0, val*8+256);
    ptr++;
    fprintf(stderr, "a : %x\n", UINTS(ptr));
    val = store_float_to_int(*ptr);
    mem_begin = val;
    unit->m_mem = (int32_t*) all_memory;
    int32_t *mem = (int32_t*) all_memory;

    ptr++;
    fprintf(stderr, "b: %x\n", UINTS(ptr));
    ptr++;
    // read constant and variables initializations
    while ( (*(uint32_t*)ptr) != 0xeeeeeeee) {
        val = store_float_to_int(*ptr);
        if ((uint32_t)val == 0xeeeeeeee) {
            fprintf(stderr, "broken bitmanip file\n");
        }
        ptr++;
        if (val == BM_TYPE_INT) {
            int loc = store_float_to_int(*ptr);
            ptr++;
            val = store_float_to_int(*ptr);
            ptr++;
            fprintf(stdout, "iconstant: %d\n", val);
            mem[loc] = val;
        } else {
            // convert number to fixed point representation
            float fval;
            int loc = store_float_to_int(*ptr);
            ptr++;
            fval = *ptr;
            ptr++;
            mem[loc] = (int32) (fval * (1<<bd));
            fprintf(stdout, "constant: %f\n", fval);
        }
    }

    fprintf(stderr, "c : %x\n", UINTS(ptr));
    ptr++;

    unit->m_num_inputs = store_float_to_int(*ptr);
    fprintf(stderr, "num inputs: %d\n", unit->m_num_inputs);
    ptr++;
    for (int k=0; k<unit->m_num_inputs; k++) {
        unit->m_inputs[k] = store_float_to_int(*ptr);
        unit->inputmem[k] = (float*) RTAlloc(unit->mWorld, 256 *  sizeof(float));
        memset(unit->inputmem[k], 0, sizeof(float) * 256);
        ptr++;
    }

    fprintf(stderr, "d : %x\n", UINTS(ptr));
    ptr++;
    unit->m_num_int_inputs = store_float_to_int(*ptr);
    fprintf(stderr, "num int inputs: %d\n", unit->m_num_int_inputs);
    ptr++;
    for (int k=0; k<unit->m_num_int_inputs; k++) {
        unit->m_int_inputs[k] = store_float_to_int(*ptr);
        unit->int_inputmem[k] = (float*) RTAlloc(unit->mWorld, 256 *  sizeof(float));
        memset(unit->int_inputmem[k], 0, sizeof(float) * 256);
        ptr++;
    }
    fprintf(stderr, "e : %x\n", UINTS(ptr));
    ptr++;
    unit->m_num_outputs = store_float_to_int(*ptr);
    fprintf(stderr, "num outputs: %d\n", unit->m_num_outputs);
    ptr++;
    for (int k=0; k<unit->m_num_outputs; k++) {
        unit->m_outputs[k] = store_float_to_int(*ptr);
        unit->outputmem[k] = (float*) RTAlloc(unit->mWorld, 256 *  sizeof(float));
        memset(unit->outputmem[k], 0, sizeof(float) * 256);
        ptr++;
    }
    fprintf(stderr, "f : %x\n", UINTS(ptr));
    ptr++;
    int i = 0;
    unit->m_num_shifts = store_float_to_int(*ptr);
    ptr++;
    fprintf(stderr, "num shifts: %d\n", unit->m_num_shifts);
    // read the shifts
    while ( (*(uint32_t*)ptr) != 0xeeeeeeee) {
        // start end and num
        unit->m_shifts[i].start = store_float_to_int(*ptr);
        fprintf(stderr, "shift: %d ", unit->m_shifts[i].start);
        ptr++;
        unit->m_shifts[i].end = store_float_to_int(*ptr);
        fprintf(stderr, "%d ", unit->m_shifts[i].end);
        ptr++;
        i++;
    }
    unit->m_num_shifts = i;
    fprintf(stderr, "g : %d, %x\n", i, UINTS(ptr));	
    ptr++;
    unit->m_num_buffers = store_float_to_int(*ptr);
    fprintf(stderr, "num buffers: %d\n", unit->m_num_buffers);	
    ptr++;
    i = 0;
    while ( (*(uint32_t*)ptr) != 0xeeeeeeee) {
        // start end and num
        int bsize = store_float_to_int(*ptr);
        fprintf(stderr, "buffer size: %d\n", bsize);	
        unit->m_buffers[i] = (int32_t*) RTAlloc(unit->mWorld, sizeof(int32_t) * bsize);
        unit->m_buffer_sizes[i] = bsize;
        memset(unit->m_buffers[i], 0, sizeof(int32_t) * bsize);
        ptr++;
        i++;
    }
    fprintf(stderr, "h : %d, %x\n", i, UINTS(ptr));	
    ptr++;
    // labels
    int numlabels = store_float_to_int(*ptr);
    ptr++;
    for (int i=0; i<numlabels; i++) {
        unit->m_labels[i] = store_float_to_int(*ptr);
        fprintf(stderr, "label: %d\n", unit->m_labels[i]);
        ptr++;
    }
    fprintf(stderr, "h : %d, %x\n", i, UINTS(ptr));	
	
    ptr++;
    bclen = store_float_to_int(*ptr);
    fprintf(stderr, "bclen: %d\n", bclen);
    ptr++;
    unit->m_bytecode = (int32_t*) (all_memory + mem_begin * sizeof(int32_t));
    //	fprintf(stderr, "bufFrames: %d\n", bufFrames);
    for (int i=0; i<bclen; i++) {
        unit->m_bytecode[i] = store_float_to_int(*ptr++);
    }
	
    fprintf(stderr, "%ld %d %d\n", ptr-buf->data, required_memory, unit->m_bytecode[0]);

    unit->m_bclen = bclen;
    unit->ip = 0;
    unit->m_time = 0.f;	
    unit->m_periodend = 0.f;
    unit->m_idx = 0;
    unit->m_idx2 = 0;
    ZOUT0(0) = 0.f;
}



static inline float read_lininterp (float *buf, float time)
{
    int itime = (int) time;
    float frac = time - (float)itime;
    return (1.f-frac) * buf[itime] + frac * buf[itime+1];
}


static inline float read_cubic (float *ptr, float pos)
{
    int ipos = (int) pos;
    float frac = pos - (float) pos;
    return cubicinterp(frac, ptr[ipos], ptr[(ipos+1)] , ptr[(ipos+2)] , ptr[(ipos+3)] );
}


static inline float read_cubic_circ(float *ptr, double pos)
{
    int ipos = (int) pos;
    float frac = (float) (pos - (double) ipos);
    return cubicinterp(frac, ptr[ipos&255], ptr[(ipos+1)&255] , ptr[(ipos+2)&255] , ptr[(ipos+3)&255] );
}

//#define PRINTOP(OP) fprintf(stderr, OP)
#define PRINTOP(OP) (NULL)


// inputs are bitdepth, rate, [inputs]

void BitManip4_next(BitManip4 *unit, int inNumSamples)
{
    float *out = OUT(0);
    int numinputs = unit->mNumInputs;
    int32 *bc = unit->m_bytecode;
    int ip = unit->ip;                  // instruction pointer
    int bclen = unit->m_bclen;
    int32 *mem = unit->m_mem;           // memory for variables, inputs and outputs
    int32 res;
    uint32 bd = unit->bd;
    int numin = unit->m_num_inputs;
    int numintin = unit->m_num_int_inputs;
    int numshifts = unit->m_num_shifts;
    double time = unit->m_time;
    float ts = ZIN0(1);
    float **inputmem = unit->inputmem;
    float **int_inputmem = unit->int_inputmem;
    float **outputmem = unit->outputmem;
    double periodend = unit->m_periodend+64.0;
    int idx = unit->m_idx;
    int *inputs = unit->m_inputs;
    int *int_inputs = unit->m_int_inputs;
    IOShifts *shifts = unit->m_shifts;
    int32_t retipstack[16];  // subroutine return address stack
    int32_t retipcur = 0;
    int idx2 = unit->m_idx2;
    double samptime = unit->m_samptime;
    float prevout = unit->m_prevout;
    float prevout2 = unit->m_prevout2;
    float starttime = time;
    int32_t **buffers = unit->m_buffers;
    int32_t *buffersizes = unit->m_buffer_sizes;

    if ((numin+numintin) != numinputs-3) {
	fprintf(stderr, "mismatch numinputs %d %d\n", numin, numinputs-3);
    }

    for (int n=0; n<numin; n++) {
	float *in = IN(n+3);

        // the next line was missing?!
        idx2 = unit->m_idx2;
	for (int k=0; k<inNumSamples; ++k) {
	    inputmem[n][idx2&255] = *in++;
	    idx2 ++;
	}
    }

    for (int n=0; n<numintin; n++) {
	float *in = IN(n+3+numin);
        idx2 = unit->m_idx2;
	for (int k=0; k<inNumSamples; ++k) {
	    int_inputmem[n][idx2&255] = *in++;
	    idx2 ++;
	}
    }

    //fprintf(stderr, "loop\n");
    while(time<periodend) {
	int op;

	// load the input cells
	for (int n=0; n<numin; n++) {
	    mem[inputs[n]] = (int32) (read_cubic_circ(inputmem[n], time-8.0) * ((float) (1<<bd)));
	}
	for (int n=0; n<numintin; n++) {
	    //mem[int_inputs[n]] = (int32) (read_cubic_circ(int_inputmem[n], time-8.0));
            mem[int_inputs[n]] = (int32) int_inputmem[n][(((int)time)-8)&255];
	}

	for(;;) {
	    if (ip >= bclen) ip=0;
	    op = bc[ip];
            
	    //fprintf(stderr, "switch: ip:%d m5:%d m4:%d m3:%d\n", ip, mem[5], mem[4], mem[3]);
	    switch (op) {

	    case BM_OP_CALL:
		retipstack[retipcur] = ip + 2; 
		retipcur++;
		ip = bc[ip+1];		
		break;

	    case BM_OP_RETURN:
		retipcur--;
		ip = retipstack[retipcur]; 
		break;

	    case BM_OP_NOT:
                PRINTOP("NOT");
		mem[bc[ip+1]] = ~mem[bc[ip+2]];
		ip += 3;
		break;

	    case BM_OP_REF: {
                int bn;
                PRINTOP("REF");
                // buffersize if found from buffer table
                bn = bc[ip+2];
                mem[bc[ip+1]] = buffers[bn] [(idx-(uint32_t) mem[bc[ip+3]]) % buffersizes[bn]];
		ip += 4;
		break;
            }
		
	    case BM_OP_BUFWR:
                PRINTOP("BUFWR");
                //fprintf(stderr, "w: %d %d %d %d\n", bc[ip+1], bc[ip+2], bc[ip+3], mem[bc[ip+3]]);
		buffers[bc[ip+1]][idx%bc[ip+2]] = mem[bc[ip+3]];
		ip += 4;
		break;

	    case BM_OP_AND:
                PRINTOP("AND");
		mem[bc[ip+1]] = mem[bc[ip+2]]&mem[bc[ip+3]];
		ip += 4;
		break;
		
	    case BM_OP_XOR:
                PRINTOP("XOR");
		mem[bc[ip+1]] = mem[bc[ip+2]]^mem[bc[ip+3]];
		ip += 4;
		break;
		
	    case BM_OP_NOR:
		mem[bc[ip+1]] = ~(mem[bc[ip+2]]|mem[bc[ip+3]]);
		ip += 4;
		break;
		
	    case BM_OP_OR:
                PRINTOP("OR");
		mem[bc[ip+1]] = mem[bc[ip+2]]|mem[bc[ip+3]];
		ip += 4;
		break;

	    case BM_OP_ASSIGN:
                PRINTOP("ASSIGN\n");
		mem[bc[ip+1]] = mem[bc[ip+2]];
		ip += 3;
		break;
		
	    case BM_OP_ADD:
                PRINTOP("ADD");
		//fprintf(stderr, "add: %d %d %d\n", bc[ip+1], bc[ip+2], bc[ip+3]);
		mem[bc[ip+1]] = mem[bc[ip+2]]+mem[bc[ip+3]];
		ip += 4;
		break;

	    case BM_OP_SUB:
                PRINTOP("ADD");
		//fprintf(stderr, "add: %d %d %d\n", bc[ip+1], bc[ip+2], bc[ip+3]);
		mem[bc[ip+1]] = mem[bc[ip+2]]-mem[bc[ip+3]];
		ip += 4;
		break;

	    case BM_OP_ROT:
		mem[bc[ip+1]] = bitrot(mem[bc[ip+2]],mem[bc[ip+3]],bd);
		ip += 4;
		break;

	    case BM_OP_SHIFT:  {
		int sh = mem[bc[ip+3]];
		if (sh < 0) { 
		    mem[bc[ip+1]] = mem[bc[ip+2]] << sh;
		} else {
		    mem[bc[ip+1]] = mem[bc[ip+2]] >> sh;
		}
		ip += 4;
            }
		break;
		
	    case BM_OP_MUL:
                PRINTOP("MUL");
		//fprintf(stderr, "mul %f\n", ((float) mem[bc[ip+2]])/(1<<15) );
		mem[bc[ip+1]] = ((mem[bc[ip+2]]>>1) * (mem[bc[ip+3]]>>1)) >> (bd-2);
		ip += 4;
		break;

	    case BM_OP_JMP:
                PRINTOP("JMP");		
		ip = bc[ip+1];
		break;

	    case BM_OP_BWRAP: {
		int32_t a = mem[bc[ip+2]];
		int32_t r;
		int const m = 1U << (bd - 1); // mask can be pre-computed if b is fixed
		a = a & ((1U << bd) - 1);  // (Skip this if bits in x above position b are already zero.)
		r = (a ^ m) - m;
		mem[bc[ip+1]] = r;
	        ip += 3;
            }
		break;

	    case BM_OP_JMPR:
                //fprintf(stderr, "JMPR %d\n", mem[bc[ip+1]]);
		// decrement counter value
		if (mem[bc[ip+1]] <= 0) {
		    ip = bc[ip+2];
		} else {
		    --mem[bc[ip+1]];
		    ip += 3;
		}
		break;

	    case BM_OP_GOTO:
		ip = unit->m_labels[(int32_t)mem[bc[ip+1]]];
		break;

	    case BM_OP_BUFREAD:
                PRINTOP("BREAD");
		mem[bc[ip+1]] = buffers[bc[ip+2]] [(idx-(uint32_t) bc[ip+3]) % bc[ip+4]];
		ip += 5;
		//BMBuffer *bmbuf = bmbufs[bc[ip+2]];
		//mem[bc[ip+1]] = bmbuf->data[(idx-mem[bc[ip+3]]) % bmbuf->length];
		break;

	    case BM_OP_YIELD:
                PRINTOP("YIELD");
		mem[bc[ip+1]] = mem[bc[ip+2]];
		res = mem[bc[ip+1]];
		ip += 3;
		goto yield_out;
		break;
		
	    default:
		fprintf(stderr, "something went wrong!\n");
		goto yield_out;
	    }	
	}
    yield_out:
	// output is in res

	// shift input and outputs
	mem[0] = res;

	for (int i=0; i<numshifts; i++) {
	    for (int n=shifts[i].end-1; n>=shifts[i].start; n--) {
                //fprintf(stderr, "%d\n", n);
		mem[n+1] = mem[n];
	    }
	}

	float outval = ((float) (res&0x80000000?(-1<<bd)+(res & ((1<<bd)-1)):(res & ((1<<bd)-1)))) / (1<<bd);	
        prevout = (0.95f * prevout + 0.05f * outval);
	outputmem[0][idx&255] = prevout;
        time = time + ts;
	idx++;
    }

    samptime = unit->m_samptime;

    // do the output

    idx2 = unit->m_idx2;
    for (int k=0; k<inNumSamples; k++) {
	float x= read_cubic_circ(outputmem[0], samptime-14.f) * 10.f ;
	*out++ = x - 0.94*prevout2;
	prevout2 = x;
	samptime += 1.f/ts;
	idx2++;
    }

    unit->m_prevout = prevout;
    unit->m_prevout2 = prevout2;
    unit->m_idx = idx;
    unit->m_time = time;
    unit->ip = ip;
    unit->m_idx2 = idx2;
    unit->m_periodend = periodend;
    unit->m_samptime = samptime;
}
    

void
load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineDtorCantAliasUnit(BitManip4);
}




