/* LowLatConv.cpp
   Copyright (C) 2020 by Bjoern Erlach

   Fast low latency fft convolver.
*/

#include <SC_PlugIn.h>

#include <string.h>
#include <stdlib.h>
#include <fftw3.h>

static InterfaceTable *ft;

#define xmemcpy memcpy

#define MAXI_BUFS 4
#define MASK 3

#define MAXI_BUFS512 8
#define MAXI_BUFS4096 256

#define is_aligned(POINTER, BYTE_COUNT) ((((uintptr_t)(const void *)(POINTER)) & (BYTE_COUNT-1)) == 0)
#define align_up(num, align) (((num) + ((align) - 1)) & ~((align) - 1))

#define FFTSIZE1 128
#define FFTSIZE1BY2 64


struct LowLatConv : public Unit
{
    float *fdbufs[MAXI_BUFS];
    float *fdbufmem;
    float *fdbufs512[MAXI_BUFS512];
    float *fdbufmem512;
    float *fdbufs4096[MAXI_BUFS4096];
    float *fdbufmem4096;
    float *m_accummem;
    float *m_accum;
    float *m_inbuf;
    float *m_accum512;
    float *m_accumother512;
    float *m_accum4096;
    float *m_accumother4096;
    int m_accidx512;
    int m_accidx4096;
    float *fd_irs[MAXI_BUFS];
    float *fd_irs512[MAXI_BUFS512];
    float *fd_irs4096[MAXI_BUFS4096];
    int m_4096irs;
    int m_512irs;
    float *m_read512;
    float *m_read4096;
    int inbufs_collected;
    int inbufs_per_fdbuf;
    int inbufs_per_fdbuf512;
    int inbufs_per_fdbuf4096;
    int num_fdbufs;
    int fd_idx;
    int fd_idx512;
    int fd_idx4096;
    int m_bypass;
    int m_4096mask;
    fftwf_plan fplan;
    fftwf_plan bplan;
    fftwf_plan fplan512;
    fftwf_plan bplan512;
    fftwf_plan fplan4096;
    fftwf_plan bplan4096;
};



extern "C"  {
int api_version(void);
void load(InterfaceTable *inTable);

void LowLatConv_Ctor(LowLatConv *unit);
void LowLatConv_next(LowLatConv *unit, int inNumSamples);
void LowLatConv_Dtor(LowLatConv *unit);
}

int api_version(void) 
{ 
    return sc_api_version; 
}


// It is important that the given buffer is exactly the right size (this has to be ensured from the language side)
void LowLatConv_Dtor(LowLatConv *unit)
{
    World *world = unit->mWorld;
    //float *fdbufsmem = (float*) RTAlloc(unit->mWorld, (fftsize+4) * unit->num_fdbufs * sizeof(float));    
    if (unit->fplan) {
	RTFree(world, unit->fdbufmem);
	RTFree(world, unit->fdbufmem512);
	RTFree(world, unit->fdbufmem4096);
	RTFree(world, unit->m_accummem);
	RTFree(world, unit->m_inbuf);
	fftwf_destroy_plan(unit->fplan);
	fftwf_destroy_plan(unit->bplan);
	fftwf_destroy_plan(unit->fplan512);
	fftwf_destroy_plan(unit->bplan512);
	fftwf_destroy_plan(unit->fplan4096);
	fftwf_destroy_plan(unit->bplan4096);
    }
}

static inline int32 nextpow2 (int32 a)
{
    int res = 4;
    while (res < a) {
        res <<= 1;
    }
    return res;
}


void LowLatConv_init(LowLatConv *unit)
{
    // get the buffer
    uint32 bufnum = (uint32)ZIN0(1);
    SndBuf *buf;
    unit->inbufs_collected = 0;
    unit->fd_idx = 0;
    unit->fd_idx512 = 0;
    unit->fd_idx4096 = 0;
    unit->inbufs_per_fdbuf = 64 / BUFLENGTH;
    unit->inbufs_per_fdbuf512 = 256 / BUFLENGTH;
    unit->inbufs_per_fdbuf4096 = 2048 / BUFLENGTH;
    fprintf(stderr, "numbufs: %d\n", unit->inbufs_per_fdbuf);

    buf = unit->mWorld->mSndBufs + bufnum;
    if (!buf->data) {
	printf("BPartConv Error: Spectral data buffer not allocated \n");
	unit->fplan = (fftwf_plan) 0;
	SETCALC(*ClearUnitOutputs);
	unit->mDone = true;
	return;
    }

    int frames = buf->frames;
    int num4096 = 1;
    int numirs4096 = 0;
    int numirs512 = 7;
    if (frames > 3968) {
        numirs4096 = (frames-3968) / 4096;
        //if (frames > (numirs4096*4096-3968)) {
        //    numirs4096 ++;
        // }
        num4096 = nextpow2(numirs4096);  
    } else {
        numirs512 = (frames-128*3)/512;
    }
    unit->m_4096irs = numirs4096;
    unit->m_4096mask = num4096-1;
    unit->m_512irs = numirs512;
    fprintf(stderr, "num4096: %d %d %d\n", numirs4096, num4096, numirs512);

    unit->fdbufmem = (float*) RTAlloc(unit->mWorld, sizeof(float)*4*128);
    memset(unit->fdbufmem, 0, sizeof(float)*128*4);
    unit->fdbufmem512 = (float*) RTAlloc(unit->mWorld, sizeof(float)*8*512);
    memset(unit->fdbufmem512, 0, sizeof(float)*8*512);
    unit->fdbufmem4096 = (float*) RTAlloc(unit->mWorld, sizeof(float)*num4096*4096);
    memset(unit->fdbufmem4096, 0, sizeof(float)*num4096*4096);

    for (int k=0; k<4; ++k) {
        unit->fdbufs[k] = unit->fdbufmem+(FFTSIZE1*k);
	if (!is_aligned(unit->fdbufs[k],16)) {
	    fprintf(stderr, "fdbuf not aligned\n");
	}
    }

    for (int k=0; k<8; ++k) {
        unit->fdbufs512[k] = unit->fdbufmem512+(512*k);
        if (!is_aligned(unit->fdbufs512[k],16)) {
            fprintf(stderr, "fdbuf not aligned\n");
        }
    }

    for (int k=0; k<num4096; ++k) {
        unit->fdbufs4096[k] = unit->fdbufmem4096+(4096*k);
        if (!is_aligned(unit->fdbufs4096[k],16)) {
            fprintf(stderr, "fdbuf not aligned\n");
        }
    }

    float *data = buf->data;
    for (int k=0; k<3; ++k) {
	unit->fd_irs[k] = (data+(FFTSIZE1)*k);
	if (!is_aligned(unit->fd_irs[k],16)) {
	    fprintf(stderr, "fd ir not aligned\n");
	}
    }

    for (int k=0; k<7; ++k) {
        unit->fd_irs512[k] = (data+(128*3)+512*k);
        if (!is_aligned(unit->fd_irs512[k],16)) {
            fprintf(stderr, "fd ir not aligned\n");
        }
    }

    for (int k=0; k<numirs4096; ++k) {
        unit->fd_irs4096[k] = (data+(128*3)+(512*7)+4096*k);
        if (!is_aligned(unit->fd_irs4096[k],16)) {
            fprintf(stderr, "fd ir not aligned\n");
        }
    }

    unit->m_accummem = (float*) RTAlloc(unit->mWorld, sizeof(float)*(128+512+512+4096+4096)+16);
    //unit->m_accum = (float*) align_up(unit->m_accummem, 16);
    unit->m_accum = unit->m_accummem;
    unit->m_accum512 = unit->m_accummem+128;
    unit->m_accumother512 = unit->m_accummem+128+512;
    unit->m_read512 = unit->m_accumother512;
    unit->m_accidx512 = 1;

    unit->m_accum4096 = unit->m_accummem+128+1024;
    unit->m_accumother4096 = unit->m_accummem+128+1024+4096;
    unit->m_read4096 = unit->m_accumother4096;
    unit->m_accidx4096 = 1;

    if (!is_aligned(unit->m_accum,16)) {
        fprintf(stderr, "accum not aligned\n");
    }
    if (!is_aligned(unit->m_accum512,16)) {
        fprintf(stderr, "accum not aligned\n");
    }

    unit->m_inbuf = (float*) RTAlloc(unit->mWorld, sizeof(float)*8192);
    memset(unit->m_accum, 0, sizeof(float)*128);
    memset(unit->m_inbuf, 0, sizeof(float)*8192);
    memset(unit->m_accum512, 0, sizeof(float)*512);
    memset(unit->m_accum4096, 0, sizeof(float)*4096);
    memset(unit->m_accumother512, 0, sizeof(float)*512);
    memset(unit->m_accumother4096, 0, sizeof(float)*4096);

    fprintf(stderr, "make plans: %d\n", unit->inbufs_per_fdbuf);
    fprintf(stderr, "make plans 512: %d\n", unit->inbufs_per_fdbuf512);
    unit->fplan = fftwf_plan_r2r_1d(FFTSIZE1, unit->fdbufs[0], unit->fdbufs[0], FFTW_R2HC, FFTW_PATIENT);
    unit->bplan = fftwf_plan_r2r_1d(FFTSIZE1, unit->m_accum, unit->m_accum, FFTW_HC2R, FFTW_PATIENT);
    unit->fplan512 = fftwf_plan_r2r_1d(512, unit->fdbufs512[0], unit->fdbufs512[0], FFTW_R2HC, FFTW_PATIENT);
    unit->bplan512 = fftwf_plan_r2r_1d(512, unit->m_accum512, unit->m_accum512, FFTW_HC2R, FFTW_PATIENT);
    unit->fplan4096 = fftwf_plan_r2r_1d(4096, unit->fdbufs4096[0], unit->fdbufs4096[0], FFTW_R2HC, FFTW_PATIENT);
    unit->bplan4096 = fftwf_plan_r2r_1d(4096, unit->m_accum4096, unit->m_accum4096, FFTW_HC2R, FFTW_PATIENT);
}


// It is important that the given buffer is exactly the right size (this has to be ensured from the language side)
void LowLatConv_Ctor(LowLatConv *unit)
{
    LowLatConv_init(unit);
    //fprintf(stderr, "num fdbufs: %d\n", unit->num_fdbufs);
    unit->m_bypass = 0;
    SETCALC(LowLatConv_next);
    //fprintf(stderr, "bbb\n");
    ZOUT0(0) = 0.f;
}


#define MIN(_x,_y) (_x<_y?_x:_y)


void LowLatConv_next(LowLatConv *unit, int inNumSamples)
{
    float *in = IN(0);
    float *out = OUT(0);
    int bypass = ZIN0(2);

    if (bypass) {
        memset(out, 0, sizeof(float) * inNumSamples);
        unit->m_bypass = bypass;
        return;
    } else if (unit->m_bypass) {
        // reinit
        LowLatConv_init(unit);        
        unit->m_bypass = 0;        
    }

    float *accum = unit->m_accum;
    float *accum512;
    float *accum4096;
    int fd_idx = unit->fd_idx;
    int fd_idx512 = unit->fd_idx512;
    int fd_idx4096 = unit->fd_idx4096;
    int mask4096 = unit->m_4096mask;
    int numirs4096 = unit->m_4096irs;
    int inbufs_collected = unit->inbufs_collected;
    
    int inpos = (inbufs_collected*inNumSamples) & 8191;
    xmemcpy(unit->m_inbuf+inpos, in, sizeof(float)*inNumSamples);
    int numirs512 = unit->m_512irs;
    int ibm = (inbufs_collected % unit->inbufs_per_fdbuf512);

    inbufs_collected++;
    inpos+=inNumSamples;

    //fprintf(stderr, "--------------------------\naccum: %d\n accum512: %d\n accum4096: %d\n idx: %d\n idx512: %d\nidx4096: %d\n, mask: %d\n inpos: %d\nnumirs512L %d\n ibm: %d\n bypass: %d\ninbufs: %d\n diff1: %d diff2: %d\n",      unit->m_accum, unit->m_accum512, unit->m_accum4096, fd_idx, fd_idx512, fd_idx4096, mask4096, inpos, numirs512, ibm, bypass, inbufs_collected, unit->m_accum4096-unit->m_accum512, unit->m_accum-unit->m_accum512);

    // check if buffer is full (and process it)
    if((inbufs_collected % unit->inbufs_per_fdbuf) == 0) {        
	//inbufs_collected = 0;
	float *current_fdbuf = unit->fdbufs[fd_idx&3];
        //fprintf(stderr, "%d\n", inpos);
        // copy the data into the fft buffer
        memcpy(current_fdbuf+64, unit->m_inbuf+inpos-64, sizeof(float) * 64);
        memcpy(current_fdbuf, unit->m_inbuf+((inpos-128) & 8191), sizeof(float) * 64);
       	fftwf_execute_r2r(unit->fplan, current_fdbuf, current_fdbuf);
        accum = (float *) __builtin_assume_aligned(unit->m_accum, 16);
        // collect the frequency domain delay line buffers into an accumulator
        float *b = (float *) __builtin_assume_aligned(current_fdbuf, 16);
        float *c = (float *) __builtin_assume_aligned(unit->fd_irs[0], 16);
        // halfcomplex format multiply
        float a;
        accum[0] = b[0] * c[0];
        for (int n = 1; n < FFTSIZE1BY2; n++) {
            a = b[n];
            accum[n] = a * c[n] - b[FFTSIZE1 - n] * c[FFTSIZE1 - n];
            accum[FFTSIZE1 - n] = a * c[FFTSIZE1 - n] + b[FFTSIZE1 - n] * c[n];
        }
        accum[FFTSIZE1BY2] = b[FFTSIZE1BY2] * c[FFTSIZE1BY2];
        // add the previous frequency domain buffers in
        for (int i=1; i<3; ++i) {
            float *b = (float *) __builtin_assume_aligned(unit->fdbufs[(fd_idx-i)&3], 16);
            float *c = (float *) __builtin_assume_aligned(unit->fd_irs[i&3], 16);
            accum[0] += b[0] * c[0];
            for (int n = 1; n < FFTSIZE1BY2; n++) {
                a = b[n];
                accum[n] += a * c[n] - b[FFTSIZE1- n] * c[FFTSIZE1 - n];
                accum[FFTSIZE1 - n] += a * c[FFTSIZE1 - n] + b[FFTSIZE1 - n] * c[n];
            }
            accum[FFTSIZE1BY2] += b[FFTSIZE1BY2] * c[FFTSIZE1BY2];
        }
	unit->fd_idx = fd_idx = (fd_idx+1) & 3;
	fftwf_execute_r2r(unit->bplan, accum, accum);
    } 

    if(ibm == 0) {        
	float *current_fdbuf = unit->fdbufs512[fd_idx512&7];
        int space1 = 8192 - ((inpos-512)&8191);
        if (space1 >= 512) {
            memcpy(current_fdbuf, unit->m_inbuf+((inpos-512) & 8191), sizeof(float) * 512);
        } else {
            memcpy(current_fdbuf, unit->m_inbuf+((inpos-512) & 8191), sizeof(float) * space1);
            memcpy(current_fdbuf+space1, unit->m_inbuf+((inpos-512+space1) & 8191), sizeof(float) * (512-space1));
        }
       	fftwf_execute_r2r(unit->fplan512, current_fdbuf, current_fdbuf);
        accum512 = (float *) __builtin_assume_aligned(unit->m_accum512, 16);
        // collect the frequency domain delay line buffers into an accumulator
        float *b = (float *) __builtin_assume_aligned(current_fdbuf, 16);
        float *c = (float *) __builtin_assume_aligned(unit->fd_irs512[0], 16);
        // halfcomplex format multiply
        float a;
        accum512[0] += b[0] * c[0];
        for (int n = 1; n < 256; n++) {
            a = b[n];
            accum512[n] += a * c[n] - b[512 - n] * c[512 - n];
            accum512[512 - n] += a * c[512 - n] + b[512 - n] * c[n];
        }
        accum512[256] += b[256] * c[256];
        // add the previous frequency domain buffers in
        unit->m_accidx512 = numirs512-1;
	unit->fd_idx512 = fd_idx512 = (fd_idx512+1) & 7;
	fftwf_execute_r2r(unit->bplan512, accum512, accum512);
        float *tmp = unit->m_accumother512;
        unit->m_read512 = unit->m_accum512 + 256;
        unit->m_accumother512 = unit->m_accum512;
        unit->m_accum512 = tmp;
    } else {
        int i = unit->m_accidx512;
        int nir;
        // i-1 is the number of buffers that need to be finished in 
        // periods-1-ibm  (we get 3 2 1)
        // 3-3
        if ((i-(unit->inbufs_per_fdbuf512-ibm)) > 0) {
            nir = 2;
        } else {
            if (i>0) {
                nir = 1;
            } else {
                nir = 0;
            }
        }
        //fprintf(stderr, "ibm: %d i: %d nir: %d numirs: %d\n", ibm, i, nir, numirs512);
        accum512 = (float *) __builtin_assume_aligned(unit->m_accum512, 16);
        // collect the frequency domain delay line buffers into an accumulator
        float *b = (float *) __builtin_assume_aligned(unit->fdbufs512[(fd_idx512-i)&7], 16);
        float *c = (float *) __builtin_assume_aligned(unit->fd_irs512[i], 16);
        // halfcomplex format multiply
        float a;
        //fprintf(stderr, "accidx: %d\n", i);
        if ((ibm == 1) && (numirs512>1)) {
            accum512[0] = b[0] * c[0];
            for (int n = 1; n < 256; n++) {
                a = b[n];
                accum512[n] = a * c[n] - b[512 - n] * c[512 - n];
                accum512[512 - n] = a * c[512 - n] + b[512 - n] * c[n];
            }
            accum512[256] = b[256] * c[256];
            nir--;
            i--;
        } else {
            if (nir > 0) {
                accum512[0] += b[0] * c[0];
                for (int n = 1; n < 256; n++) {
                    a = b[n];
                    accum512[n] += a * c[n] - b[512 - n] * c[512 - n];
                    accum512[512 - n] += a * c[512 - n] + b[512 - n] * c[n];
                }
                accum512[256] += b[256] * c[256];
            }
            nir--;
            i--;
        }
        if (nir>0) {
            //accum512 = (float *) __builtin_assume_aligned(unit->m_accum512, 16);
            // collect the frequency domain delay line buffers into an accumulator
            b = (float *) __builtin_assume_aligned(unit->fdbufs512[(fd_idx512-i)&7], 16);
            c = (float *) __builtin_assume_aligned(unit->fd_irs512[i], 16);
            accum512[0] += b[0] * c[0];
            for (int n = 1; n < 256; n++) {
                a = b[n];
                accum512[n] += a * c[n] - b[512 - n] * c[512 - n];
                accum512[512 - n] += a * c[512 - n] + b[512 - n] * c[n];
            }
            accum512[256] += b[256] * c[256];
            i--;
        }
        unit->m_accidx512 = i;
    }

    // check if buffer is full (and process it)
    if((numirs4096 > 0) && ((inbufs_collected % unit->inbufs_per_fdbuf4096) == 0)) { 
        //fprintf(stderr, "fdidx: %d\n", fd_idx4096);
	float *current_fdbuf = unit->fdbufs4096[fd_idx4096&mask4096];
        //fprintf(stderr, "%d %d\n", inpos, inbufs_collected);
        // copy the data into the fft buffer
        int space1 = 8192 - ((inpos-4096)&8191);
        if (space1 >= 4096) {
            memcpy(current_fdbuf, unit->m_inbuf+((inpos-4096) & 8191), sizeof(float) * 4096);
        } else {
            memcpy(current_fdbuf, unit->m_inbuf+((inpos-4096) & 8191), sizeof(float) * space1);
            memcpy(current_fdbuf+space1, unit->m_inbuf+((inpos-4096+space1) & 8191), sizeof(float) * (4096-space1));
        }
        //fprintf(stderr, "%d %d %d\n", fd_idx4096, space1, unit->m_accidx4096);
       	fftwf_execute_r2r(unit->fplan4096, current_fdbuf, current_fdbuf);
        accum4096 = (float *) __builtin_assume_aligned(unit->m_accum4096, 16);
        // collect the frequency domain delay line buffers into an accumulator
        float * b = (float *) __builtin_assume_aligned(current_fdbuf, 16);
        float * c = (float *) __builtin_assume_aligned(unit->fd_irs4096[0], 16);
        // halfcomplex format multiply
        float a;
        if (numirs4096 == 1) {
            accum4096[0] = b[0] * c[0];
            for (int n = 1; n < 2048; n++) {
                a = b[n];
                accum4096[n] = a * c[n] - b[4096 - n] * c[4096 - n];
                accum4096[4096 - n] = a * c[4096 - n] + b[4096 - n] * c[n];
            }
            accum4096[2048] = b[2048] * c[2048];
        } else {
            accum4096[0] += b[0] * c[0];
            for (int n = 1; n < 2048; n++) {
                a = b[n];
                accum4096[n] += a * c[n] - b[4096 - n] * c[4096 - n];
                accum4096[4096 - n] += a * c[4096 - n] + b[4096 - n] * c[n];
            }
            accum4096[2048] += b[2048] * c[2048];
        }
        // add the previous frequency domain buffers in
        unit->m_accidx4096 = 1;
	unit->fd_idx4096 = fd_idx4096 = (fd_idx4096+1) & mask4096;
	fftwf_execute_r2r(unit->bplan4096, accum4096, accum4096);
        float *tmp = unit->m_accumother4096;
        unit->m_read4096 = unit->m_accum4096 + 2048;
        unit->m_accumother4096 = unit->m_accum4096;
        unit->m_accum4096 = tmp;
    } else {
        int i = unit->m_accidx4096;
        int ibm = (inbufs_collected % unit->inbufs_per_fdbuf4096);
        int nir = (numirs4096-i+1) / (unit->inbufs_per_fdbuf4096-ibm);
        if ((ibm == 1) && (numirs4096>1)) {
            //fprintf(stderr, "nir 1: %d %d %d %d\n", nir, ibm, fd_idx4096, mask4096);
            accum4096 = (float *) __builtin_assume_aligned(unit->m_accum4096, 16);
            // collect the frequency domain delay line buffers into an accumulator
            float *b = (float *) __builtin_assume_aligned(unit->fdbufs4096[(fd_idx4096-i)&mask4096], 16);
            float *c = (float *) __builtin_assume_aligned(unit->fd_irs4096[i], 16);
            // halfcomplex format multiply
            float a;
            //fprintf(stderr, "accidx: %d\n", i);
            accum4096[0] = b[0] * c[0];
            for (int n = 1; n < 2048; n++) {
                a = b[n];
                accum4096[n] = a * c[n] - b[4096 - n] * c[4096 - n];
                accum4096[4096 - n] = a * c[4096 - n] + b[4096 - n] * c[n];
            }
            accum4096[2048] = b[2048] * c[2048];
            i++;
            nir--;
        }
        // here
        while ((nir > 0) && (i < numirs4096)) {
            accum4096 = (float *) __builtin_assume_aligned(unit->m_accum4096, 16);
            // collect the frequency domain delay line buffers into an accumulator
            float *b = (float *) __builtin_assume_aligned(unit->fdbufs4096[(fd_idx4096-i)&mask4096], 16);
            float *c = (float *) __builtin_assume_aligned(unit->fd_irs4096[i], 16);
            // halfcomplex format multiply
            float a;
            accum4096[0] += b[0] * c[0];
            for (int n = 1; n < 2048; n++) {
                a = b[n];
                accum4096[n] += a * c[n] - b[4096 - n] * c[4096 - n];
                accum4096[4096 - n] += a * c[4096 - n] + b[4096 - n] * c[n];
            }
            accum4096[2048] += b[2048] * c[2048];
            //fprintf(stderr, "nir 2: %d %d %d %d\n", nir, ibm, fd_idx4096, mask4096);
            i++;
            nir--;
        }
        unit->m_accidx4096 = i;
    }

    accum += 64;
    float *read512 = unit->m_read512;
    float *read4096 = unit->m_read4096;
    // copy inNumSamples into the time domain output of the unit
    if (numirs4096 > 0) {
        for (int k=0; k<inNumSamples; ++k) {
            *out++ = *accum++ + *read512++ + *read4096++;
        }
    } else {
        for (int k=0; k<inNumSamples; ++k) {
            *out++ = *accum++ + *read512++;
        }
    }
    
    unit->m_read512 = read512;
    unit->m_read4096 = read4096;
    unit->inbufs_collected = inbufs_collected % 1024;
}



//channels not used- should just be mono, num frames= num samples
//buffer preparation
void BPrepareLowLatConv(World *world, struct SndBuf *buf, struct sc_msg_iter *msg)
{
    float *data1 = buf->data;
    uint32 frombufnum = msg->geti();
    //output size must be frombuf->frames*2
    if (frombufnum >= world->mNumSndBufs) frombufnum = 0;
    SndBuf* frombuf = world->mSndBufs + frombufnum;
    int frames2 = frombuf->frames;
    float *data2 = frombuf->data;
    int numirs512 = 7;
    int extrair4096 = 0;
    int extrair512 = 0;
    memset(data1, 0, sizeof(float)*buf->frames);

    int numirs4096 = 0;
    if (frames2 > 1984) {
        numirs4096 = (frames2-1984) / 2048;
        if (frames2 > (1984+numirs4096*2048)) {
            extrair4096 = frames2 - (1984+numirs4096*2048);
        }
        //num4096 = nextpow2(numirs4096);  
    } else {
        numirs512 = (frames2-64*3) / 256;
        if (frames2 > (64*3+numirs512*256)) {
            extrair512 = frames2-(64*3+numirs512*256);
        }
    }

    // make sure we can hold the additional size for the spectrum as required by fftw
    float * inputbuf = (float*)RTAlloc(world, 4096 * sizeof(float));
    fftwf_plan theplan = fftwf_plan_r2r_1d(FFTSIZE1, inputbuf, inputbuf, FFTW_R2HC, FFTW_ESTIMATE);

    for (int i=0; i<3; ++i) {
	//int indexnow= nover2*i;
        memset(inputbuf, 0, sizeof(float) * FFTSIZE1);
        memcpy(inputbuf, data2+(i*FFTSIZE1BY2), sizeof(float) * FFTSIZE1BY2);
       	fftwf_execute(theplan);
        for (int k=0; k<FFTSIZE1; k++) {
            inputbuf[k] /= ((float) FFTSIZE1);
        }
	memcpy(data1+(i*FFTSIZE1), inputbuf, sizeof(float) * FFTSIZE1);
    }

    fftwf_destroy_plan(theplan);
    theplan = fftwf_plan_r2r_1d(512, inputbuf, inputbuf, FFTW_R2HC, FFTW_ESTIMATE);

    for (int i=0; i<numirs512; ++i) {
	//int indexnow= nover2*i;
        memset(inputbuf, 0, sizeof(float) * 512);
        memcpy(inputbuf, data2+192+(i*256), sizeof(float) * 256);
       	fftwf_execute(theplan);
        for (int k=0; k<512; k++) {
            inputbuf[k] /= ((float) 512);
        }
	memcpy(data1+384+(i*512), inputbuf, sizeof(float) * 512);
    }

    if (extrair512 > 0) {
        memset(inputbuf, 0, sizeof(float) * 512);
        memcpy(inputbuf, data2+192+(numirs512*256), sizeof(float) * extrair512);
       	fftwf_execute(theplan);
        for (int k=0; k<512; k++) {
            inputbuf[k] /= ((float) 512);
        }
	memcpy(data1+384+(numirs512*512), inputbuf, sizeof(float) * 512);        
        numirs512++;
    }

    fftwf_destroy_plan(theplan);
    theplan = fftwf_plan_r2r_1d(4096, inputbuf, inputbuf, FFTW_R2HC, FFTW_ESTIMATE);

    for (int i=0; i<numirs4096; ++i) {
        memset(inputbuf, 0, sizeof(float) * 4096);
        memcpy(inputbuf, data2+192+(7*256)+(i*2048), sizeof(float) * 2048);
       	fftwf_execute(theplan);
        for (int k=0; k<4096; k++) {
            inputbuf[k] /= ((float) 4096);
        }
	memcpy(data1+384+(7*512)+(i*4096), inputbuf, sizeof(float) * 4096);
    }

    if (extrair4096 > 0) {
        memset(inputbuf, 0, sizeof(float) * 4096);
        memcpy(inputbuf, data2+192+(7*256)+(numirs4096*2048), sizeof(float) * extrair4096);
       	fftwf_execute(theplan);
        for (int k=0; k<4096; k++) {
            inputbuf[k] /= ((float) 4096);
        }
	memcpy(data1+384+(7*512)+(numirs4096*4096), inputbuf, sizeof(float) * 4096);        
        numirs4096++;
    }

    // fprintf(stderr, "numirs: %d %d %d %d %d %d\n", 384+(numirs512*512)+(numirs4096*4096), 
    //        +192+(numirs512*256)+(numirs4096*2048),
    //        buf->frames, frombuf->frames, numirs4096, numirs512);

    fftwf_destroy_plan(theplan);
    RTFree(world, inputbuf);
}


void load(InterfaceTable *inTable)
{
    ft = inTable;
    DefineDtorUnit(LowLatConv);
    DefineBufGen("BPrepareLowLatConv", BPrepareLowLatConv);
}
