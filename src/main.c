#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define M_PI 3.14159265358979323846

#define N_FFT 128
#define N_SYM 62
#define N_HALF 31
#define N_SAMPLES 19200
#define RUNS 100
#define FILE_NAME "signal.pcm"

void gen_pss(double pss[3][N_SYM]) {
    int u[] = {25, 29, 34};
    for (int nid2 = 0; nid2 < 3; nid2++) {
        for (int n = 0; n < N_SYM; n++) {
            pss[nid2][n] = -M_PI * u[nid2] * (n + 31) * (n + 32) / 63;
        }
    }
}

int detect_pss(const double* signal, double pss[3][N_SYM], double* max_corr) {
    double best_corr = 0;
    int best_nid2 = 0;
    
    for (int nid2 = 0; nid2 < 3; nid2++) {
        double corr_real = 0, corr_imag = 0;
        for (int n = 0; n < N_SYM; n++) {
            double phase = pss[nid2][n];
            corr_real += signal[2*n] * cos(phase) - signal[2*n+1] * sin(phase);
            corr_imag += signal[2*n] * sin(phase) + signal[2*n+1] * cos(phase);
        }
        double corr = sqrt(corr_real * corr_real + corr_imag * corr_imag);
        if (corr > best_corr) {
            best_corr = corr;
            best_nid2 = nid2;
        }
    }
    *max_corr = best_corr;
    return best_nid2;
}

void fft(double* re, double* im, int n) {
    if (n <= 1) return;
    
    int half = n / 2;
    double re_even[64], im_even[64], re_odd[64], im_odd[64];
    
    for (int i = 0; i < half; i++) {
        re_even[i] = re[2*i];
        im_even[i] = im[2*i];
        re_odd[i] = re[2*i+1];
        im_odd[i] = im[2*i+1];
    }
    
    fft(re_even, im_even, half);
    fft(re_odd, im_odd, half);
    
    for (int k = 0; k < half; k++) {
        double angle = -2 * M_PI * k / n;
        double wr = cos(angle), wi = sin(angle);
        double tr = wr * re_odd[k] - wi * im_odd[k];
        double ti = wr * im_odd[k] + wi * re_odd[k];
        re[k] = re_even[k] + tr;
        im[k] = im_even[k] + ti;
        re[k + half] = re_even[k] - tr;
        im[k + half] = im_even[k] - ti;
    }
}

void gen_m_sequence(int* s, int len) {
    int x[31] = {0, 0, 0, 0, 1};
    for (int j = 0; j < 26; j++) {
        x[j + 5] = (x[j + 2] + x[j]) % 2;
    }
    for (int i = 0; i < len; i++) {
        s[i] = 1 - 2 * x[i];
    }
}

void cyclic_shift(int* s, int m, int* out, int len) {
    for (int n = 0; n < len; n++) {
        out[n] = s[(n + m) % 31];
    }
}

void gen_scrambling_c(int nid2, int* c0, int* c1) {
    int x[31] = {0, 0, 0, 0, 1};
    for (int j = 0; j < 26; j++) {
        x[j + 5] = (x[j + 3] + x[j]) % 2;
    }
    for (int i = 0; i < 31; i++) {
        c0[i] = 1 - 2 * x[(i + nid2) % 31];
        c1[i] = 1 - 2 * x[(i + nid2 + 3) % 31];
    }
}

void gen_scrambling_z(int m0, int* z1) {
    int x[31] = {0, 0, 0, 0, 1};
    for (int j = 0; j < 26; j++) {
        x[j + 5] = (x[j + 4] + x[j + 2] + x[j + 1] + x[j]) % 2;
    }
    for (int i = 0; i < 31; i++) {
        z1[i] = 1 - 2 * x[(i + (m0 % 8)) % 31];
    }
}

void gen_sss(int nid1, int nid2, int* d) {
    int s[31], s0[31], s1[31];
    int c0[31], c1[31], z1[31];
    
    gen_m_sequence(s, 31);
    
    int m0 = 15 * (nid1 / 112) + 5 * (nid1 % 8);
    int m1 = (nid1 % 112) + 5 * (nid1 / 112);
    
    cyclic_shift(s, m0, s0, 31);
    cyclic_shift(s, m1, s1, 31);
    
    gen_scrambling_c(nid2, c0, c1);
    gen_scrambling_z(m0, z1);
    
    for (int n = 0; n < 31; n++) {
        d[2 * n]     = s0[n] * c0[n];
        d[2 * n + 1] = s1[n] * c1[n] * z1[n];
    }
}

int detect_pci(const double* signal, int nid2) {
    double re[N_FFT], im[N_FFT];
    for (int i = 0; i < N_FFT; i++) {
        re[i] = signal[2*i];
        im[i] = signal[2*i+1];
    }
    fft(re, im, N_FFT);
    
    int start = N_FFT / 2 - N_SYM / 2;
    double max_corr = 0;
    int best_nid1 = 0;
    
    int d[62];
    
    for (int nid1 = 0; nid1 < 168; nid1++) {
        gen_sss(nid1, nid2, d);
        
        double corr_real = 0, corr_imag = 0;
        for (int k = 0; k < N_SYM; k++) {
            corr_real += re[start + k] * d[k];
            corr_imag += im[start + k] * d[k];
        }
        
        double corr = sqrt(corr_real * corr_real + corr_imag * corr_imag);
        
        if (corr > max_corr) {
            max_corr = corr;
            best_nid1 = nid1;
        }
    }
    
    return 3 * best_nid1 + nid2;
}

double* read_iq_file(const char* filename, int* samples_read) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        printf("Warning: Cannot open '%s'\n", filename);
        return NULL;
    }
    
    short* buf = (short*)malloc(N_SAMPLES * 2 * sizeof(short));
    if (!buf) {
        fclose(f);
        return NULL;
    }
    
    size_t read = fread(buf, sizeof(short), N_SAMPLES * 2, f);
    fclose(f);
    *samples_read = (int)read / 2;
    
    if (*samples_read < N_SAMPLES) {
        printf("Warning: Only %d complex samples (expected %d)\n", *samples_read, N_SAMPLES);
    }
    
    double* signal = (double*)malloc(N_SAMPLES * 2 * sizeof(double));
    if (!signal) {
        free(buf);
        return NULL;
    }
    
    for (int i = 0; i < N_SAMPLES * 2 && i < read; i++) {
        signal[i] = buf[i] / 32768.0;
    }
    
    free(buf);
    return signal;
}

int main() {
    printf("\nLTE PCI Detector\n");
    printf("FFT: %d, Symbols: %d, Samples: %d\n\n", N_FFT, N_SYM, N_SAMPLES);
    
    int samples_read = 0;
    double* signal = read_iq_file(FILE_NAME, &samples_read);
    
    double pss[3][N_SYM];
    gen_pss(pss);
    
    printf("Running %d iterations...\n", RUNS);
    clock_t start_time = clock();
    
    int results[RUNS];
    double times[RUNS];
    
    for (int run = 0; run < RUNS; run++) {
        clock_t run_start = clock();
        
        double max_corr;
        int nid2 = detect_pss(signal, pss, &max_corr);
        results[run] = detect_pci(signal + N_FFT * 2, nid2);
        
        clock_t run_end = clock();
        times[run] = (double)(run_end - run_start) / CLOCKS_PER_SEC * 1000;
    }
    
    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    
    printf("\nResults\n");
    int pci = results[0];
    printf("Physical Cell ID (PCI): %d\n", pci);
    printf("  NID(2): %d (from PSS)\n", pci % 3);
    printf("  NID(1): %d (from SSS)\n", pci / 3);
    
    double min_time = times[0], max_time = times[0];
    double sum_time = 0;
    for (int i = 0; i < RUNS; i++) {
        sum_time += times[i];
        if (times[i] < min_time) min_time = times[i];
        if (times[i] > max_time) max_time = times[i];
    }
    double avg_time = sum_time / RUNS;
    
    printf("Performance (%d runs):\n", RUNS);
    printf("  Average: %.3f ms\n", avg_time);
    printf("  Min:     %.3f ms\n", min_time);
    printf("  Max:     %.3f ms\n", max_time);
    printf("  Total:   %.3f ms\n\n", total_time);
    
    free(signal);
    return 0;
}
