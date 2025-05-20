#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h> // strtok, strstr
#define SAMPLES_PER_FILE 6000
// typedef __fp16 float16;
#define MAX_N 256
#define PI 3.14159265358979323846
#define CUSTOM_FLOAT16_MAX_EXP 127  // Max exponent for custom float16
#define CUSTOM_FLOAT16_MIN_EXP -126 // Min exponent for custom float16
#define CUSTOM_FLOAT_SIGN(x) ((x.value >> 15) & 0x1)
#define CUSTOM_FLOAT_EXP(x) ((x.value >> 7) & 0xFF)
#define CUSTOM_FLOAT_MANT(x) (x.value & 0x7F)
#define CUSTOM_FLOAT16_EXP_BIAS 127
#define CUSTOM_FLOAT_BIAS 127
#include <stdbool.h>

typedef struct
{
    uint16_t value;
} custom_float16_t;

typedef struct
{
    custom_float16_t re;
    custom_float16_t im;
} complex32_t;
typedef struct
{
    float re;
    float im;
} complex32;
custom_float16_t float_to_custom_float16(float f)
{
    custom_float16_t out = {0};

    if (f == 0.0f)
    {
        return out;
    }

    int sign = (f < 0);
    float abs_f = fabsf(f);

    int exp;
    float frac = frexpf(abs_f, &exp); // abs_f = frac * 2^exp, 0.5 ≤ frac < 1

    // Shift to get 8-bit mantissa (1 implicit + 7 bits stored)
    // Round to nearest instead of truncating
    uint16_t mantissa = (uint16_t)(frac * (1 << 8) + 0.5f); // Rounding added here

    // Handle rounding overflow
    if (mantissa >= (1 << 8))
    {
        mantissa >>= 1;
        exp += 1;
    }

    int biased_exp = exp + (CUSTOM_FLOAT16_EXP_BIAS - 1); // same as IEEE style

    // Clamp
    if (biased_exp <= 0)
    {
        out.value = 0;
        return out;
    }
    if (biased_exp >= 255)
    {
        // Clamp to max representable
        biased_exp = 254;
        mantissa = (1 << 8) - 1;
    }

    // Pack
    out.value = (sign << 15) | ((biased_exp & 0xFF) << 7) | ((mantissa - (1 << 7)) & 0x7F);
    return out;
}

// Convert custom float16 to float
float custom_float16_to_float(custom_float16_t val)
{
    if (val.value == 0)
        return 0.0f;

    int sign = (val.value >> 15) & 0x1;
    int exponent = (val.value >> 7) & 0xFF;
    int mantissa = val.value & 0x7F;

    float f_mant = 1.0f + mantissa / 128.0f; // restore implicit 1
    int real_exp = exponent - CUSTOM_FLOAT16_EXP_BIAS;

    float result = ldexpf(f_mant, real_exp);
    return sign ? -result : result;
}

bool get_sign(custom_float16_t a)
{
    return (a.value >> 15) & 1;
}

uint8_t custom_float_get_mantissa(custom_float16_t f)
{
    return f.value & 0x7F;
}

uint8_t custom_float_get_exp(custom_float16_t f)
{
    return (f.value >> 7) & 0xFF;
}
bool is_nan(custom_float16_t f)
{
    return custom_float_get_exp(f) == 0xFF && custom_float_get_mantissa(f) != 0;
}

bool is_inf(custom_float16_t f)
{
    return custom_float_get_exp(f) == 0xFF && custom_float_get_mantissa(f) == 0;
}
bool is_custom_float_zero(custom_float16_t f)
{
    // If exponent and mantissa are zero (sign doesn't matter)
    return (f.value & 0x7FFF) == 0;
}

// // Main Addition Function
// float add_custom_float(custom_float16_t a, custom_float16_t b)
// {
//     // Zero handling
//     if (is_custom_float_zero(a) && is_custom_float_zero(b))
//         return 0.0f;

//     if (is_custom_float_zero(a))
//         return custom_float16_to_float(b); // we'll define this
//     if (is_custom_float_zero(b))
//         return custom_float16_to_float(a);

//     uint8_t exp1 = custom_float_get_exp(a) - CUSTOM_FLOAT_BIAS;
//     uint8_t exp2 = custom_float_get_exp(b) - CUSTOM_FLOAT_BIAS;
//     uint8_t mantissa1 = custom_float_get_mantissa(a);
//     uint8_t mantissa2 = custom_float_get_mantissa(b);
//     bool sign1 = get_sign(a);
//     bool sign2 = get_sign(b);
//     // Convert to actual float values
//     float m1 = 1.0f + ((float)mantissa1 / 128.0f); // 2^7 = 128
//     float m2 = 1.0f + ((float)mantissa2 / 128.0f);

//     float num1 = (sign1 ? -1.0f : 1.0f) * m1 * powf(2.0f, exp1);
//     float num2 = (sign2 ? -1.0f : 1.0f) * m2 * powf(2.0f, exp2);

//     return num1 + num2;
// }
// // ... [Keep all remaining code unchanged]
// float sub_custom_float(custom_float16_t a, custom_float16_t b)
// {
//     custom_float16_t neg_b = {(uint16_t)(b.value ^ (1 << 15))};
//     // custom_float16_t neg_b = {b.value ^ 0x8000};
//     return add_custom_float(a, neg_b);
// }

// custom_float16_t add_custom_float_pure(custom_float16_t a, custom_float16_t b)
// {
//     // Handle zero cases
//     if (is_custom_float_zero(a))
//         return b;
//     if (is_custom_float_zero(b))
//         return a;

//     // Extract sign, exponent, and mantissa (with implicit 1)
//     int sign_a = get_sign(a);
//     int sign_b = get_sign(b);
//     int exp_a = custom_float_get_exp(a);
//     int exp_b = custom_float_get_exp(b);
//     int mant_a = (1 << 7) | custom_float_get_mantissa(a); // implicit 1
//     int mant_b = (1 << 7) | custom_float_get_mantissa(b); // implicit 1

//     // Align exponents
//     int exp_diff = exp_a - exp_b;
//     if (exp_diff > 0)
//     {
//         mant_b >>= exp_diff;
//         exp_b = exp_a;
//     }
//     else if (exp_diff < 0)
//     {
//         mant_a >>= -exp_diff;
//         exp_a = exp_b;
//     }

//     // Perform operation
//     int result_sign;
//     int result_mant;
//     if (sign_a == sign_b)
//     {
//         result_mant = mant_a + mant_b;
//         result_sign = sign_a;
//     }
//     else
//     {
//         if (mant_a > mant_b)
//         {
//             result_mant = mant_a - mant_b;
//             result_sign = sign_a;
//         }
//         else
//         {
//             result_mant = mant_b - mant_a;
//             result_sign = sign_b;
//         }
//     }

//     // Normalize
//     int result_exp = exp_a;
//     while (result_mant && result_mant < (1 << 7))
//     {
//         result_mant <<= 1;
//         result_exp--;
//     }

//     // Rounding (cut to 7 bits)
//     int rounded_mant = result_mant & 0x7F;

//     // Handle underflow or overflow
//     if (result_exp <= 0)
//         return (custom_float16_t){0};
//     if (result_exp >= 0xFF)
//         return (custom_float16_t){(result_sign << 15) | (0xFE << 7) | 0x7F};

//     // Pack result
//     custom_float16_t result;
//     result.value = (result_sign << 15) | ((result_exp & 0xFF) << 7) | rounded_mant;
//     return result;
// }

// custom_float16_t sub_custom_float_pure(custom_float16_t a, custom_float16_t b)
// {
//     custom_float16_t neg_b = {b.value ^ (1 << 15)}; // Flip sign bit
//     return add_custom_float_pure(a, neg_b);
// }
typedef struct
{
    int sign;
    int exponent;
    int mantissa; // includes implicit 1
} unpacked_custom_float_t;

unpacked_custom_float_t unpack_custom_float16(custom_float16_t f)
{
    unpacked_custom_float_t result;
    result.sign = (f.value >> 15) & 0x1;
    result.exponent = (f.value >> 7) & 0xFF;
    int raw_mantissa = f.value & 0x7F;

    if (result.exponent == 0)
    {
        result.mantissa = raw_mantissa; // subnormal
    }
    else
    {
        result.mantissa = (1 << 7) | raw_mantissa; // add implicit 1
    }

    return result;
}

custom_float16_t pack_custom_float16(int sign, int exponent, int mantissa)
{
    custom_float16_t result;

    // Normalize mantissa and adjust exponent
    while (mantissa && mantissa < (1 << 7))
    {
        mantissa <<= 1;
        exponent--;
    }

    // Round mantissa to 7 bits
    int final_mantissa = mantissa & 0x7F;

    // Clamp exponent to valid range
    if (exponent <= 0)
    {
        return (custom_float16_t){0};
    }
    else if (exponent >= 0xFF)
    {
        return (custom_float16_t){(uint16_t)((sign << 15) | (0xFE << 7) | 0x7F)}; // max
    }

    result.value = (sign << 15) | ((exponent & 0xFF) << 7) | final_mantissa;
    return result;
}
custom_float16_t add_custom_float(custom_float16_t a, custom_float16_t b)
{
    if (is_custom_float_zero(a))
        return b;
    if (is_custom_float_zero(b))
        return a;

    unpacked_custom_float_t ua = unpack_custom_float16(a);
    unpacked_custom_float_t ub = unpack_custom_float16(b);

    // Align exponents
    while (ua.exponent > ub.exponent)
    {
        ub.mantissa >>= 1;
        ub.exponent++;
    }
    while (ub.exponent > ua.exponent)
    {
        ua.mantissa >>= 1;
        ua.exponent++;
    }

    int result_mantissa;
    int result_sign;
    if (ua.sign == ub.sign)
    {
        result_mantissa = ua.mantissa + ub.mantissa;
        result_sign = ua.sign;
    }
    else
    {
        if (ua.mantissa >= ub.mantissa)
        {
            result_mantissa = ua.mantissa - ub.mantissa;
            result_sign = ua.sign;
        }
        else
        {
            result_mantissa = ub.mantissa - ua.mantissa;
            result_sign = ub.sign;
        }
    }

    return pack_custom_float16(result_sign, ua.exponent, result_mantissa);
}
custom_float16_t sub_custom_float(custom_float16_t a, custom_float16_t b)
{
    b.value ^= (1 << 15); // Flip sign bit
    return add_custom_float(a, b);
}

// Convert float to custom float16
// Convert float to custom float16
custom_float16_t mul_custom_float(custom_float16_t a, custom_float16_t b)
{
    float fa = custom_float16_to_float(a);
    float fb = custom_float16_to_float(b);
    float product = fa * fb;
    return float_to_custom_float16(product);
}

custom_float16_t div_custom_float(custom_float16_t a, custom_float16_t b)
{
    float fa = custom_float16_to_float(a);
    float fb = custom_float16_to_float(b);
    float quotient = fa / fb;
    return float_to_custom_float16(quotient);
}

float get_re(complex32_t c)
{
    return custom_float16_to_float(c.re);
}

float get_im(complex32_t c)
{
    return custom_float16_to_float(c.im);
}

void set_re(complex32_t *c, float val)
{
    c->re = float_to_custom_float16(val);
}

void set_im(complex32_t *c, float val)
{
    c->im = float_to_custom_float16(val);
}

/* 2-point FFT “butterfly” */
void fft2(complex32_t x[2])
{
    custom_float16_t a_re = x[0].re;
    custom_float16_t a_im = x[0].im;
    custom_float16_t b_re = x[1].re;
    custom_float16_t b_im = x[1].im;
    /* X0 = a + b */
    x[0].re = add_custom_float(a_re, b_re);
    x[0].im = add_custom_float(a_im, b_im);
    /* X1 = a - b */
    x[1].re = sub_custom_float(a_re, b_re);
    x[1].im = sub_custom_float(a_im, b_im);
}

/* Recursive radix-2 FFT using only fft2 */
void fft(complex32_t *x, int N)
{
    if (N == 2)
    {
        fft2(x);
        return;
    }
    int M = N / 2;

    /* split into even/odd (on the stack via VLA) */
    complex32_t even[M], odd[M];
    for (int i = 0; i < M; i++)
    {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    /* recurse till base case is hit*/
    fft(even, M);
    fft(odd, M);

    /* combine butterflies via 2-point FFT */
    for (int k = 0; k < M; k++)
    {
        float angle = -2.0 * PI * k / N;
        float w_re = cos(angle);
        float w_im = sin(angle);
        /* twiddle * odd[k] */
        float t_re = w_re * custom_float16_to_float(odd[k].re) - w_im * custom_float16_to_float(odd[k].im);
        float t_im = w_re * custom_float16_to_float(odd[k].im) + w_im * custom_float16_to_float(odd[k].re);
        /* form pair */
        complex32_t pair[2];
        pair[0].re = even[k].re;
        pair[0].im = even[k].im;
        /* cast twiddled values back to 16-bit */
        pair[1].re = float_to_custom_float16(t_re);
        pair[1].im = float_to_custom_float16(t_im);
        /* one fft2 on this pair */
        fft2(pair);
        /* store results */
        x[k].re = pair[0].re;
        x[k].im = pair[0].im;
        x[k + M].re = pair[1].re;
        x[k + M].im = pair[1].im;
    }
}
/* 2-point FFT “butterfly” */
void fft2_32(complex32 x[2])
{
    float a_re = x[0].re;
    float a_im = x[0].im;
    float b_re = x[1].re;
    float b_im = x[1].im;
    /* X0 = a + b */
    x[0].re = a_re + b_re;
    x[0].im = a_im + b_im;
    /* X1 = a - b */
    x[1].re = a_re - b_re;
    x[1].im = a_im - b_im;
}

/* Recursive radix-2 FFT using only fft2 */
void fft32(complex32 *x, int N)
{
    if (N == 2)
    {
        fft2_32(x);
        return;
    }
    int M = N / 2;

    /* split into even/odd (on the stack via VLA) */
    complex32 even[M], odd[M];
    for (int i = 0; i < M; i++)
    {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    /* recurse till base case is hit*/
    fft32(even, M);
    fft32(odd, M);

    /* combine butterflies via 2-point FFT */
    for (int k = 0; k < M; k++)
    {
        float angle = -2.0 * PI * k / N;
        float w_re = cos(angle);
        float w_im = sin(angle);
        /* twiddle * odd[k] */
        float t_re = w_re * odd[k].re - w_im * odd[k].im;
        float t_im = w_re * odd[k].im + w_im * odd[k].re;
        /* form pair */
        complex32 pair[2];
        pair[0].re = even[k].re;
        pair[0].im = even[k].im;
        /* cast twiddled values back to 16-bit */
        pair[1].re = t_re;
        pair[1].im = t_im;
        /* one fft2 on this pair */
        fft2_32(pair);
        /* store results */
        x[k].re = pair[0].re;
        x[k].im = pair[0].im;
        x[k + M].re = pair[1].re;
        x[k + M].im = pair[1].im;
    }
}
void apply_hamming_window(float *data)
{
    for (int i = 0; i < MAX_N; i++)
    {
        data[i] *= 0.54 - 0.46 * cos(2 * PI * i / (MAX_N - 1)); // Hamming window
    }
}

int main(void)
{
    FILE *fp = fopen("min_max_results.txt", "r");
    if (!fp)
    {
        perror("fopen");
        return EXIT_FAILURE;
    }

    char header[512], data_line[200000], filename[MAX_N];
    float samples[SAMPLES_PER_FILE];
    float window[MAX_N];
    complex32_t buf[MAX_N];
    complex32 float_arr[MAX_N];
    float abs_diff = 0;
    int total_bins = 0;
    float fmin = INT_MAX;
    float fmax = INT_MIN;
    float cmin = INT16_MAX;
    float cmax = INT16_MIN;

    while (fgets(header, sizeof(header), fp))
    {
        // skip any non-.wav lines
        if (!strstr(header, ".wav"))
            continue;

        // extract filename
        if (sscanf(header, "%255[^,], Min: %*f, Max: %*f", filename) != 1)
            continue;
        printf(">>> Processing file: %s\n", filename);

        // read exactly one line of up to 6000 floats
        if (!fgets(data_line, sizeof(data_line), fp))
            break;
        int idx = 0;
        for (char *tok = strtok(data_line, " \t\n"); tok && idx < SAMPLES_PER_FILE; tok = strtok(NULL, " \t\n"))
        {
            samples[idx++] = strtof(tok, NULL);
        }
        // pad missing
        for (; idx < SAMPLES_PER_FILE; idx++)
            samples[idx] = 0.0f;

        // 24 windows
        int windows = (SAMPLES_PER_FILE + MAX_N - 1) / MAX_N;
        for (int w = 0; w < windows; w++)
        {
            int start = w * MAX_N;
            int cnt = (SAMPLES_PER_FILE - start >= MAX_N ? MAX_N : SAMPLES_PER_FILE - start);
            // copy+pad
            for (int i = 0; i < cnt; i++)
                window[i] = samples[start + i];
            for (int i = cnt; i < MAX_N; i++)
                window[i] = 0.0f;

            apply_hamming_window(window);

            // load into custom‐float16 complex buffer
            for (int i = 0; i < MAX_N; i++)
            {
                buf[i].re = float_to_custom_float16(window[i]);
                buf[i].im = float_to_custom_float16(0.0f);
                float_arr[i].re = window[i];
                float_arr[i].im = 0.0f;
            }

            // do custom‐float16 FFT
            fft(buf, MAX_N);
            fft32(float_arr, MAX_N);

            printf("  Window %2d/%2d FFT results:\n", w + 1, windows);
            for (int i = 0; i < MAX_N; i++)
            {
                float r = custom_float16_to_float(buf[i].re);
                float m = custom_float16_to_float(buf[i].im);
                float f1 = float_arr[i].re;
                if (r > cmax)
                {
                    cmax = r;
                }
                if (r < cmin)
                {
                    cmin = r;
                }
                if (m > cmax)
                {
                    cmax = m;
                }
                if (m < cmin)
                {
                    cmin = m;
                }
                if (f1 > fmax)
                {
                    fmax = f1;
                }
                if (f1 < fmin)
                {
                    fmin = f1;
                }
                float f2 = float_arr[i].im;
                if (f2 > fmax)
                {
                    fmax = f2;
                }
                if (f2 < fmin)
                {
                    fmin = f2;
                }
                float mag_float = pow(pow(f1, 2.0) + pow(f2, 2.0), 1.0 / 2);
                float mag_custom = pow(pow(r, 2.0) + pow(m, 2.0), 1.0 / 2);
                const float EPS = 1e-6f;
                float diff = fabs(mag_custom - mag_float) / (EPS + mag_float);
                // diff /= mag_float;
                // diff *= 100.0f;
                printf("16 bit :\n");
                printf("[%3d] %f + j%f\n", i + 1, r, m);
                printf("32 bit :\n");
                printf("[%3d] %f + j%f\n", i + 1, f1, f2);
                printf("Diff : %f\n", diff);
                total_bins++;
                abs_diff += diff;
            }
        }
    }

    printf("\n=== GLOBAL SUMMARY ===\nTotal FFT bins: %d\n", total_bins);
    abs_diff /= total_bins;
    printf("Abs diff %f\n", abs_diff);
    printf("%f %f\n", cmin, cmax);
    printf("%f %f\n", fmin, fmax);
    fclose(fp);
    float tot = 0;
    // int c = 0;
    float a = 38.3532;
    float b = 2.9316;
    custom_float16_t A = float_to_custom_float16(a);
    custom_float16_t B = float_to_custom_float16(b);
    custom_float16_t C = add_custom_float(A, B);
    float c1 = custom_float16_to_float(C);
    float c0 = a + b;
    printf("add(%.3f,%.3f): float=%.3f  custom=%.3f  err=%.3f%%\n",
           a, b, c0, c1, 100 * fabsf(c1 - c0) / (fabsf(c0) + 1e-6f));

    return 0;
}
