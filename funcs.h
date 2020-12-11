#ifndef FUNCS_H
#define FUNCS_H

#include <stdio.h>
#include <gmp.h>
#include <string.h>

//#define par_b "105598554324686606448243902967547806710936194579460793109384049876194488547338"
#define par_b "9774"
#define par_mod "115792089237316195423570985008687907853269984665640564039457584007913111864739"
#define par_q "115792089237316195423570985008687907852907080286716537199505774922796921406320"

struct point{
    mpz_t x, y, z;
};

struct curve_params{
    mpz_t a, d, mod;
};

int mpz_sqrtm(mpz_t x, const mpz_t y , const mpz_t z);

void get_d(mpz_t x, mpz_t mod);

void mpz_negm(mpz_t rop, mpz_t x, mpz_t mod);

void mpz_addm(mpz_t rop,mpz_t x, mpz_t y, mpz_t mod);

void mpz_mulm(mpz_t rop,mpz_t x, mpz_t y, mpz_t mod);

void mpz_mulm_ui(mpz_t rop, mpz_t x, unsigned long int y, mpz_t mod);

void mpz_addm_ui(mpz_t rop, mpz_t x, unsigned long int y, mpz_t mod);

void mpz_negm_ui(mpz_t rop, unsigned long int x, mpz_t mod);

void get_nu(mpz_t rop,  mpz_t d, mpz_t mod,  mpz_t u,  mpz_t v);

void init_point(struct point a);

void print_point(struct point a, char *name);

void get_ep(mpz_t rop, mpz_t d, mpz_t x, mpz_t y, mpz_t mod);

void get_xy(mpz_t x, mpz_t y, mpz_t u, mpz_t v, mpz_t d, mpz_t mod);

struct point rot_sum(struct point a, struct point b, struct curve_params c);

struct point double_point(struct point a, struct curve_params b);

struct point sum_points(struct point a, struct point b, struct curve_params params);

struct point to_affine(struct point proj, struct curve_params params);

struct point gen_mult_point(struct point a, struct curve_params b, mpz_t k);

int is_on_curve(struct point a, struct curve_params b);

struct point test_q(struct point a, struct curve_params b);

void test_next_prev(struct point a, struct curve_params b);

void test_random_orders(struct point a, struct curve_params b);

struct point neg(struct point a, struct curve_params b);

#endif//FUNCS_H
