#include "funcs.h"

int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p) {
    mpz_t w, n_inv, y;
    unsigned int i, s;


    if(mpz_divisible_p(n, p)) {         /* Is n a multiple of p?            */
        mpz_set_ui(q, 0);               /* Yes, then the square root is 0.  */
        return 1;                       /* (special case, not caught        */
    }                                   /* otherwise)                       */
    if(mpz_legendre(n, p) != 1)         /* Not a quadratic residue?         */
        return 0;                       /* No, so return error              */
    if(mpz_tstbit(p, 1) == 1) {         /* p = 3 (mod 4) ?                  */
        mpz_set(q, p);
        mpz_add_ui(q, q, 1);
        mpz_fdiv_q_2exp(q, q, 2);
        mpz_powm(q, n, q, p);           /* q = n ^ ((p+1) / 4) (mod p)      */
        return 1;
    }

    mpz_set(q, p);
    mpz_sub_ui(q, q, 1);                /* q = p-1                          */
    s = 0;                              /* Factor out 2^s from q            */
    while(mpz_tstbit(q, s) == 0) s++;
    mpz_fdiv_q_2exp(q, q, s);           /* q = q / 2^s                      */
    mpz_set_ui(w, 2);                   /* Search for a non-residue mod p   */
    while(mpz_legendre(w, p) != -1)     /* by picking the first w such that */
        mpz_add_ui(w, w, 1);            /* (w/p) is -1                      */
    mpz_powm(w, w, q, p);               /* w = w^q (mod p)                  */
    mpz_add_ui(q, q, 1);
    mpz_fdiv_q_2exp(q, q, 1);           /* q = (q+1) / 2                    */
    mpz_powm(q, n, q, p);               /* q = n^q (mod p)                  */
    mpz_invert(n_inv, n, p);
    for(;;) {
        mpz_powm_ui(y, q, 2, p);        /* y = q^2 (mod p)                  */
        mpz_mul(y, y, n_inv);
        mpz_mod(y, y, p);               /* y = y * n^-1 (mod p)             */
        i = 0;
        while(mpz_cmp_ui(y, 1) != 0) {
            i++;
            mpz_powm_ui(y, y, 2, p);    /*  y = y ^ 2 (mod p)               */
        }
        if(i == 0) {                    /* q^2 * n^-1 = 1 (mod p), return   */

            return 1;
        }
        if(s-i == 1) {                  /* In case the exponent to w is 1,  */
            mpz_mul(q, q, w);           /* Don't bother exponentiating      */
        } else {
            mpz_powm_ui(y, w, 1 << (s-i-1), p);
            mpz_mul(q, q, y);
        }
        mpz_mod(q, q, p);               /* r = r * w^(2^(s-i-1)) (mod p)    */
    }
}

void get_d(mpz_t x, mpz_t mod){
    mpz_add_ui(x, x, 8);
    mpz_mod(x, x, mod);
    mpz_invert(x, x, mod);
    mpz_mul_ui(x, x, 105);
    mpz_mod(x, x, mod);
}

void mpz_negm(mpz_t rop, mpz_t x, mpz_t mod){
    mpz_neg(rop, x);
    mpz_mod(rop, rop, mod);
}

void mpz_addm(mpz_t rop,mpz_t x, mpz_t y, mpz_t mod){
    mpz_add(rop, x, y);
    mpz_mod(rop, rop, mod);
}

void mpz_mulm(mpz_t rop,mpz_t x, mpz_t y, mpz_t mod){
    mpz_mul(rop, x, y);
    mpz_mod(rop, rop, mod);
}

void mpz_mulm_ui(mpz_t rop, mpz_t x, unsigned long int y, mpz_t mod){
    mpz_mul_ui(rop, x, y);
    mpz_mod(rop, rop, mod);
}

void mpz_addm_ui(mpz_t rop, mpz_t x, unsigned long int y, mpz_t mod){
    mpz_add_ui(rop, x, y);
    mpz_mod(rop, rop, mod);
}

void mpz_negm_ui(mpz_t rop, unsigned long int x, mpz_t mod){
    mpz_t temp;
    mpz_init_set_ui(temp, x);
    mpz_negm(rop, temp, mod);
    mpz_clear(temp);
}

void get_nu(mpz_t rop,  mpz_t d, mpz_t mod,  mpz_t u,  mpz_t v){
    mpz_t temp1, temp2, temp3, d2, d3;
    mpz_inits(temp1, temp2, temp3, d2, d3, NULL);

    mpz_powm_ui(d2, d, 2, mod);
    mpz_powm_ui(d3, d, 3, mod);

    mpz_mulm_ui(temp1, d, 3, mod);
    mpz_mulm(temp1, temp1, u, mod);
    mpz_negm(temp1, temp1, mod); // -3Du

    mpz_negm_ui(temp2, 36, mod); // -36

    mpz_addm(temp1, temp1, temp2, mod); //-36 -3Du

    mpz_mulm_ui(temp2, d3, 9, mod); // 9D^3

    mpz_addm(temp1, temp1, temp2, mod); // 9D^3 -36 - 3Du

    mpz_addm(temp1, temp1, v, mod); //v + 9D^3 - 3Du - 36

    mpz_negm_ui(temp2, 1, mod);
    mpz_addm(temp2, temp2, d3, mod); // D^3 - 1

    mpz_mulm(temp1, temp1, temp2, mod);
    mpz_mulm_ui(temp1, temp1, 6, mod); // 6(D^3 - 1)(v + 9D^3 - 3Du - 36)

    mpz_mulm_ui(temp2, d3, 3, mod); // 3D^3

    mpz_mulm(temp3, d, u, mod);
    mpz_negm(temp3, temp3, mod); // -Du

    mpz_addm(temp2, temp2, temp3, mod); // 3D^3 - Du

    mpz_negm_ui(temp3, 12, mod);
    mpz_addm(temp2, temp2, temp3, mod); // 3D^3 - Du - 12
    mpz_powm_ui(temp2, temp2, 3, mod); // (3D^3 - Du - 12)^3

    mpz_mulm_ui(temp3, d2, 9, mod); // 9D^2
    mpz_addm(temp3, temp3, u, mod); // u + 9D^2
    mpz_powm_ui(temp3, temp3, 3, mod); // (u + 9D^2)^3

    mpz_addm(temp2, temp2, temp3, mod); // (u + 9D^2)^3 + (3D^3 - Du - 12)^3

    mpz_invert(temp2, temp2, mod); // 1/((u + 9D^2)^3 + (3D^3 - Du - 12)^3)
    mpz_mulm(rop, temp1, temp2, mod);

    mpz_clears(temp1, temp2, temp3, d2, d3, NULL);

}

void init_point(struct point a){
    mpz_inits(a.x, a.y, a.z, NULL);
}

void print_point(struct point a, char *name){
    gmp_printf("%s x: %Zd\n", name, a.x);
    gmp_printf("%s y: %Zd\n", name, a.y);
    gmp_printf("%s z: %Zd\n", name, a.z);
}

void get_ep(mpz_t rop, mpz_t d, mpz_t x, mpz_t y, mpz_t mod){
    mpz_t temp1, temp2, d3;
    mpz_inits(temp1, temp2, d3, NULL);

    mpz_powm_ui(d3, d, 3, mod);

    mpz_negm_ui(temp1, 1, mod);
    mpz_addm(temp1, temp1, d3, mod); // D^3 - 1
    mpz_mulm_ui(temp1, temp1, 12, mod); // 12(D^3 - 1)

    mpz_mulm(temp2, d, x, mod);
    mpz_addm(temp2, temp2, y, mod);
    mpz_addm_ui(temp2, temp2, 1, mod); // Dx + y + 1
    mpz_invert(temp2, temp2, mod); // 1/(Dx + y + 1)

    mpz_mulm(rop, temp1, temp2, mod);

    mpz_clears(temp1, temp2, d3, NULL);

}

void get_xy(mpz_t x, mpz_t y, mpz_t u, mpz_t v, mpz_t d, mpz_t mod){
    mpz_t temp1, temp2, temp3, d2, d3, nu;
    mpz_inits(temp1, temp2, temp3, d2, d3, nu, NULL);

    get_nu(nu, d, mod, u, v);

    mpz_powm_ui(d2, d, 2, mod);
    mpz_powm_ui(d3, d, 3, mod);

    //X
    mpz_mulm_ui(temp1, d2, 9, mod); // 9D^2
    mpz_addm(temp1, temp1, u, mod); // 9D^2 + u
    mpz_mulm(x, temp1, nu, mod); // nu*(9D^2 + u)


    //Y
    mpz_mulm(temp2, d, u, mod); // D*x

    mpz_negm(temp2, temp2, mod); // -D*x
    mpz_negm_ui(temp3, 12, mod); // -12
    mpz_addm(temp2, temp2, temp3, mod); // -D*x - 12

    mpz_mulm_ui(temp3, d3, 3, mod); // 3D^3
    mpz_addm(temp1, temp2, temp3, mod); // 3*D^3 -D*x - 12
    mpz_mulm(temp2, temp1, nu, mod); // nu*(3*D^3 -D*x - 12)

    mpz_negm_ui(temp3, 1, mod); // -1
    mpz_addm(y, temp2, temp3, mod); // - 1 + nu*(3*D^3 -D*x - 12)

    mpz_clears(temp1, temp2, temp3, d2, d3, nu, NULL);
}

struct point rot_sum(struct point a, struct point b, struct curve_params c){
    struct point res, tempa, tempb;
    mpz_t temp1, temp2, temp3, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F;
    mpz_inits(temp1, temp2,temp3,  temp_A, temp_B, temp_C, temp_D, temp_E, temp_F, res.x, res.y, res.z,tempa.x, tempa.y, tempa.z,tempb.x, tempb.y, tempb.z, NULL);
    mpz_negm_ui(temp1, 1, c.mod);

    mpz_mulm(temp_A, a.x, b.z, c.mod);
    mpz_mulm(temp_B, a.z, b.z, c.mod);
    mpz_mulm(temp_C, a.y, b.x, c.mod);
    mpz_mulm(temp_D, a.y, b.y, c.mod);
    mpz_mulm(temp_E, a.z, b.y, c.mod);
    mpz_mulm(temp3, a.x, b.x, c.mod);
    mpz_mulm(temp_F, temp3, c.a, c.mod);

    mpz_mulm(temp1, temp_A, temp_B, c.mod);
    mpz_mulm(temp2, temp_C, temp_D, c.mod);
    mpz_negm(temp3, temp2, c.mod);
    mpz_addm(res.x, temp1, temp3, c.mod); //X3 = AB - CD

    mpz_mulm(temp1, temp_D, temp_E, c.mod);
    mpz_mulm(temp2, temp_F, temp_A, c.mod);
    mpz_negm(temp3, temp2, c.mod);
    mpz_addm(res.y, temp1, temp3, c.mod); //Y3 = DE - FA

    mpz_mulm(temp1, temp_F, temp_C, c.mod);
    mpz_mulm(temp2, temp_B, temp_E, c.mod);
    mpz_negm(temp3, temp2, c.mod);
    mpz_addm(res.z, temp1, temp3, c.mod);

    mpz_clears(temp1, temp2,temp3, temp_A, temp_B, temp_C, temp_D, temp_E, temp_F, NULL);

    return res;
}

struct point double_point(struct point a, struct curve_params b){
    struct point res;
    mpz_t temp1, temp2, temp_P, temp_2P, temp_S, temp_A, temp_C, temp_D, temp_E, temp_F, temp_G, temp_minustwo;
    mpz_inits(temp1, temp2, temp_P, temp_2P, temp_S, temp_A, temp_C, temp_D, temp_E,temp_F, temp_G, temp_minustwo, res.x, res.y, res.z, NULL);

    struct point temp;
    mpz_inits(temp.x, temp.y, temp.z, NULL);
    temp = a;
    res = rot_sum(a, temp, b);

    /*mpz_addm(temp_P, a.y, a.z, b.mod); //P = Y*Z

    mpz_mulm_ui(temp_2P, temp_P, 2, b.mod); //2*P

    mpz_addm(temp_S, a.y, a.z, b.mod); // S = Y*Z

    mpz_powm_ui(temp1, temp_S, 2, b.mod);
    mpz_negm(temp2, temp_P, b.mod);
    mpz_addm(temp_A, temp1, temp2, b.mod); // A = S^2 - P

    mpz_negm(temp1, temp_2P, b.mod);
    mpz_addm(temp1, temp_A, temp1, b.mod);
    mpz_mulm(temp_C, temp1, temp_S, b.mod); // C = (A - 2*P)*S

    mpz_negm(temp1, a.y, b.mod);
    mpz_addm(temp1, temp1, a.z, b.mod);
    mpz_mulm(temp_D, temp1, temp_A, b.mod); // D = A*(Z - Y)

    mpz_mulm(temp1, a.x, temp_2P, b.mod);
    mpz_mulm(temp1, temp1, b.d, b.mod);
    mpz_negm(temp1, temp1, b.mod);
    mpz_mulm_ui(temp2, temp_C, 3, b.mod);
    mpz_addm(temp_E, temp1, temp2, b.mod); // E = 3C - d*X*2P

    mpz_negm_ui(temp_minustwo, 2, b.mod); // minustwo = -2 (???)

    mpz_mulm(temp1, a.x, b.d, b.mod);
    mpz_mulm(res.x, temp_minustwo, temp1, b.mod); // X3 = minustwo*X*D

    mpz_negm(temp1, temp_E, b.mod);
    mpz_addm(temp1, temp_D, temp_E, b.mod);
    mpz_mulm(res.y, temp1, a.z, b.mod); // Y3 = (D - E)*Z1

    mpz_addm(temp1, temp_D, temp_E, b.mod);
    mpz_mulm(res.z, temp1, a.y, b.mod); // Z3 = (D + E)*Y1*/


    /*mpz_powm_ui(temp_D, a.x, 3, b.mod); // D = X^3
    mpz_powm_ui(temp_E, a.y, 3, b.mod); // E = Y^3
    mpz_powm_ui(temp_F, a.z, 3, b.mod); // F = Z^3
    mpz_mulm(temp_G, temp_D, b.a, b.mod); // G = a*D

    gmp_printf("D :%Zd\n", temp_D);
    gmp_printf("E :%Zd\n", temp_E);
    gmp_printf("F :%Zd\n", temp_F);
    gmp_printf("G :%Zd\n", temp_G);

    mpz_negm(temp1, temp_F, b.mod);
    mpz_addm(temp1, temp1, temp_E, b.mod);
    mpz_mulm(res.x, temp1, a.x, b.mod); // new X = X(E - F)

    //gmp_printf("X :%Zd\n", res.x);

    mpz_negm(temp1, temp_E, b.mod);
    mpz_addm(temp1, temp1, temp_G, b.mod);
    mpz_mulm(res.y, temp1, a.z, b.mod); // new Y = Z(G - E)

    //gmp_printf("Y :%Zd\n", res.y);

    mpz_negm(temp1, temp_G, b.mod);
    mpz_addm(temp1, temp1, temp_F, b.mod);
    mpz_mulm(res.z, temp1, a.y, b.mod); // new Z = Y(F - G)*/

    //gmp_printf("Z :%Zd\n", res.z);

    mpz_clears(temp1, temp2, temp_P, temp_2P, temp_S, temp_A, temp_C, temp_D, temp_E,temp_F, temp_G,  temp_minustwo, NULL);
    return res;
}

struct point sum_points(struct point a, struct point b, struct curve_params params){
    mpz_t temp1, temp2;
    struct point res;

    /* X3 */

    mpz_inits(temp1, temp2, res.x, res.y, res.z, NULL);

    mpz_powm_ui(temp1, a.x, 2, params.mod);   // temp1 = X1^2
    mpz_mulm(temp1, temp1, b.y, params.mod);  // temp1 = X1^2 * Y2
    mpz_mulm(temp1, temp1, b.z, params.mod);  // temp1 = X1^2 * Y2 * Z2

    mpz_powm_ui(temp2, b.x, 2, params.mod);   // temp2 = X2^2
    mpz_mulm(temp2, temp2, a.y, params.mod);  // temp2 = X2^2 * Y1
    mpz_mulm(temp2, temp2, a.z, params.mod);  // temp2 = X2^2 * Y1 * Z1
    mpz_negm(temp2, temp2, params.mod);       // temp2 = -( X2^2 * Y1 * Z1 )

    mpz_addm(res.x, temp1, temp2, params.mod);// res.x = X1^2 * Y2 * Z2 -( X2^2 * Y1 * Z1 )

    /* Y3 */

    mpz_powm_ui(temp1, a.z, 2, params.mod);   // temp1 = Z1^2
    mpz_mulm(temp1, temp1, b.y, params.mod);  // temp1 = Z1^2 * Y2
    mpz_mulm(temp1, temp1, b.x, params.mod);  // temp1 = Z1^2 * Y2 * X2

    mpz_powm_ui(temp2, b.z, 2, params.mod);   // temp2 = Z2^2
    mpz_mulm(temp2, temp2, a.y, params.mod);  // temp2 = Z2^2 * Y1
    mpz_mulm(temp2, temp2, a.x, params.mod);  // temp2 = Z2^2 * Y1 * X1
    mpz_negm(temp2, temp2, params.mod);       // temp2 = -( Z2^2 * Y1 * X1 )

    mpz_addm(res.y, temp1, temp2, params.mod);// res.y = Z1^2 * Y2 * X2 -( Z2^2 * Y1 * X1 )

    /* Z3 */

    mpz_powm_ui(temp1, a.y, 2, params.mod);   // temp1 = Y1^2
    mpz_mulm(temp1, temp1, b.z, params.mod);  // temp1 = Y1^2 * Z2
    mpz_mulm(temp1, temp1, b.x, params.mod);  // temp1 = Y1^2 * Z2 * X2

    mpz_powm_ui(temp2, b.y, 2, params.mod);   // temp2 = Y2^2
    mpz_mulm(temp2, temp2, a.z, params.mod);  // temp2 = Y2^2 * Z1
    mpz_mulm(temp2, temp2, a.x, params.mod);  // temp2 = Y2^2 * Z1 * X1
    mpz_negm(temp2, temp2, params.mod);       // temp2 = -( Y2^2 * Z1 * X1 )

    mpz_addm(res.z, temp1, temp2, params.mod);// res.y = Y1^2 * Z2 * X2 -( Y2^2 * Z1 * X1 )

    mpz_clears(temp1, temp2, NULL);

    return(res);
}

struct point to_affine(struct point proj, struct curve_params params){
    mpz_t temp;
    mpz_init(temp);

    struct point res;
    mpz_inits(res.x, res.y, res.z, NULL);

    mpz_invert(temp, proj.z, params.mod);

    mpz_mulm(res.x, proj.x, temp, params.mod);
    mpz_mulm(res.y, proj.y, temp, params.mod);
    mpz_mulm(res.z, proj.z, temp, params.mod);

    mpz_clear(temp);
    return res;
}

struct point gen_mult_point(struct point a, struct curve_params b, mpz_t k){
    struct point q, temp1, temp2, r;
    mpz_inits(temp1.x, temp1.y, temp1.z,temp2.x, temp2.y, temp2.z, r.x, r.y, r.z, NULL);
    r = a;
    mpz_init_set_ui(q.x, 0);
    mpz_init_set_str(q.y, "115792089237316195423570985008687907853269984665640564039457584007913111864738", 10);
    mpz_init_set_ui(q.z, 1);

    int i = strlen(mpz_get_str(NULL, 2, k))-1;
    //printf("%s\n", mpz_get_str(NULL, 2, k));

    for(;i >= 0; i--){ //MOntgomery algorythmS


        if(mpz_tstbit(k, i) == 1){
            q = sum_points(q, r, b);
            r = rot_sum(r, r, b);
        }
        else{
            r = sum_points(r, q, b);
            q = rot_sum(q, q, b);
        }
        /*printf("%d\n", i);
        temp1 = to_affine(r, b);
        temp2 = to_affine(q, b);
        print_point(temp1, "r");
        print_point(temp2, "q");*/
    }
    mpz_clears(temp1.x, temp1.y, temp1.z,temp2.x, temp2.y, temp2.z, r.x, r.y, r.z, NULL);
    return q;
}

int is_on_curve(struct point a, struct curve_params b){
    mpz_t temp1, temp2;
    mpz_inits(temp1, temp2, NULL);

    mpz_powm_ui(temp1, a.x, 3, b.mod);

    mpz_mulm(temp1, temp1, b.a, b.mod);// a*X^3

    mpz_powm_ui(temp2, a.y, 3, b.mod);
    mpz_addm(temp1, temp1, temp2, b.mod);// a*X^3 + Y^3

    mpz_powm_ui(temp2, a.z, 3, b.mod);
    mpz_addm(temp1, temp1, temp2, b.mod);// a*X^3 + Y^3 + z^3

    mpz_mulm(temp2, a.x, a.y, b.mod);
    mpz_mulm(temp2, temp2, a.z, b.mod);
    mpz_mulm(temp2, temp2, b.d, b.mod); // d*X*Y*Z


    int res = mpz_cmp(temp1, temp2);

    if(res == 0){
        printf("\n%s\n", "Point is on the curve");

    }
    else{
        printf("\n%s\n", "Point is NOT on the curve");
    }
    mpz_clears(temp1, temp2, NULL);
    return res;
}

struct point test_q(struct point a, struct curve_params b){
    mpz_t q;
    struct point res;

    mpz_inits(res.x, res.y, res.z, NULL);
    mpz_init_set_str(q, par_q, 10);

    res = gen_mult_point(a, b, q);

    res = to_affine(res, b);

    print_point(res, "[q]P");

    mpz_clear(q);
    return res;
}

void test_next_prev(struct point a, struct curve_params b){
    mpz_t q, q1, q2;
    struct point next, prev, test;
    mpz_init_set_str(q, par_q, 10);
    mpz_init_set_str(q1, par_q, 10);
    mpz_init_set_str(q2, par_q, 10);

    mpz_sub_ui(q1, q1, 1);
    mpz_add_ui(q2, q2, 1);

    prev = gen_mult_point(a, b, q1);
    next = gen_mult_point(a, b, q2);

    test = a;

    a = to_affine(a, b);
    print_point(a, "P");

    next = to_affine(next, b);
    print_point(next, "[q + 1]P");

    prev = to_affine(prev, b);
    print_point(prev, "[q - 1]P");

    test = neg(test, b);
    test = to_affine(test, b);
    print_point(test, "neg P");

    printf("%s\n", "[q + 1]P == P, [q - 1]P == -P - TRUE");

}

void test_random_orders(struct point a, struct curve_params b){
    mpz_t k1, k2, k1k2;
    mpz_inits(k1, k2, k1k2, NULL);

    struct point point_k1, point_k2, test, point_k1k2;
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);

    mpz_urandomm(k1, rstate, b.mod);
    mpz_urandomm(k2, rstate, b.mod);

    mpz_addm(k1k2, k1, k2, b.mod);

    point_k1 = gen_mult_point(a, b, k1);
    point_k2 = gen_mult_point(a, b, k2);
    point_k1k2 = gen_mult_point(a, b, k1k2);

    test = rot_sum(point_k1, point_k2, b);
    test = to_affine(test, b);
    print_point(test, "[k1]P + [k2]P");

    point_k1k2 = to_affine(point_k1k2, b);

    print_point(point_k1k2, "[k1 + k2]P   ");

    printf("%s\n", "[k1]P + [k2]P == [k1 + k2]P - TRUE");
}

struct point neg(struct point a, struct curve_params b){
    struct point res;
    mpz_inits(res.x, res.y, res.z, NULL);
    /*mpz_negm(res.x, a.x, b.mod);
    mpz_negm(res.y, a.y, b.mod);
    mpz_negm(res.z, a.z, b.mod);*/
    mpz_swap(a.y, a.z);

    return a;
}
