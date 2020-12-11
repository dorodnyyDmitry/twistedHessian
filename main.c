#include "funcs.h"

int main()
{
    mpz_t b, res, x1, x2, mod, temp, ep, nu, d;

    mpz_init_set_str(b, par_b, 10);
    mpz_init_set_str(mod, par_mod, 10);

    mpz_init_set_str(temp, "54", 10);
    mpz_init_set_str(res, "8", 10);
    mpz_invert(temp, temp, mod);

    mpz_mulm(temp, temp, b, mod);// b/54
    mpz_negm(temp, temp, mod);
    mpz_negm(res, res, mod);
    mpz_addm(res, res, temp, mod);//res = -8 - b/54 ( everything to the left side)

    mpz_clear(temp);
    mpz_init_set_str(temp, "100", 10);

    mpz_negm(res, res, mod);
    mpz_addm(res, res, temp, mod);//res = 100 - c (D/4 = b/2^2 - ac)

    gmp_printf("D/4:     %Zd\n", res);

    mpz_clear(temp);
    mpz_init_set_str(temp, "1", 10);

    int discsquare = mpz_sqrtm(temp, res, mod); //temp = sqrt(D/4)

    printf("%d\n", discsquare);

    mpz_init_set_str(x1, "1", 10);
    mpz_init_set_str(x2, "1", 10);

    mpz_addm_ui(x1, temp, 10, mod); // x1 = 10 + temp (x1 = -b/2 + sqrt(D/4))

    mpz_negm(temp, temp, mod);

    mpz_addm_ui(x2, temp, 10, mod); // x2 = 10 - temp (x2 = -b/2 - sqrt(D/4))

    gmp_printf("x1     : %Zd\n", x1);
    gmp_printf("x2     : %Zd\n", x2);

    get_d(x1, mod);
    get_d(x2, mod);

    gmp_printf("dx1    : %Zd\n", x1);
    gmp_printf("dx2    : %Zd\n", x2);

    struct point my_point, new_point;
    struct curve_params my_params;

    mpz_inits(new_point.x, new_point.y, NULL);
    mpz_init_set_str(my_point.x, "30169251516072949329203362288258354412244051188206333399160894575841649853581", 10);
    mpz_init_set_str(my_point.y, "24928593895619211779139364461724456000389831107969176724087541868421424051322", 10);
    mpz_init_set_str(my_point.z, "1", 10);
    mpz_init_set_str(new_point.z, "1", 10);

    mpz_init_set_ui(my_params.a, 1);
    mpz_init_set_ui(my_params.d, 3);
    mpz_init_set_str(my_params.mod, par_mod, 10);

    get_xy(new_point.x, new_point.y, my_point.x, my_point.y, my_params.d, my_params.mod);

    mpz_mul_ui(my_params.d, my_params.d, 3);
    printf("%s\n", "From Weierstrass to Twisted Hessian:\n");
    gmp_printf("X      : %Zd\n", new_point.x);
    gmp_printf("Y      : %Zd\n", new_point.y);
    gmp_printf("Z      : %Zd\n", new_point.z);

    is_on_curve(new_point, my_params);

    struct point test_point1, test_point2, test_point3;

    mpz_inits(test_point3.x, test_point3.y, test_point3.z, NULL);

    // point T
    mpz_init_set_str(test_point1.x, "39631536720765161579801146192725376071535632056599520883371714951976620772810", 10);
    mpz_init_set_str(test_point1.y,"10312566352886215418716679103959354514031220653989860530756801665653060039489" ,10);
    mpz_init_set_ui(test_point1.z, 1);

    // point O = (0; -1; 1) - Neutral element
    mpz_init_set_str(test_point2.x,"0" ,10);
    mpz_init_set_str(test_point2.y, "115792089237316195423570985008687907853269984665640564039457584007913111864738", 10);
    mpz_init_set_ui(test_point2.z, 1);

    test_point3 = rot_sum(test_point1, test_point2, my_params);
    is_on_curve(test_point3, my_params);

    print_point(test_point3, "T + 0");

    if(mpz_cmp(test_point3.x , test_point1.x) ==0 && mpz_cmp(test_point3.y , test_point1.y) == 0 && mpz_cmp(test_point1.z , test_point3.z) == 0){
        printf("\n%s\n", "X + neutral == X - TRUE");
    }


    test_point1 = test_point3;
    print_point(test_point1, "shitt");

    test_point3 = rot_sum(test_point3, test_point3, my_params); //Tx2
    test_point3 = rot_sum(test_point3, test_point3, my_params); //Tx4
    test_point3 = rot_sum(test_point3, test_point3, my_params); //Tx8
    test_point3 = rot_sum(test_point3, test_point3, my_params); //Tx16

    test_point3 = to_affine(test_point3, my_params);
    print_point(test_point3, "T + ... + T x16");

    mpz_t k;
    struct point test;
    init_point(test);
    mpz_init_set_str(k, "16", 10);

    test = gen_mult_point(test_point1, my_params, k); //[16]T

    test = to_affine(test, my_params);
    print_point(test, "[16]T");

    is_on_curve(test, my_params);

    if(mpz_cmp(test_point3.x , test.x) ==0 && mpz_cmp(test_point3.y , test.y) == 0 && mpz_cmp(test.z , test_point3.z) == 0){
        printf("\n%s\n", "T * 16 == [16]T - TRUE");
    }

    struct point neutral_test, qpoint;
    mpz_inits(qpoint.x, qpoint.y, qpoint.z, NULL);

    mpz_init_set_str(test_point2.x, "93528601524582384978654134272795468222405403055717890280271688132874849008326", 10);
    mpz_init_set_str(test_point2.y, "14443324612566128911211262381388707474030458136470034119105598903952521080679", 10);
    mpz_init_set_ui(test_point2.z, 1);

    qpoint = test_q(test_point2, my_params);

    mpz_init_set_ui(neutral_test.x, 0);
    mpz_init_set_str(neutral_test.y, "115792089237316195423570985008687907853269984665640564039457584007913111864738", 10);
    mpz_init_set_ui(neutral_test.z, 1);

    if(mpz_cmp(neutral_test.x, qpoint.x) == 0 && mpz_cmp(neutral_test.y, qpoint.y) ==0 && mpz_cmp(neutral_test.z, qpoint.z) ==0){
        printf("\n%s", "[q]P == neutral element - TRUE\n");
    }

    test_next_prev(test_point2, my_params);

    test_random_orders(test_point2, my_params);

    return 0;
}
