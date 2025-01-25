/**********************************************************
 *                                                        *
 *    gcc -g -lm -I. Frobenius.c libgmp.a -o Frobenius    *
 *                                                        *
 **********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>

mpz_t *calc_s_t(mpz_t a, mpz_t T, mpz_t n)
{
    static mpz_t out[2];
    /**************************************
     * out[0] = s_n                       *
     * out[1] = t_n                       *
     **************************************/
    mpz_t t0, s0, h1, h2;
    mpz_set_ui(out[0], 1);
    mpz_set(out[1], T);
    mpz_init(s0);
    mpz_init(t0);
    mpz_init(h1);
    mpz_init(h2);

    for (int j = mpz_sizeinbase(n, 2) - 2; j > -1; j--)
    {
        mpz_set(s0, out[0]);
        mpz_set(t0, out[1]);

        if (mpz_tstbit(n, j) == 1)
        {
            mpz_mul(h1, s0, a);
            mpz_mul_ui(h2, t0, 2);
            mpz_add(h1, h1, h2);
            mpz_mul(out[0], h1, s0);
            mpz_mod(out[0], out[0], n);

            mpz_add(h1, t0, s0);
            mpz_sub(h2, t0, s0);
            mpz_mul(out[1], h1, h2);
            mpz_mod(out[1], out[1], n);

            mpz_set(s0, out[0]);
            mpz_set(t0, out[1]);

            mpz_mul(h1, a, s0);
            mpz_mul(h2, T, s0);
            mpz_add(h2, h1, h2);
            mpz_add(out[0], h2, t0);
            mpz_mod(out[0], out[0], n);

            mpz_mul(h1, t0, T);
            mpz_sub(out[1], h1, s0);
            mpz_mod(out[1], out[1], n);
        }
        else
        {
            mpz_mul(h1, s0, a);
            mpz_mul_ui(h2, t0, 2);
            mpz_add(h1, h1, h2);
            mpz_mul(out[0], h1, s0);
            mpz_mod(out[0], out[0], n);

            mpz_add(h1, t0, s0);
            mpz_sub(h2, t0, s0);
            mpz_mul(out[1], h1, h2);
            mpz_mod(out[1], out[1], n);
        }
    }

    mpz_clear(s0);
    mpz_clear(t0);
    mpz_clear(h1);
    mpz_clear(h2);

    return out;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("\nUsage: ./Frobenius <lower bound> <upper bound>\n");
        exit(1);
    }

    mpz_t lower_bound, upper_bound, n, n_1, n_1_2, a, aa_4, abs_aa_4, T, Q, P, PQ, gcd;
    mpz_t *s_t;

    int j;
    mpz_init(n);
    mpz_init(n_1);
    mpz_init(n_1_2);
    mpz_init(a);
    mpz_init(aa_4);
    mpz_init(abs_aa_4);
    mpz_init(T);
    mpz_init(Q);
    mpz_init(P);
    mpz_init(PQ);
    mpz_init(gcd);

    mpz_init_set_str(lower_bound, argv[1], 10);
    mpz_init_set_str(upper_bound, argv[2], 10);

    if (mpz_cmp_ui(lower_bound, 15) == -1)
    {
        mpz_set_ui(lower_bound, 15);
    }
    if (mpz_even_p(lower_bound))
    {
        mpz_add_ui(lower_bound, lower_bound, 1);
    }

    for (mpz_set(n, lower_bound); mpz_cmp(n, upper_bound) <= 0; mpz_add_ui(n, n, 2))
    {
        if (!mpz_probab_prime_p(n, 1) && !mpz_perfect_square_p(n))
        {
            mpz_set_ui(a, 1);
            mpz_sub_ui(aa_4, a, 4); // aa_4 = a^2 - 4
            mpz_gcd(gcd, aa_4, n);
            while ((mpz_cmp_ui(gcd, 1) == 0) && (mpz_jacobi(aa_4, n) != -1))
            {
                mpz_add_ui(a, a, 2);
                mpz_mul(aa_4, a, a);
                mpz_sub_ui(aa_4, aa_4, 4);
                mpz_gcd(gcd, aa_4, n);
            }
            mpz_abs(abs_aa_4, aa_4);

            mpz_set_ui(T, 1);
            mpz_add(Q, T, a);
            mpz_mul(Q, Q, T);
            mpz_add_ui(Q, Q, 1);
            mpz_mod(Q, Q, n); // Q = T^2 + aT + 1 mod n

            mpz_mul_ui(P, T, 2);
            mpz_add(P, P, a); // P = 2T + a

            mpz_mul(PQ, P, Q);
            mpz_gcd(gcd, PQ, n);
            while ((mpz_cmp_ui(gcd, 1) == 0) && ((mpz_jacobi(Q, n) != -1) || (mpz_cmp(abs_aa_4, Q) == 0)))
            {
                mpz_add_ui(T, T, 1); // T = T + 1
                mpz_add(Q, T, a);
                mpz_mul(Q, Q, T);
                mpz_add_ui(Q, Q, 1);
                mpz_mod(Q, Q, n); // Q = T^2 + aT + 1 mod n
                mpz_mul_ui(P, T, 2);
                mpz_add(P, P, a); // P = 2T + a
                mpz_mul(PQ, P, Q);
                mpz_gcd(gcd, PQ, n);
            }
            mpz_mul(PQ, PQ, aa_4);
            mpz_gcd(gcd, PQ, n); // gcd(a^2-4, P, Q)

            if (mpz_cmp_ui(gcd, 1) == 0)
            {
                mpz_sub_ui(n_1, n, 1);
                mpz_div_ui(n_1_2, n_1, 2);
                mpz_powm(Q, Q, n_1_2, n);
                if (mpz_cmp(Q, n_1) == 0) // Q^((n-1)/2) = -1 mod n
                {
                    gmp_printf("Q %Zd\n", n);
                }
                s_t = calc_s_t(a, T, n);
                if (mpz_cmp(s_t[0], n_1) == 0) // s_n = -1 mod n
                {
                    gmp_printf("sn %Zd\n", n);
                }
                mpz_add(s_t[0], a, T);
                mpz_sub(s_t[0], s_t[0], s_t[1]);
                mpz_mod(s_t[0], s_t[0], n);
                if (mpz_cmp_ui(s_t[0], 0) == 0) // t_n = a + T mod n
                {
                    gmp_printf("tn %Zd\n", n);
                }
            }
        }
    }

    fflush(stdout);
    exit(0);
}
