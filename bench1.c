#include "arb_hypgeom.h"
#include "flint/profiler.h"

static double log2_bin_uiui_fast(ulong n, ulong k)
{
    static const float htab[] = {0.2007, 0.3374, 0.4490, 0.5437, 0.6254, 0.6963,
        0.7580, 0.8114, 0.8572, 0.8961, 0.9285, 0.9545, 0.9746,
        0.9888, 0.9973, 1.0};

    if (k == 0 || k >= n)
        return 0;
    if (k > n / 2)
        k = n - k;
    k = (32.0 * k) / n;
    return n * htab[FLINT_MIN(k, 15)];
}

void
arb_hypgeom_legendre_p_ui_deriv_bound(mag_t dp, mag_t dp2, ulong n, const arb_t x, const arb_t x2sub1)
{
    mag_t t, u, xm, x2sub1m;

    mag_init(t);
    mag_init(u);
    mag_init(xm);
    mag_init(x2sub1m);

    arb_get_mag(xm, x);
    arb_get_mag_lower(x2sub1m, x2sub1);

    /* |P'(x)| <= min(n(n+1)/2, n/sqrt(1-x^2)) */
    mag_set_ui(t, n);
    mag_add_ui(t, t, 1);
    mag_mul_2exp_si(t, t, -1);
    mag_rsqrt(u, x2sub1m);
    mag_max(t, t, u);
    mag_mul_ui(dp, t, n);

    /* |P''(x)| <= min((n+2)(n+1)n(n-1)/8, (2x|P'(x)| + n(n+1))/(1-x^2)) */
    mag_mul(dp2, dp, xm);
    mag_mul_2exp_si(dp2, dp2, 1);

    mag_set_ui(t, n);
    mag_add_ui(t, t, 1);
    mag_mul_ui(t, t, n);
    mag_add(dp2, dp2, t);
    mag_div(dp2, dp2, x2sub1m);

    mag_set_ui(t, n);
    mag_add_ui(t, t, 2);
    mag_set_ui(u, n);
    mag_add_ui(u, u, 1);
    mag_mul(t, t, u);
    mag_mul_ui(t, t, n);
    mag_mul_ui(t, t, n - 1);
    mag_mul_2exp_si(t, t, -3);
    mag_min(dp2, dp2, t);


    mag_clear(t);
    mag_clear(u);
    mag_clear(xm);
    mag_clear(x2sub1m);
}

/* copy of the code, but with basecase disabled */
void
arb_hypgeom_legendre_p_ui_norec(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)
{
    arb_t xsub1, x2sub1;
    double xx, xxsub1, cancellation_zero, cancellation_one;
    double cost_zero, cost_one, cost_asymp;
    double log2x, log2u, tolerance, asymp_error;
    double yy, log2nsy, log2k, size;
    slong wp;
    slong d, k, K_zero, K_one, K_asymp;
    int basecase_ok;

    if (!arb_is_finite(x) || n > UWORD_MAX / 4)
    {
        if (res != NULL)
            arb_indeterminate(res);
        if (res_prime != NULL)
            arb_indeterminate(res_prime);
        return;
    }

    if (arf_sgn(arb_midref(x)) < 0)
    {
        arb_t t;
        arb_init(t);
        arb_neg(t, x);
        arb_hypgeom_legendre_p_ui(res, res_prime, n, t, prec);
        if (n % 2 == 1 && res != NULL)
            arb_neg(res, res);
        if (n % 2 == 0 && res_prime != NULL)
            arb_neg(res_prime, res_prime);
        arb_clear(t);
        return;
    }

    if (arb_is_one(x) && n < UWORD_MAX)
    {
        if (res != NULL)
            arb_set(res, x);
        if (res_prime != NULL)
        {
            arb_set_ui(res_prime, n);
            arb_mul_ui(res_prime, res_prime, n + 1, prec);
            arb_mul_2exp_si(res_prime, res_prime, -1);
        }
        return;
    }

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    if (n == 1)
    {
        if (res != NULL) arb_set_round(res, x, prec);
        if (res_prime != NULL) arb_one(res_prime);
        return;
    }

    xx = arf_get_d(arb_midref(x), ARF_RND_UP);

    /* Use basecase recurrence? */
    /* The following tests are not very elegant, and not completely accurate
       either, but they are fast in the common case. */
    if (res_prime != NULL)
    {
        basecase_ok = ((xx < 0.999999 && n < 10 && prec < 2000) ||
                       (xx < 0.999999 && n < 50 && prec < 1000) ||
                       (xx < 0.9999 && n < 100 && prec < 1000) ||
                       (xx < 0.999  && n < 350 && prec < 1000) ||
                       (xx < 0.9 && n < 400 && prec < 1000))
                   && ((xx > 0.00001 && n < 10 && prec < 2000) ||
                       (xx > 0.00001 && n < 60 && prec < 1000) ||
                       (xx > 0.01 && n < 200 && prec < 1000) ||
                       (xx > 0.1 && n < 400 && prec < 1000));
    }
    else if (prec < 500)
    {
        basecase_ok = ((xx < 0.999999 && n < 20) ||
                       (xx < 0.999 && n < 60) ||
                       (xx < 0.9 && n < 100))
                   && ((xx > 0.00001 && n < 20) ||
                       (xx > 0.01 && n < 60) ||
                       (xx > 0.1 && n < 100));
    }
    else
    {
        basecase_ok = 0;
    }

    basecase_ok = 0;

    if (basecase_ok)
    {
        mag_t t;
        mag_init(t);
        arb_get_mag(t, x);
        if (mag_cmp_2exp_si(t, 0) >= 0)
            basecase_ok = 0;
        mag_clear(t);
    }

    if (basecase_ok)
    {
        arb_hypgeom_legendre_p_ui_rec(res, res_prime, n, x, prec);
        return;
    }

    arb_init(xsub1);
    arb_init(x2sub1);

    arb_sub_ui(xsub1, x, 1, prec + 10);

    arb_mul(x2sub1, x, x, 2 * prec);
    arb_sub_ui(x2sub1, x2sub1, 1, prec + 10);
    arb_neg(x2sub1, x2sub1);

    /* use series at 1 unless |x| < 1-eps */
    if (!arb_is_negative(xsub1) ||
        arf_cmp_d(arb_midref(xsub1), ldexp(1.0, -2 * FLINT_BIT_COUNT(n))) > 0)
    {
        if (arf_cmp_d(arb_midref(xsub1), 2.0) >= 0)
        {
            if (n < 10000.0 * prec && n < UWORD_MAX / 4)
                K_one = n + 1;
            else
                K_one = 1;
        }
        else  /* check for early convergence */
        {
            xxsub1 = arf_get_d(arb_midref(xsub1), ARF_RND_UP);
            log2u = log(fabs(xxsub1) * 0.5) * 1.44269504088896;
            if (log2u < -30)
                log2u = arf_abs_bound_lt_2exp_si(arb_midref(xsub1)) - 1.0;

            K_one = n + 1;
            K_one = FLINT_MIN(K_one, 100000.0 * prec);
            K_one = FLINT_MIN(K_one, UWORD_MAX * 0.25);

            size = 0.0;

            if (n * (2.0 + log2u) < -prec)
            {
                for (k = 1; k < K_one; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, k)
                        + log2_bin_uiui_fast(n + k, k) + k * log2u;

                    if (size < -prec)
                    {
                        K_one = k;
                        break;
                    }
                }
            }
        }

        arb_hypgeom_legendre_p_ui_one(res, res_prime, n, x, K_one, prec);
    }
    else   /* guaranteed to have |x| < 1 */
    {
        cost_zero = 1e100;
        cost_one = 1e100;
        cost_asymp = 1e100;

        xx = FLINT_MAX(xx, 1e-50);
        xxsub1 = arf_get_d(arb_midref(xsub1), ARF_RND_UP);

        /* Estimate cancellation for series expansion at 0. */
        /* |P_n(xi)| ~= (x+sqrt(1+x^2))^n. */
        cancellation_zero = n * log(xx + sqrt(1.0 + xx * xx)) * 1.44269504088896;
        cancellation_zero = FLINT_MIN(cancellation_zero, 1.272 * n);
        cancellation_zero = FLINT_MAX(cancellation_zero, 0.0);

        /* Estimate cancellation for series expansion at 1. */
        /* For x >= 1, P_n(x) ~= I_0(n*sqrt(2(x-1))) ~= exp(n*sqrt(2(x-1))) */
        if (xxsub1 >= 0.0)
        {
            cancellation_one = 0.0;
        }
        else
        {
            cancellation_one = n * sqrt(2.0*fabs(xxsub1)) * 1.44269504088896;
            cancellation_one = FLINT_MIN(cancellation_one, 2.0 * n);
            cancellation_one = FLINT_MAX(cancellation_one, 0.0);
        }

        d = n / 2;
        K_zero = d + 1;
        K_one = n + 1;
        K_asymp = 1;
        asymp_error = 0.0;

        wp = 1.01 * prec + FLINT_BIT_COUNT(n);
        tolerance = -wp;

        /* adjust for relative tolerance near 0 */
        if (n % 2)
        {
            tolerance += arf_abs_bound_lt_2exp_si(arb_midref(x));
        }

        if (n > 10)
        {
            /* look for early truncation of series at 1 */
            log2u = log(fabs(xxsub1) * 0.5) * 1.44269504088896;
            if (log2u < -30)
                log2u = arf_abs_bound_lt_2exp_si(arb_midref(xsub1)) - 1.0;

            log2x = log(fabs(xx)) * 1.44269504088896;
            if (log2x < -30)
                log2x = arf_abs_bound_lt_2exp_si(arb_midref(x));

            if (n * (2.0 + log2u) < tolerance)
            {
                for (k = 1; k < K_one; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, k)
                        + log2_bin_uiui_fast(n + k, k) + k * log2u;

                    if (size < tolerance)
                    {
                        K_one = k;
                        break;
                    }
                }
            }

            /* look for early truncation of series at 0 */
            if (n * (1.0 + log2x) < tolerance)
            {
                for (k = 1; k < K_zero; k = FLINT_MAX(k+1, k*1.05))
                {
                    size = log2_bin_uiui_fast(n, d - k)
                        + log2_bin_uiui_fast(n+1+2*k, n) - n + 2.0 * k * log2x;

                    if (size < tolerance)
                    {
                        K_zero = k;
                        break;
                    }
                }
            }

            /* look for possible convergence of asymptotic series */
            /* requires information about y = sqrt(1-x^2) */
            yy = arf_get_d(arb_midref(x2sub1), ARF_RND_DOWN);
            yy = sqrt(FLINT_MAX(0.0, yy));
            log2nsy = log(2.0 * n * yy) * 1.44269504088896;

            cost_zero = (prec + cancellation_zero) * K_zero;
            cost_one = (prec + cancellation_one) * K_one;

            for (k = 1; k < n &&
                    prec * k < FLINT_MIN(cost_zero, cost_one);
                        k = FLINT_MAX(k + 1, k * 1.05))
            {
                /* todo: better account for prefactor in the asymptotic series? */
                log2k = log(k) * 1.44269504088896;
                size = 3.0 + k * (log2k - 1.43);  /* estimate K! */
                size -= k * log2nsy;              /* 1/(2n sin(theta))^K */

                if (size < asymp_error)
                {
                    asymp_error = size;
                    K_asymp = k;
                }

                if (size < tolerance)
                {
                    break;
                }
            }
        }

        cost_zero = (prec + cancellation_zero) * K_zero;
        cost_one = (prec + cancellation_one) * K_one;
        cost_asymp = (prec + 0.0) * K_asymp * 2.0;

#if 0
        printf("zero:  K = %ld, cost = %g, cancel %g\n",                K_zero, cost_zero, cancellation_zero);
        printf("one:   K = %ld, cost = %g, cancel %g\n",                K_one, cost_one, cancellation_one);
        printf("asymp: K = %ld, cost = %g, error = %f (tol = %f)\n", K_asymp, cost_asymp, asymp_error, tolerance);
#endif

        if (asymp_error < tolerance && cost_asymp < FLINT_MIN(cost_zero, cost_one))
        {
            arb_hypgeom_legendre_p_ui_asymp(res, res_prime, n, x, K_asymp, wp);
        }
        else if (FLINT_MIN(cost_zero, cost_one) < (1e6 * prec) * prec && n < UWORD_MAX / 4)
        {
            mag_t err1, err2, xrad;
            arb_t xmid;

            mag_init(err1);
            mag_init(err2);
            mag_init(xrad);
            arb_init(xmid);

            arf_set(arb_midref(xmid), arb_midref(x));
            mag_zero(arb_radref(xmid));
            mag_set(xrad, arb_radref(x));

            arb_hypgeom_legendre_p_ui_deriv_bound(err1, err2, n, x, x2sub1);

            if (cost_zero < cost_one)
                arb_hypgeom_legendre_p_ui_zero(res, res_prime, n, xmid, K_zero, wp + cancellation_zero);
            else
                arb_hypgeom_legendre_p_ui_one(res, res_prime, n, xmid, K_one, wp + cancellation_one);

            if (res != NULL)
            {
                mag_mul(err1, err1, xrad);
                arb_add_error_mag(res, err1);
                arb_set_round(res, res, prec);
            }

            if (res_prime != NULL)
            {
                mag_mul(err2, err2, xrad);
                arb_add_error_mag(res_prime, err2);
                arb_set_round(res_prime, res_prime, prec);
            }

            mag_clear(err1);
            mag_clear(err2);
            mag_clear(xrad);
            arb_clear(xmid);
        }
        else if (asymp_error < -2.0)
        {
            /* todo -- clamp to [-1,1]? */
            arb_hypgeom_legendre_p_ui_asymp(res, res_prime, n, x, K_asymp, wp);
        }
        else
        {
            if (res != NULL)
            {
                arf_zero(arb_midref(res));
                mag_one(arb_radref(res));
            }

            if (res_prime != NULL)
            {
                arf_zero(arb_midref(res_prime));
                mag_set_ui(arb_radref(res_prime), n);
                mag_add_ui(arb_radref(res_prime), arb_radref(res_prime), 1);
                mag_mul_ui(arb_radref(res_prime), arb_radref(res_prime), n);
                mag_mul_2exp_si(arb_radref(res_prime), arb_radref(res_prime), -1);
            }
        }
    }

    arb_clear(xsub1);
    arb_clear(x2sub1);
}

#define TIMEIT_PRINT1(__var, __timer, __reps) \
    __var = __timer->cpu*0.001/__reps;

#define TIMEIT_REPEAT1(__timer, __reps) \
    do \
    { \
        slong __timeit_k; \
        __reps = 1; \
        while (1) \
        { \
            timeit_start(__timer); \
            for (__timeit_k = 0; __timeit_k < __reps; __timeit_k++) \
            {

#define TIMEIT_END_REPEAT1(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 30) \
                break; \
            __reps *= 10; \
        } \
    } while (0);

#define TIMEIT_START1 \
    do { \
        timeit_t __timer; slong __reps; \
        TIMEIT_REPEAT1(__timer, __reps)

#define TIMEIT_STOP1(__var) \
        TIMEIT_END_REPEAT1(__timer, __reps) \
        TIMEIT_PRINT1(__var, __timer, __reps) \
    } while (0);

int main()
{
    slong n, d, prec, wp, rep;
    double t1, t2, t3, t4, t5;
    slong cas;

    arb_t x, y, z, nn, mm;
    arb_init(x);
    arb_init(y);
    arb_init(z);
    arb_init(nn);
    arb_init(mm);

    for (cas = 3; cas < 4; cas++)
    {
    for (n = 10; n <= 10000 * 1.2; n = FLINT_MAX(n+1, n*1.2))
    {
        fmpq_poly_t pol;
        fmpq_poly_init(pol);

        if (n % 2 == 1)
            n++;

        prec = n * (1 / 16.);
        prec = 53;

        if (cas == 0)
            prec = 64;
        else if (cas == 1)
            prec = n / 10;
        else if (cas == 2)
            prec = n;
        else
            prec = 10 * n;

        prec = FLINT_MAX(prec, 2);

        fmpq_poly_legendre_p(pol, n);

        arb_poly_t f;
        arb_poly_init(f);

        slong k;
        for (k = 0; k <= n / 2; k++)
        {
            arb_set_fmpz(x, pol->coeffs + 2 * k);
            arb_poly_set_coeff_arb(f, k, x);
        }

        arb_ptr v, w;

        v = _arb_vec_init(n / 2);
        w = _arb_vec_init(n / 2);

        for (k = 0; k < n / 2; k++)
        {
            arb_set_ui(v + k, 4 * (k + 1) - 1);
            arb_div_ui(v + k, v + k, 4 * n + 2, prec);
            arb_cos_pi(v + k, v + k, prec);
            mag_zero(arb_radref(v + k));
        }

        t2 = 1e300;
        for (rep = 0; rep < 10; rep++)
        {
            TIMEIT_START1
            for (k = 0; k < n / 2; k++)
                arb_hypgeom_legendre_p_ui(z, y, n, v + k, prec);
            TIMEIT_STOP1(t3)
            t2 = FLINT_MIN(t2, t3);
            if (t2 > 1.0) break;
        }

        t4 = 1e300;
        for (rep = 0; rep < 10; rep++)
        {
            TIMEIT_START1
            for (k = 0; k < n / 2; k++)
                arb_hypgeom_legendre_p_ui_rec(z, y, n, v + k, prec);
            TIMEIT_STOP1(t3)
            t4 = FLINT_MIN(t4, t3);
            if (t4 > 1.0) break;
        }

        t5 = 1e300;
        for (rep = 0; rep < 10; rep++)
        {
            TIMEIT_START1
            for (k = 0; k < n / 2; k++)
                arb_hypgeom_legendre_p_ui_norec(z, y, n, v + k, prec);
            TIMEIT_STOP1(t3)
            t5 = FLINT_MIN(t5, t3);
            if (t5 > 1.0) break;
        }

        for (k = 0; k < n / 2; k++)
        {
            arb_mul(v + k, v + k, v + k, prec);
            mag_zero(arb_radref(v + k));
        }

        wp = prec + 1.27 * n + 1.6 * n;
        wp = prec + 2.9 * n;

        t1 = 1e300;
if (n < 15000)
{
        for (rep = 0; rep < 10; rep++)
        {
            TIMEIT_START1
//            arb_poly_evaluate_vec_fast(w, f, v, n / 2, wp);
// will total recall 

            arb_ptr * tree;
            tree = _arb_poly_tree_alloc(n);
            _arb_poly_tree_build(tree, v, n / 2, wp);
            _arb_poly_evaluate_vec_fast_precomp(w, f->coeffs, f->length, tree, n / 2, wp);
            _arb_poly_evaluate_vec_fast_precomp(w, f->coeffs, f->length, tree, n / 2, wp);
            _arb_poly_tree_free(tree, n);

            TIMEIT_STOP1(t3)
            t1 = FLINT_MIN(t1, t3);
            if (t1 > 1.0) break;
        }
}

        slong acc = WORD_MAX;
        for (k = 0; k < n / 2; k++)
            acc = FLINT_MIN(acc, arb_rel_accuracy_bits(w + k));

        printf("%ld %ld   %g     %g    %g    %g\n", n, prec, t4, t1, t2, t5);

        fmpq_poly_clear(pol);
        arb_poly_clear(f);
        _arb_vec_clear(v, n / 2);
        _arb_vec_clear(w, n / 2);
    }
    }

    return 0;
}

