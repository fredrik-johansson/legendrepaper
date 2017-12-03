#include "acb.h"
#include "acb_hypgeom.h"
#include "flint/profiler.h"

void
integrand(acb_t res, const acb_t x, slong prec)
{
/*
    acb_add_ui(res, x, 2, prec);
    acb_log(res, res, prec);
*/

    acb_mul_ui(res, x, 10, prec);
    acb_hypgeom_airy(res, NULL, NULL, NULL, res, prec);

/*
    acb_mul_onei(res, x);
    acb_add_ui(res, res, 1, prec);
    acb_gamma(res, res, prec);
*/
}

void
dequad(acb_t res, slong hexp, int doit, slong prec)
{
    arb_t m, r, h, eh, ekh, skh, ckh, emskh, th, ch, w;
    acb_t x, y, sum;
    slong k, eval;
    arb_t tol;

    arb_init(m); arb_init(r); arb_init(h);
    arb_init(eh); arb_init(ekh); arb_init(skh); arb_init(ckh);
    arb_init(emskh); arb_init(th); arb_init(ch); arb_init(w);
    acb_init(x); acb_init(y); acb_init(sum);
    arb_init(tol);

    arb_zero(m);
    arb_one(r);

    arb_one(h);
    arb_mul_2exp_si(h, h, -hexp);

    arb_exp(eh, h, prec);
    arb_one(ekh);

    arb_one(tol);
    arb_mul_2exp_si(tol, tol, -prec);

    eval = 0;

    for (k = 0; ; k++)
    {
/*
        skh = (ekh - 1/ekh)/2
        ckh = (ekh + 1/ekh)/2
        emskh = exp(-skh)
        th = (1 - emskh**2) / (1 + emskh**2)
        ch = 2 * emskh / (1 + emskh**2)
*/

        /* skh, ckh = sinh(kh), cosh(kh) */
        arb_inv(skh, ekh, prec);
        arb_add(ckh, ekh, skh, prec);
        arb_sub(skh, ekh, skh, prec);
        arb_mul_2exp_si(skh, skh, -1);
        arb_mul_2exp_si(ckh, ckh, -1);
        /* emskh = exp(-sinh(kh)) */
        arb_neg(emskh, skh);
        arb_exp(emskh, emskh, prec);

        /* th = tanh(sinh(kh)) = (1 - emskh^2) / (1 + emskh^2) */
        /* ch = 1/cosh(sinh(kh)) = 2 emskh / (1 + emskh^2) */
        arb_mul(ch, emskh, emskh, prec);
        arb_sub_ui(th, ch, 1, prec);
        arb_neg(th, th);
        arb_add_ui(ch, ch, 1, prec);
        arb_div(th, th, ch, prec);
        arb_div(ch, emskh, ch, prec);
        arb_mul_2exp_si(ch, ch, 1);

        /* x = m + r * th */
        arb_mul(acb_realref(x), r, th, prec);
        arb_add(acb_realref(x), acb_realref(x), m, prec);
        arb_zero(acb_imagref(x));

        /* w = r * h * ckh * ch^2 */
        arb_mul(w, ch, ch, prec);
        arb_mul(w, w, ckh, prec);
        arb_mul(w, w, h, prec);
        arb_mul(w, w, r, prec);

        if (arb_lt(w, tol))
        {
            break;
        }

        if (doit)
        {

            /* y = f(x) */
            integrand(y, x, prec);
            eval++;

            if (k != 0)
            {
                /* x = m - r * th */
                arb_mul(acb_realref(x), r, th, prec);
                arb_sub(acb_realref(x), m, acb_realref(x), prec);
                arb_zero(acb_imagref(x));

                integrand(x, x, prec);
                eval++;

                acb_add(y, y, x, prec);
            }

            acb_addmul_arb(sum, y, w, prec);
        }

        arb_mul(ekh, ekh, eh, prec);
    }

    printf("eval %ld\n", eval);

    acb_set(res, sum);

    arb_clear(m); arb_clear(r); arb_clear(h);
    arb_clear(eh); arb_clear(ekh); arb_clear(skh); arb_clear(ckh);
    arb_clear(emskh); arb_clear(th); arb_clear(ch); arb_clear(w);
    acb_clear(x); acb_clear(y); acb_clear(sum);
    arb_clear(tol);
}

int main()
{
    acb_t res;
    arb_t ans;
    slong prec, n;

    acb_init(res);
    arb_init(ans);

    for (n = 7; n <= 9; n++)
    {

    prec = 3408;
    TIMEIT_START
    dequad(res, n, 0, prec);
    TIMEIT_STOP
    TIMEIT_START
    dequad(res, n, 1, prec);
    TIMEIT_STOP

/*
    arb_log_ui(ans, 3, prec);
    arb_mul_ui(ans, ans, 3, prec);
    arb_sub_ui(ans, ans, 2, prec);
*/

        arb_set_str(ans, "0.10990317364333819333679412354778428474599847954105073385203878829379756653063638802899923350903098290084544000716092300224638014883498907909937547546057918360586542802015929771412842672793896412651154877388183163187749381201110304596029598068840079719994820587716417997769733215720194885772488875417932147835404272225778303818962923684370056841150538028772951375462244933214536295145531084365176443776486676351535997642529430579371872216283291598287525272042617476846433468682148951903566637658293851754737801011119019162275127454432182631765432267748173425851255992132423475567929072633263742734845842276918217159364133165355168570493428995229209085852419126876888871376665686942983220517100227411698708452503609940208631578244381606423328810877071411807925895716108143834636019373295465493569894268921520321195392834866522818779193384033150976179363433593935378938110641049936027823447608014746852226229433683825915803439741613598966991065510411914080495661160745803784953047805643683571543228900105332339572465588026252979000187988555465644352567604295538355301635939443229838264248794469171429981231388851167778351752122785408885847085519220934635677105282815636327923794647435105886331453250882348880611176955516100177606258423516911454806476403256001290371362487214192705664950702887631953984962346391844450983915682320636646379664651210824968646372509052689791599418031105995623107445261096849420676909604518164738242363482685184669375156867170110478007201914411897820640561291858825709444063", prec);

/*
        arb_set_str(ans, "1.57239266949806559075198297278478209135231172387993819069313546393443758953396973642297613604098814617137347054301950568564423437784143953198884437745463958850066838746298693743121116367225999188039745778437224356721527775098778780058196835359347901900868493176978966522723290073721765442260140089887879310539559950005889031808529424106615471409357523859689924174656147994835184657673238252910791012308069100786753213003217184671360430728995240063246077713231739446079241460198006838693324424960823832455701711152398258931166025007655286284356510603191038504006607212832774762593165004467470061268686736145526122777106666140133340996660135637161281616253060971586067494983744673906624698662874339216422262555653389520600023731245476749363887592211814621370371108413066160817506304770001995416228554114898162596946474284126931797360268771390029347960491211725589036254309986007171293083829916490396339046076547249721182317246635532646722118439474298919066925642617797525030972148714872870635023412791856194084485749502528315843639197700686261344676540411537163287340007610732800972884902590813135765869419114842682402323743664713107217253840523193880617104648177800811523870470975724255071784420478657063", prec);
*/

    acb_printn(res, 30, 0); printf("\n");
    arb_printn(ans, 30, 0); printf("\n");

    arb_sub(ans, ans, acb_realref(res), prec);
    arb_printn(ans, 30, 0); printf("\n");

    }
}

