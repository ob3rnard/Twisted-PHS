from sage.all import *
import fp
from lattice import *
from number_field import *

# ------------------------------------------------------------------------------------------------
# This implements Tw/Opt/()-PHS strategies
# Generically, there is a method flag:
__METHOD_TYPE      = ['PHS', 'OPT', 'TW', 'NONE'];
# If 'NONE' is specified, there are overrides by specific options, but it is not guaranteed to work.
# Functions are named: "twphs_..." for addressing the three methods, and calls specifics:
#    phs_... / opt_... / tw_...
# whenever needed depending on __METHOD_TYPE passed to the function.
# 'PHS' is (almost) iso-functional to [PHS19] code.

# There are several steps
# 1. twphs_get_fb:         choose Factor base.
# 2. twphs_get_matrix:     lattice construction.
# 3. twphs_get_target:     compute the target from the cldl according to chosen method & above lattice.
# 4. twphs_build_solution: build solution to idSVP



# -------------------------------------------------------------------------------------------------
# 1. FACTOR BASE CHOICE
#    Dimensioning depends on:
#          ClK (FB must generate),
#          Disc K (phs),
#          Scaling 'c' (PHS/OPT),
#          Targeted root volume (OPT): there is an *unused* __LNROOTVOL_TARGET set to 0.3.
#          RAND (PHS)
#          Class_number hK and regulator RK (TW). Note: for cyclotomics, HR can be found using zeta_K.
#              (for NTRU Prime, it seems more complicated, which is weird)

# Roughly averages the log reduced volume of [PHS19] on small examples.
__LNROOTVOL_TARGET = 0.3;

# Scaling
# --------------------------------------------------
# [PHS19]: c is chosen in  as n^1.5/#fb
def phs_get_c(n, n_fb, b_prec=fp.BIT_PREC_DEFAULT):
    return (RealField(b_prec)(n)**RealField(b_prec)(3/2))/RealField(b_prec)(n_fb);

# Opt-PHS scaling
# But isometry induces a gap of 1+ln(n), not sqrt(n)
def opt_get_c(n, fb, b_prec=fp.BIT_PREC_DEFAULT):
    _n = RealField(b_prec)(n);
    _s_fb = sum(ln(RealField(b_prec)(_p.norm())) for _p in fb);
    return RealField(b_prec)(max(1, (1+ln(_n))*_n/_s_fb)); # Almost always 1.

# Dispatch function Tw-/Opt-/- PHS for scaling.
def twphs_get_c(n, fb, method, b_prec=fp.BIT_PREC_DEFAULT):
    assert (method in __METHOD_TYPE);

    if (method == 'PHS'): return phs_get_c(n, len(fb), b_prec=b_prec);
    if (method == 'OPT'): return opt_get_c(n, fb, b_prec=b_prec);
    if (method == 'TW') or (method == 'NONE'): return RealField(b_prec)(1);


# Dimension FB
# -------------------------------------------------
# PHS dim_FB
#   As in [PHS19] code
#   Rk: According to the article, it should be sthg like ln(D)+n.lnln(D)-rkU (almost twice as large)
def phs_get_dimfb(K):
    _rk_u = get_rank_units(K);
    _DK   = K.discriminant().abs();
    _r    = RDF(ln(_DK)).floor() - _rk_u;
    return _r;

# Opt-PHS dim_FB
#   This version targets same root volume as PHS (Sorry for the implementation quality)
#   Could also aim for root volume of exp(0.3)
def opt_get_dimfb(K, fb_cl, fb_oth):
    n      = K.degree();
    r1, r2 = K.signature();
    kphs   = phs_get_dimfb(K);
    cphs   = phs_get_c(n, kphs);
    hR     = K.class_number()*K.regulator();
    un_rk  = get_rank_units(K);
    lnVL   = RDF(ln(hR)+ln(2)*(-r2/2.)+ln(sqrt(n)));
    lnVphs = RDF(1/(un_rk+kphs)*(lnVL+un_rk*ln(cphs))); # or, I guess, exp(__LNROOT_TARGET)
    cst_target = RDF(lnVL - un_rk*lnVphs);

    # /!\ Warn: fb_cl can be empty (principal fields like generic NTRU)
    #     As implemented here, the starting point cannot be k=0.
    # >>> We return tha max (k, #fb_cl) instead.
    fb_all = fb_cl + fb_oth;
    k=1; S=RDF(ln(fb_all[0].norm()));
    while (k*lnVphs - un_rk*ln( max(1.,n*(1+ln(n))/S) )) < cst_target:
        k=k+1;
        S=S+RDF(ln(fb_all[k-1].norm()));
    
    return max(k, len(fb_cl));


# Tw-PHS dim_FB
#   Add ideals until the reduced volume rises.
#   Rk: we don't cut in the middle of iso-norm ideals.
def tw_get_dimfb(K, fb_cl, fb_oth):
    n      = K.degree();
    r1, r2 = K.signature();
    un_rk  = get_rank_units(K);
    hR     = K.class_number()*K.regulator();

    fb_all = fb_cl + fb_oth;
    k      = len(fb_cl); # Start with at least the generators of ClK
    while (fb_all[k].norm().log().log() < 1/(un_rk+k)*(ln(sqrt(n+k))-ln(2)*r2/2.+ln(hR)+sum(log(log(pp.norm())) for pp in fb_all[:k])) - ln(sqrt(1+1/(n+k)))):
        k = k+1;

    # Include all isonorm ideals of the last one.
    pmax = fb_all[k-1].norm();
    while (fb_all[k].norm() == pmax):
        k=k+1;
    
    return k;


# Dispatch functions
def twphs_get_dimfb(K, method, fb_cl, fb_oth):
    assert (method in __METHOD_TYPE and method != 'FORCED');

    if (method == 'PHS'): return phs_get_dimfb(K);
    if (method == 'OPT'): return opt_get_dimfb(K, fb_cl, fb_oth);
    if (method == 'TW'):  return tw_get_dimfb(K, fb_cl, fb_oth);


# Compute Factor Base
# --------------------------------------------------
#  K: number field
#  If method == 'NONE', then we look at the following:
#  [opt] r = size:     size of the factor base (NB: only the generating set is extracted from ClK_FB)
#                      By default (0), compute r with pmhs18_get_r()
#                      Note that r must be >= number of generators of cl_K
#  [opt] rand = False: take ideals from smallest norm to highest norm in the FB until r is reached,
#            (= True : take random ideals in the set of the first 3r prime ideals)
#
#  For any method:
#  a. Compute ClK (if not cached in Sage)
#  b. Compute the first 3r ideals of smallest norms, for r = floor(ln(disc_K)) - rk_U
#  c. Compute from this the final dimension (for OPT / TW)
#  d. Extract the first ones (OPT, TW) or sample at random in this (PHS or rand option)
import random;
def twphs_get_fb(K, method='TW', r=0, rand=False):
    assert(method in __METHOD_TYPE);

    # ------- List sufficiently many ideals (3*dim)
    # First we constitute the panel of first 3*ln(DK) ideals
    _r    = r if (method == 'NONE') else phs_get_dimfb(K); # Meth. 'OPT' and 'TW' are upper bdd by 'PHS'

    # ------- Ensure we generate ClK 
    # So, the first ideals in the factor base are generators for cl_K. Will be feasible if r >= log hK
    # This is not really necessary, but let's assume that we insist on being able to generate all cl_K
    t = cputime();
    _cl_K  = K.class_group(proof=False); # Normally, this computation is cached if already done for K
    t = cputime(t); print("Cl_K: {:.2f}s".format(t));
    _fb_cl = list(_cl_K.gens_ideals()); assert(len(_fb_cl) <= _r); # NB: Necessarily *not* principal
    
    # Elaborate a set of at least 3*r ideals (including ClK generators)
    # --------------------------------------
    # At the end the set contains 3r+whatever is needed to have all ideals of the last norm
    # Use a silly heuristic to avoid multiple tries: for integers, pk \leq k(\ln k + \lln k)
    # Adjust the output by 1.5 to compensate for the "Prime Ideal" setting
    _B       = ceil(1.5*(3*_r)*(ln(3*_r)+ln(ln(3*_r))));
    # Exclude also ramified primes.
    _bdd_ids = [ pid for pid in K.primes_of_bounded_norm(_B) if (pid.ramification_index() == 1) ];
    while (len(_bdd_ids) < 3*_r): # In practice, never triggered
        print("fb: Ooops... Bound {} insufficient, increasing. {}".format(_B,len(_bdd_ids)));
        _B       = ceil(1.2*_B);
        _bdd_ids = [ pid for pid in K.primes_of_bounded_norm(_B) if (pid.ramification_index() == 1) ];

    # Restrict to 3r + eps (don't truncate btw ideals of same norm)
    _bdd_ids = _bdd_ids[:3*_r]+[_id for _id in _bdd_ids[3*_r:] if (_id.norm()==_bdd_ids[3*_r-1].norm())];
    # Now remove ideals already in _fb_cl (note id in _fb_cl =\> id in _oth_ids, though likely)
    _oth_ids = [ _id for _id in _bdd_ids if not (_id in _fb_cl)];
    
    # Obtain final FB dimensions
    # --------------------------------------
    _k = r if (method == 'NONE') else twphs_get_dimfb(K, method, _fb_cl, _oth_ids);

    # Enlarging the factor base.
    # --------------------------
    _n_adds = _k - len(_fb_cl); # assert above: this is >= 0
    if (method == 'PHS') or (method == 'NONE' and rand == True): # rand ignored if method != 'FORCED'
        _fb = _fb_cl + random.sample(_oth_ids, _n_adds); # Ensures all distinct
    else:
        _fb = _fb_cl + _oth_ids[:_n_adds];

    return _fb;


# Original interface
def phs_get_fb(K):
    return twphs_get_fb(K, method='PHS');

def opt_get_fb(K):
    return twphs_get_fb(K, method='OPT');

def tw_get_fb(K):
    return twphs_get_fb(K, method='TW');



# -------------------------------------------------------------------------------------------------
# 2. GET LOG S-UNIT LATTICE BASIS
#    Depends on:
#        - Scaling (c)
#        - Chosen logarithmic embeddings (see number_field.py for option details)
#        - Isometry/Projection choices: get_twfHcE_matrix()

# Projection/isometry thing.
# --------------------------------------------------
# a. Always project on some H: first n coordinates for PHS/OPT, all n+#fb for TW
#    Note that currently, method 'NONE' is projecting like PHS/OPT.
# b. Scaling is done *before* the isometry to avoid branching on dimensions if 'NONE'.
#    It commutes anyway, modulo numeric stability.
# c. Full-rank strategy:
#        PHS:  Aggressive pruning of the r2 (duplicated) coordinates, + the last one (arbitrary choice)
#        OPT:  Isometry sending first n coordinates to the space orthogonal to 1, and duplicates to 0.
#        TW:   Isometry on the whole space, sending 1,...,1 and duplicates to 0
#        NONE: Do nothing.
def get_twfHcE_matrix(r1, r2, fb, method, b_prec=fp.BIT_PREC_DEFAULT):
    assert (method in __METHOD_TYPE);

    _n  = r1+2*r2;
    _nu = r1+r2-1;
    _k  = len(fb);
    
    # Always project on H or H0
    _span_H = (_n+_k) if (method == 'TW') else (_n);
    _pH     = get_projection_H(_n+_k, _span_H, b_prec=b_prec);
    assert (fp.fp_check_zero("PrH(1)", (vector([1]*(_span_H)+[0]*(_n+_k-_span_H))*_pH).coefficients(), target=b_prec));
    
    # Scaling by c (done before because of the 'NONE' iso type
    _c      = twphs_get_c(_n, fb, method, b_prec=b_prec);
    _c_Id   = block_diagonal_matrix(_c*identity_matrix(_n), identity_matrix(len(fb)));
    assert (_c_Id.base_ring().precision() >= b_prec);
    
    # Full-rank strategy
    if (method == 'NONE'): # Do nothing.
        _fH = matrix(RealField(b_prec), identity_matrix(_n+k));
    elif (method == 'PHS'): # Prune r2 coordinates + last remaining
        _fH = block_diagonal_matrix([get_ext_pruning(r1,r2), identity_matrix(_k)],
                                    subdivide=False);
        _fH = _fH.change_ring(RealField(b_prec));
    elif (method == 'OPT'): # MK-isometry on _n first
        _fH = block_diagonal_matrix([get_minkH(r1,r2,0,b_prec=b_prec), identity_matrix(_k)],
                                    subdivide=False);
        assert (fp.fp_check_zero("Vol(fH)==1", [vol(_fH.transpose())-RealField(b_prec)(1)],
                                 target=b_prec, sloppy=True)); # /!\ Sloppy
    elif (method == 'TW'): # MK-isometry on _n+_k (whole space)
        _fH = get_minkH(r1, r2, _k, b_prec=b_prec);
        assert (fp.fp_check_zero("Vol(fH)==1", [vol(_fH.transpose())-RealField(b_prec)(1)],
                                 target=b_prec, sloppy=True)); # /!\ Sloppy
        
    # Finally, proj -> scaling -> full-rank
    _fHcE = (_pH*_c_Id*_fH);
    assert (_fH.base_ring().precision() >= b_prec);
    assert ((_fHcE.nrows() == _n+_k) and (_fHcE.ncols() == _nu+_k));
    assert (_fHcE.base_ring().precision() >= b_prec);

    return _fHcE;


# Provided for one shots, but inefficient as soon as several images are needed.
def twfHcE(v, r1, r2, fb, method, b_prec=fp.BIT_PREC_DEFAULT):
    assert(v.base_ring().precision() >= b_prec);
    _fHcE   = get_twfHcE_matrix(r1, r2, fb, method, b_prec=b_prec);
    return v*_fHcE;


# Get Log S-unit lattice basis 
# --------------------------------------------------
# un:     units
# s_un:   S-units
# p_inf:  Set of infinite places.
#             The aim is to have a consistent order between several calls or Sage sessions
#             [It seems consistent always, but there is no guarantee in the documentation]
#             Note that p_inf[].codomain().precision() should be sufficient to handle sun coefficients.
#             [This is handled internally, but it saves some time.]
# fb:     List of finite places (factor base of prime ideals)
# b_prec: Output bit precision (of course, it should not exceed the precision of p_inf).
def twphs_get_matrix(un, s_un, p_inf, fb, method, b_prec=fp.BIT_PREC_DEFAULT):
    assert ((len(fb) == len(s_un))
            and (len(un) == get_rank_units(fb[0].number_field()) and len(p_inf) == len(un)+1)
            and (p_inf[0].codomain().precision() >= b_prec));
    _K       = fb[0].number_field();
    _r1, _r2 = _K.signature();
    _n_inf   = get_nb_inf_places(_K);
    
    # Transformation to full rank
    _fHcE  = get_twfHcE_matrix(_r1, _r2, fb, method, b_prec=b_prec);

    # Work precision for the log(su) part
    _w_prec = max([max([RealField(1000)(log(_coef.abs())) for _coef in _su.list()]) for _su in un + s_un]); # All coefficients are integers, if s_un is S-units(fb).
    p_inf   = extend_inf_places(_K, p_inf, to_prec=_w_prec);
    
    # Matrix with expected logs + val (one row for each [S-]unit)
    _inf_type = 'EXPANDED'; # s1,..,sr1,sr1+1,sr1+1,...,sr1+r2,sr1+r2.
    _fb_type  = 'TWISTED' if (method == 'TW') else 'FLAT'; # flat is just -vp(elt).
    _Log_sun  = matrix([ log_embedding(_su, p_inf, fb=fb,
                                       inf_type=_inf_type, fb_type=_fb_type, b_prec=b_prec)
                         for _su in un + s_un ]);
    
    # Apply pH.fH.c
    _twphs = _Log_sun*_fHcE; # Note: probably worth to apply fHcE at each step of the _Log_sun line.
    assert (_twphs.base_ring().precision() >= b_prec);
    return _twphs;


# Original interface
def phs_get_matrix(un, s_un, p_inf, fb, b_prec=fp.BIT_PREC_DEFAULT):
    return twphs_get_matrix(un, s_un, p_inf, fb, b_prec=b_prec, method='PHS');

def opt_get_matrix(un, s_un, p_inf, fb, b_prec=fp.BIT_PREC_DEFAULT):
    return twphs_get_matrix(un, s_un, p_inf, fb, b_prec=b_prec, method='OPT');

def tw_get_matrix(un, s_un, p_inf, fb, b_prec=fp.BIT_PREC_DEFAULT):
    return twphs_get_matrix(un, s_un, p_inf, fb, b_prec=b_prec, method='TW');



# -------------------------------------------------------------------------------------------------
# 3. GET TARGET AFTER CLDL

# Guess the value of beta (the constant to apply to each coordinate)
# In the twisted case, the final drift is <cst> - ln N(p).
#     For the moment, we apply <cst> - ln N(pmin) everywhere, and this '-ln N(pmin)' is put in the guess
#     --> This works already well.
def twphs_guess_beta(a, bkzlat, K, fb, method):
    if ((method == 'OPT') or (method == 'PHS')):
        return (vol_reduced(bkzlat)-bkzlat.base_ring()(0.8));
    if (method == 'TW'):
        _n = K.degree();
        _k = len(fb);
        _sum = sum(bkzlat.base_ring()(ln(_pid.norm())) for _pid in fb);
        _vol = vol_reduced(bkzlat);
        return bkzlat.base_ring()( (_k+_n)/_n*_vol + ln(a.norm())/_n - _sum/_n - ln(fb[0].norm()));


# Construct the target for idCVP in the log-S-unit lattice
# - eta, a, fb are such that, there exists some vp for p in fb st <eta> = a.\prod_{p in fb} p^vp
#   So morally, the target is fHcE( [Ks:R] log |s(eta)|, -vp log Np )
# - beta is the drift parameter to ensure that the output solution will have all positive valuations
#   the target will be shifted according to method.
#   Note that the drift (even in TW case) is always <cst> + something that at worst depend on fb[i].
# - b_prec should be the precision of the lattice (BKZ) matrix.
# - method must correspond to the targeted lattice from phs_get_matrix()
#       If multiple targets from same lattice, precomputing fHcE is useful
#       (the output of get_twfHcE_matrix()
#   --> opt: _pcmp_fhce = matrix
def twphs_get_target(eta, a, p_inf, fb, method, beta=0.0, b_prec=fp.BIT_PREC_DEFAULT, _pcmp_fhce=0):
    assert (p_inf[0].codomain().precision() >= b_prec);
    _K        = fb[0].number_field();
    _n        = _K.degree();
    _r1, _r2  = _K.signature();

    # Work prec: enough to handle (rational) coefficients of eta
    #_w_prec = max([max([RealField(1000)(log(_coef.abs())) for _coef in _su.list()]) for _su in un + s_un]); # All coefficients are integers, if s_un is S-units(fb).
    #p_inf   = extend_inf_places(_K, p_inf, to_prec=_w_prec);
    _work_prec = p_inf[0].codomain().precision(); # Above takes too long, assume p_inf HAS enough prec
    
    # Log embedding
    _inf_type = 'EXPANDED'; # s1,..,sr1,sr1+1,sr1+1,...,sr1+r2,sr1+r2.
    _fb_type  = 'TWISTED' if (method == 'TW') else 'FLAT'; # flat is just -vp(elt).
    # Log of eta, and compensation for the p|a for p in FB
    _log_eta = log_embedding(eta, p_inf, fb=fb, inf_type=_inf_type, fb_type=_fb_type, b_prec=_work_prec);
    _log_a   = log_embedding(a,   p_inf, fb=fb, inf_type='NONE',    fb_type=_fb_type, b_prec=_work_prec); # inf_type='NONE' triggers "Ideal" code -> only FB
    _log_t   = _log_eta - _log_a; assert (_log_t.base_ring().precision() >= _work_prec);
    _log_t.change_ring(RealField(b_prec));
    if (_fb_type == 'TWISTED'): # Otherwise, there is nothing we could check
        assert (fp.fp_check_zero("sum(log_t)-N(a)", [sum(_log_t)-a.norm().log(prec=b_prec)], target=b_prec));

    # Construct drift vector: lambda_p = beta - ln N(p).
    # Twisted case:
    #     We removed the ln N(p) for now: this brings more complications for "exact_cvp" (infinite loop)
    #     Maybe we should try "max (0, beta - ln N(p))".
    #     This works already well using a constant "- ln N(pmin)" (see twphs_guess_beta('TW')).
    _beta    = RealField(b_prec)(beta);
    _drift_v = vector([RealField(b_prec)(0)]*_n + [_beta]*len(fb));

    # Projection wizardry
    if (_pcmp_fhce == 0):
        _log_t   = twfHcE(_log_t + _drift_v, _r1, _r2, fb, method, b_prec=b_prec);
    else:
        _log_t   = (_log_t + _drift_v)*_pcmp_fhce;

    assert (_log_t.base_ring().precision() >= b_prec);
    return _log_t;



# -------------------------------------------------------------------------------------------------
# 4. SOLUTION BUILDING IN K
#    Note this last step does *not* depend on the method.

# Construct the candidate solution corresponding to log_..(eta) - log_cvp(BL, eta)
# - BL  is the lattice description, each column corresponding to one S-unit in [u cat su]
#   --> Use the *original* lattice basis, not the BKZ one !
# - cldl is the target st. <cldl> = target_ideal . prod_{p\in fb} p^vp
# - log_s_cvp is some close vector of log_..(cldl) in lattice BL
# - The output is eta/s, where s is st log_..(s) = log_s_cvp.
def twphs_build_solution(cldl, log_s_cvp, BL, u_su):
    _b_prec = log_s_cvp.base_ring().precision(); # log_s_cvp can come from a BKZ basis with less precision than BL

    # Carefully compute _y \in ZZ^dim st _y*BL = log_s_cvp
    _y_R = BL.solve_left(log_s_cvp, check=False); # check=True only works on exact rings
    _y   = vector(map(round, _y_R));
    assert (fp.fp_check_zero("y_R-ZZ", (_y_R-_y).coefficients(), target=_b_prec)); # Should be fine is BL has enough precision.
    
    # S-unit corresponding to log_s_cvp
    assert (len(u_su) == len(_y));
    # This can be (very) long !! Other methods (dividing step by step or ...) seem even longer.
    _s   = prod(map(pow, u_su, _y)); 

    return (cldl/_s);

