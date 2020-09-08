from sage.all import *
import fp
from lattice import *


# Hopefully, always use the (e)GRH hypothesis and useful unproven conjectures for number fields
# We will still try to specify it each time it is needed, for debug and checks
proof.number_field(False);


# -------------------------------------------------------------------------------------------
# Some labelling
#     Cyclotomic fields of conductor m: 'z<m>' (ex.: 'z23') [deg: 22]
#     NTRU Prime fields x^p-x-1:        'n<p>' (ex.: 'n23') [deg: 23]
#     (Default):                        'g<d>' (ex.: 'g23') [deg: 23]
def nf_get_tag(K):
    _Zx  = PolynomialRing(Integers(), names=('x',));
    (x,) = _Zx._first_ngens(1);

    if (K.defining_polynomial().is_cyclotomic() == True): # Cyclotomic
        _tag = "z{}".format(K.conductor());
    elif (K.defining_polynomial() == (x**K.degree()-x-1)): # NTRU Prime
        _tag = "n{}".format(K.degree());
    else: # Generic
        _tag = "g{}".format(K.degree());
    
    return _tag;


def nf_set_tag(tag, eq=[]):
    _typ = tag[0];
    _val = sage_eval(tag[1:]);

    if (_typ == "z"):
        _K      = CyclotomicField(_val);
    elif (_typ == "n"):
        _Zx  = PolynomialRing(Integers(), names=('x',));
        (x,) = _Zx._first_ngens(1);
        _K   = NumberField(x**_val - x - 1, name='zp');
    else:
        _K = NumberField(PolynomialRing(Integers(), 'x')(eq), name='a');
    
    return _K;



# -------------------------------------------------------------------------------------------
# Algebraic numbers things

# There seemed to exist a precision drift somewhere (?).
# At the time, doubling precision on infinite places corrected the phenomenon.
__BIT_PREC_INFP_SCALE = 2.0;


# Returns the number of infinite places of K (r1+r2)
def get_nb_inf_places(K):
    return sum(K.signature());

# Returns infinite places of suitable precision to reach b_prec. [In practice, 2*b_prec]
# Rk: At some point, there was some precision issues with complex embeddings (?).
def get_inf_places(K, b_prec=fp.BIT_PREC_DEFAULT):
    return K.places(prec=__BIT_PREC_INFP_SCALE*b_prec);

# Extends the precision of infinite places
# This is done by recomputing new places, and verifying the order matches.
def extend_inf_places(K, p_inf, to_prec=fp.BIT_PREC_DEFAULT):
    assert(len(p_inf) == get_nb_inf_places(K));
    if (to_prec*__BIT_PREC_INFP_SCALE <= p_inf[0].codomain().precision()):
        return p_inf;
    _new_p_inf = K.places(prec=to_prec*__BIT_PREC_INFP_SCALE);
    assert (fp.fp_check_zero("phi++-phi", [ _new_p_inf[_k](K.gen())-p_inf[_k](K.gen()) for _k in range(len(p_inf)) ], target=p_inf[0].codomain().precision()));
    return _new_p_inf;


# T2-norm
# This is the l2-norm of some complex embedding of alpha.
# Unlike Magma, we don't return the square of T2-norm, though it would be an integer (or rational)
def t2_norm(alpha, b_prec=fp.BIT_PREC_DEFAULT):
    # Beware of precision issues with complex_embeddings. Twice should compensate.
    return vector(ComplexField(b_prec),
                  alpha.complex_embeddings(prec=__BIT_PREC_INFP_SCALE*b_prec)).norm();


# Logarithmic embedding.
# When conjugates are present:
#     the order of infinite embeddings is as follows : 1...r_1,  ..., j, \conj(j), ...
#     the order of finite embeddings follows the same principle: -vp ln p, ..(\times f=[OK/P:Z/p]), -vq etc
# For finite/infinite parts, there are 3 options: TWISTED, EXPANDED and FLAT: (dimensions given if applied on both parts)
#     TWISTED:     [Ks:Qs] ln |a|_s                                      Dim: r_1+r_2  / k
#     FLAT:        ln |a|_s on infinite places, -vp(a) on finite places  Dim: r_1+r_2  / k
#     EXPANDED:    ln |a|_s, [Ks:Qs] times                               Dim: r_1+2r_2 / sum([OK/p:Z/p])
# For ideals, we add the following on infinite places:
#     NONE:        Put 0 everywhere.
# Note the "FLAT" version is designed to fit more or less PHS's choice on finite valuations;
# we would naturally expect instead ln |a|_s on all places (= vp ln(p) for finite places).
# /!\ p_inf[].codomain().precision() shuld be large enough to handle coefficients of elt.
__LOG_TYPE = ['TWISTED', 'FLAT', 'EXPANDED', 'NONE'];
def log_embedding(elt, p_inf, fb=[], inf_type='TWISTED', fb_type='TWISTED', b_prec=fp.BIT_PREC_DEFAULT):
    _K          = elt.number_field() if (inf_type == 'NONE') else elt.parent();
    _p_inf_prec = p_inf[0].codomain().precision();
    assert(len(p_inf) == get_nb_inf_places(_K));
    assert((inf_type in __LOG_TYPE) and (fb_type in __LOG_TYPE));
    assert(_p_inf_prec >= b_prec*__BIT_PREC_INFP_SCALE,
           "Precision of infinite places ({}) not sufficient (expected:{})".format(_p_inf_prec, b_prec*__BIT_PREC_INFP_SCALE));

    _r1, _r2 = _K.signature();
    # Several remarks:
    # 1. abs_val return |s(eta)|^[Ks:R] on infinite places, N(p)^(-v_p(eta)) on finite places
    # 2. the log function has some weird precision issues (?)
    # 3. going through finite places with abs_val() takes too much time in Sage;
    #    computing it from scratch (valuation, norm) can be up to 3 times faster
    if (inf_type == 'FLAT'):
        _log_inf = [ RealField(b_prec)(_phi(elt).abs().log()) for _phi in p_inf ];
    elif (inf_type == 'TWISTED'):
        _log_inf = \
            [ RealField(b_prec)(_phi(elt).abs().log()) for _phi in p_inf[:_r1] ] + \
            [ 2*RealField(b_prec)(_phi(elt).abs().log()) for _phi in p_inf[_r1:] ];
    elif (inf_type == 'EXPANDED'):
        _log_inf = \
            [ RealField(b_prec)(_phi(elt).abs().log()) for _phi in p_inf[:_r1] ] + \
            sum((2*[RealField(b_prec)(_phi(elt).abs().log())] for _phi in p_inf[_r1:]), [] );
    elif (inf_type == 'NONE'):
        _log_inf = [ RealField(b_prec)(0) for _i in range(_K.degree())];
        
    if (fb_type == 'FLAT'):
        _log_fb  = [ -elt.valuation(_pid) for _pid in fb ]; # PHS's version, but logically this would be -vp(eta) ln(p) (vs. -vp(eta) ln(norm(pid)) for TWISTED)
    elif (fb_type == 'TWISTED'):
        _log_fb  = [ -elt.valuation(_pid) * _pid.norm().log(prec=b_prec) for _pid in fb ];
    elif (fb_type == 'EXPANDED'):
        _log_fb  = sum((_pid.residue_class_degree() * [-elt.valuation(_pid)*_pid.smallest_integer().log(prec=b_prec)] for _pid in fb), []);

    _log_embedding = vector(_log_inf + _log_fb);
    assert (_log_embedding.base_ring().precision() >= b_prec);
    return _log_embedding;



# ------------------------------------------------------------------------------------------
# Units: particular units (eg. cyclotomic) and S-Units

# Obtain the rank of the (fundamental) unit group
def get_rank_units(K):
    return sum(K.signature())-1;


# O is an order of some number field, fb is a list of prime ideals of O
def get_S_units(K, fb=[]):
    # In some versions of Sage, S_units() is taking forever,
    # while S_unit_group() then gens_values() works fine. Hence the workaround.
    _sU_gp  = K.S_unit_group(S=fb, proof=False); 
    _sU     = _sU_gp.gens_values();

    # Since they are computed anyway, the second return value are fundamental units.
    _s_un   = [ _su for _su in _sU if (_su.norm().abs() != 1) ];
    _un     = [ _su for _su in _sU if (_su.norm().abs() == 1) and (_su.multiplicative_order() == +Infinity) ];
    assert ((len(_un) == get_rank_units(K))
            and (len(_s_un) == len(fb)));
    
    return _sU_gp, (_un, _s_un);



# ------------------------------------------------------------------------------------------
# Solving ideal-SVP in NumberFields

# Input: an ideal of some number field
# Output: a shortest element of ideal
# Shortest means the output of fplll-SVP shortest_vector function on the Minkowski embedding of ideal
def idsvp_exact(ideal, b_prec=fp.BIT_PREC_DEFAULT, approx_bkz=False):
    _K       = ideal.number_field();
    _r1, _r2 = _K.signature();

    # Embedding lattice in Minkowski space
    # minkowski(x1,x2,..) returns M = ( S(x1) | S(x2 | ... ) in columns, hence the transposition.
    _id_ZB   = ideal.basis();
    _id_L    = Matrix(RealField(b_prec), _K.minkowski_embedding(_id_ZB, prec=2*b_prec).transpose());
    if (approx_bkz == False):
        _s_L     = svp_exact(_id_L);
    else:
        _id_bkz, _ = bkz(_id_L, block_size=40, bkzmaxloops=300);
        _s_L       = _id_bkz[0];
    
    # Going back... is just inverting the Minkowski matrix (on the Z-basis of ideal)
    _s_ZB    = _s_L*rc_mat_inverse(_id_L);
    #_s_ZB      = _id_L.solve_left(_s_L, check=False); # (check=True) work only for exact rings
    _dist_s_ZZ = [ _coef - round(_coef) for _coef in _s_ZB ];
    assert (fp.fp_check_zero("_s_ZB-ZZ", _dist_s_ZZ, target=b_prec, sloppy=True)); # /!\ Sloppy

    _s       = sum(round(_s_ZB[_k])*_id_ZB[_k] for _k in range(len(_id_ZB)));
    assert (_s in ideal), "[Err] Minkowski pre-image is not in the target ideal (pre-image space)";
    assert (fp.fp_check_zero("T2(x)-l2(MK.x)", [_s_L.norm()-t2_norm(_s, b_prec=b_prec)], target=b_prec, sloppy=True)); # /!\ Sloppy
    
    return _s;



# ------------------------------------------------------------------------------------------
# Input/Output (to files) functions

# Want to save heavy computations related to K
# - infinite places + factor base
# - S-units(fb)

# With all these:
# + can construct freely *any* matrix for this specific factor base
# + should be sufficient to study log-su geometry (cov radius, volume, orthogonality) [only the matrix is needed]
# + solving idSVP : we still need the clDL (so we need s-units comp at each step ? or we can reuse the knowledge of su(fb) ?)


# Infinite places
# ---------------------------------------------
# Format of ( phi : K.gen() -> a) is (base 10):
#     a\n           for real places
#     Re(a) Im(a)\n for complex places
def __inf_place_out_stream(stream, K, phi):
    assert (phi.domain() == K);
    _z = phi(K.gen());

    if (_z.is_real() == True):
        stream.write("{}\n".format(_z.str(base=10)));
    else:
        stream.write("{} {}\n".format(_z.real_part().str(base=10), _z.imag_part().str(base=10)));

    return;


def __inf_place_in_stream(stream, K, to_prec=0):
    _z_reim  = [sage_eval(_s) for _s in stream.readline().rstrip().split(' ')];

    # Determine precision
    _in_prec = min(_z_reim_p.parent().precision() for _z_reim_p in _z_reim);
    assert (to_prec <= _in_prec);
    _prec   = _in_prec if (to_prec == 0) else to_prec;
    _RC     = RealField(_prec) if (len(_z_reim) == 1) else ComplexField(_prec);

    # Map input strings into RR or CC, verify it is indeed a root
    _z      = _RC(*_z_reim);
    assert (fp.fp_check_zero("K.eq(z)", [K.gen().minpoly()(_z)], target=_RC.precision()));
    
    return K.hom(_z, codomain=_RC, check=False); # Rk: Same as code of K.places()


# The first line must be: "# inf_places: nb=<n> prec=<bp>\n"
def inf_places_out_stream(stream, K, p_inf):
    stream.write("# inf_places: nb={} prec={}\n".format(len(p_inf), p_inf[0].codomain().precision()));

    for _phi in p_inf:
        __inf_place_out_stream(stream, K, _phi);
    stream.flush();

    return;


def inf_places_in_stream(stream, K):
    # Control first line, extract r1+r2, precision
    _inf_head = stream.readline().rstrip(); assert (_inf_head.startswith("# inf_places: nb="));
    _inf_head = _inf_head[len("# inf_places: "):];
    _kv   = dict(_s.split('=',1) for _s in _inf_head.split(' ')); assert (len(_kv) == 2);
    _nb   = sage_eval(_kv.get("nb")); assert (_nb == get_nb_inf_places(K));
    _prec = sage_eval(_kv.get("prec"));

    # Read _nb lines
    _p_inf = [];
    for _i in range(_nb):
        _p_inf.append(__inf_place_in_stream(stream, K, to_prec=_prec));

    assert (len(_p_inf) == _nb);
    return _p_inf;



# Finite places (factor base)
# ---------------------------------------------
# Format ideal <p, g> is, for x=K.gen(), g=a_0 + a_1 x + ... + a_d x^d:
#     p [a_0,a_1,...,a_d]\n (no padding with trailing zeroes from d to deg(K))
def __fb_pid_out_stream(stream, K, pid):
    _g1, _g2 = pid.gens_two();
    stream.write("{} {}\n".format(_g1, str(_g2.polynomial().list()).replace(' ', '')));
    
    return;


def __fb_pid_in_stream(stream, K):
    _g1, _g2 = [sage_eval(_s) for _s in stream.readline().rstrip().split(' ')];
    # Little 'Sagerie', it works for cyclotomic fields, but generic fields need padding with [0]'s.
    _g2  = _g2 + [0]*(K.degree()-len(_g2));
    _pid = K.ideal(K(_g1),K(_g2));

    assert (_pid.is_prime());
    return _pid;


# The first line must be: # fb_places: k=<n>\n
def fb_out_stream(stream, K, fb):
    stream.write("# fb_places: k={}\n".format(len(fb)));

    for _pid in fb:
        __fb_pid_out_stream(stream, K, _pid);
    stream.flush();

    return;


def fb_in_stream(stream, K):
    # Control first line, extract r1+r2, precision
    _fb_head = stream.readline().rstrip(); assert (_fb_head.startswith("# fb_places: k="));
    _fb_head = _fb_head[len("# fb_places: "):];
    _kv   = dict(_s.split('=',1) for _s in _fb_head.split(' ')); assert (len(_kv) == 1); # Overkill
    _k   = sage_eval(_kv.get("k"));

    # Read _nb lines
    _fb = [];
    for _i in range(_k):
        _fb.append(__fb_pid_in_stream(stream, K));

    assert (len(_fb) == _k);
    return _fb;


# Same functions to output/input challenges
def chal_out_stream(stream, K, chals, bsz):
    stream.write("# chals: k={} bsz={}\n".format(len(chals), bsz));

    for _chal in chals:
        __fb_pid_out_stream(stream, K, _chal);
    stream.flush();

    return;


def chal_in_stream(stream, K):
    # Control first line, extract r1+r2, precision
    _chal_head = stream.readline().rstrip(); assert (_chal_head.startswith("# chals: k="));
    _chal_head = _chal_head[len("# chals: "):];
    _kv   = dict(_s.split('=',1) for _s in _chal_head.split(' ')); assert (len(_kv) == 2); # Overkill
    _k   = sage_eval(_kv.get("k")); # ignore key "bsz"

    # Read _nb lines
    _chal = [];
    for _i in range(_k):
        _chal.append(__fb_pid_in_stream(stream, K));

    assert (len(_chal) == _k);
    return _chal;



# S-units (corresponding to previous places)
# ---------------------------------------------
# Format for x=K.gen(), su=a_0 + a_1 x + ... + a_d x^d:
#     [a_0,a_1,...,a_d]\n (no padding with trailing zeroes from d to deg(K))
def __su_out_stream(stream, K, su):
    stream.write("{}\n".format(str(su.polynomial().list()).replace(' ', '')));
    return;


def __su_in_stream(stream, K):
    _line = stream.readline();
    if not _line:
        return K(0); # for uncomplete cldl files
    _su = sage_eval(_line.rstrip());
    # Little 'Sagerie', it works for cyclotomic fields, but generic fields need padding.
    _su  = K(_su + [0]*(K.degree()-len(_su)));
    return _su;


def sunits_out_stream(stream, K, u, su):
    stream.write("# sunits(fb): nu={} k={}\n".format(len(u), len(su)));
    assert (len(u) == get_rank_units(K));
    
    for _su in u+su:
        __su_out_stream(stream, K, _su);
    stream.flush();

    return;


def sunits_in_stream(stream, K):
    # Control first line, extract r1+r2, precision
    _su_head = stream.readline().rstrip(); assert (_su_head.startswith("# sunits(fb): nu="));
    _su_head = _su_head[len("# sunits(fb): "):];
    _kv   = dict(_s.split('=',1) for _s in _su_head.split(' ')); assert (len(_kv) == 2);
    _nu   = sage_eval(_kv.get("nu")); assert (_nu == get_rank_units(K));
    _k    = sage_eval(_kv.get("k"));

    # Read units
    _u = [];
    for _i in range(_nu):
        _u.append(__su_in_stream(stream, K));
        assert (_u[-1].norm().abs() == 1);
    # Read S-unitso
    _su = [];
    for _i in range(_k):
        _su.append(__su_in_stream(stream, K));

    return _u, _su;


# Reading ClDL solution works the same way as reading S-units
def cldl_in_stream(stream, K):
    _cldl_head = stream.readline().rstrip(); assert (_cldl_head.startswith("# cldl(../data/"));
    _cldl_head = _cldl_head[len("# cldl("):];
    _fchal     = _cldl_head.split("): ")[0];
    _kv   = dict(_s.split('=',1) for _s in _cldl_head.split(' ')[1:]); assert (len(_kv) == 3);
    _k    = sage_eval(_kv.get("k")); # Ignore bsz
    _typ  = _kv.get("typ");
    
    # Read cldl elts
    _cldl = [];
    for _i in range(_k):
        _c = __su_in_stream(stream, K);
        if (_c == K(0)):
            break;
        _cldl.append(_c);
    
    return _cldl, _fchal, _typ;



# Number field 
# ---------------------------------------------------
def nf_out_stream(stream, K):
    stream.write("# number_field: tag={} eq={}\n".format(nf_get_tag(K), str(K.defining_polynomial().list()).replace(' ', '')));
    return;


def nf_in_stream(stream):
    # Control first line, extract r1+r2, precision
    _nf_head = stream.readline().rstrip(); assert (_nf_head.startswith("# number_field: tag="));
    _nf_head = _nf_head[len("# number_field: "):];
    _kv   = dict(_s.split('=',1) for _s in _nf_head.split(' ')); assert (len(_kv) == 2);
    _tag  = _kv.get("tag");
    _eq   = sage_eval(_kv.get("eq"));

    _K    = nf_set_tag(_tag, _eq);
    return _K;


# Precomputation files
# -----------------------------------------------------
# Complete file (K / FB / SU / INF)
def pcmp_out_data(filename, K, p_inf, fb, u, su):
    _f_out = open(filename, "w");
    nf_out_stream(_f_out, K);

    # Pcmp: places (inf+fb), sunits, in whatever order
    fb_out_stream(_f_out, K, fb);
    sunits_out_stream(_f_out, K, u, su);
    inf_places_out_stream(_f_out, K, p_inf); # At the end because it has no interest and is very verbose

    _f_out.close();
    return;


def pcmp_read_data(filename):
    _f_in = open(filename, "r");

    _K      = nf_in_stream(_f_in);
    _fb     = fb_in_stream(_f_in, _K);
    _u, _su = sunits_in_stream(_f_in, _K);
    _p_inf  = inf_places_in_stream(_f_in, _K);

    _f_in.close();
    return _K, _p_inf, _fb, (_u, _su);


# FB only file
#  nf / fb (head + list)
def fb_read_data(filename, K):
    _f_in = open(filename, "r");
    _K    = nf_in_stream(_f_in); assert (_K == K);
    _fb   = fb_in_stream(_f_in, K);
    
    _f_in.close();
    return _fb;


# SU only file
# nf / su (head + list)
def sunits_read_data(filename, K):
    _f_in   = open(filename, "r");
    _K      = nf_in_stream(_f_in); assert (_K == K);
    _u, _su = sunits_in_stream(_f_in, K);

    _f_in.close();
    return _u, _su;


# inf places only file
# nf / p_inf (head+list)
def inf_places_read_data(filename, K):
    _f_in = open(filename, "r");
    _K    = nf_in_stream(_f_in); assert (_K == K);
    _p_inf = inf_places_in_stream(_f_in, K);
    assert (len(_p_inf) == get_nb_inf_places(K));

    _f_in.close();
    return _p_inf;


# challenges only
# nf / chals (head + list)
def chal_read_data(filename, K):
    _f_in  = open(filename, "r");
    _K     = nf_in_stream(_f_in); assert (_K == K);
    _chals = chal_in_stream(_f_in, K);

    _f_in.close();
    return _chals;


# ClDL solutions
# nf / cldls (head+list)
def cldl_read_data(filename, K):
    _f_in  = open(filename, "r");
    _K     = nf_in_stream(_f_in); assert (_K == K);
    _cldls, _chname, _FBtyp = cldl_in_stream(_f_in, K);

    _f_in.close();
    return _cldls, _chname, _FBtyp;
