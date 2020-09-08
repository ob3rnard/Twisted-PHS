from sage.all import *
import fp


# ----------------------------------------------------------------------------------------------
# General "matrix on real" caveats handling
def rc_mat_inverse(M):
    _R      = M.base_ring(); # Real or Complex
    _b_prec = _R.precision();
    _RR     = _R.to_prec(4*_b_prec); # Should be enough, 1 is not sufficient, 2/3 maybe

    _iM     = ((M.change_ring(_RR)).inverse()).change_ring(_R);

    # Check there is no precision issue (NB: compute the product in high precision)
    _chk_M  = (M.change_ring(_RR)*_iM.change_ring(_RR)).change_ring(_R) - identity_matrix(M.nrows());
    assert (fp.fp_check_zero("M*iM-In", _chk_M.coefficients(), target=_b_prec, sloppy=True)); # /!\ Sloppy

    return _iM;


# Conversions RR <--> ZZ
# Scaling factors are given in log_2 to match "precision" of RR

# Get integer equivalent of (M) [ returns floor(2^prec*M) on R ]
# Returns also the log_2 of the applied scaling factor.
def matvec_real_to_ZZ(Mv, work_prec=0):
    _R         = Mv.base_ring().fraction_field();
    _log_scale = _R.precision() if (work_prec == 0) else work_prec;
    _scale     = Integer(2)**(_log_scale);
    _MZ        = Mv.apply_map(lambda _mij:Integer(floor(_scale*_mij)));

    return _MZ, _log_scale;


# Returns to real numbers
def matvec_ZZ_to_real(Mv, log_scale, to_prec=0):
    assert ((to_prec == 0) or log_scale <= to_prec);

    _prec = log_scale if (to_prec == 0) else to_prec;
    _R    = RealField(_prec);
    _MR   = Mv.apply_map(lambda _mij:_R(_mij)/_R(2)**(log_scale));
    
    assert (_MR.base_ring().precision() >= _prec);
    return _MR;



# ----------------------------------------------------------------------------------------------
# Gram-Schmidt Orthogonalization
# Input: M is m*n matrix on any Ring
# Return G,P: G is GSO of M (normalized if normalize=True), and permutation matrix (of the rows)
#             P is the transition matrix such that M = P*G
# NB: annoyingly, sage does not implement this outside of ZZ/QQ, RDF and CDF.
def gram_schmidt_ortho(M, normalize=False):
    _R     = M.base_ring().fraction_field();
    _n     = M.nrows();
    
    _G = M.change_ring(_R);
    _P = identity_matrix(_R, _n);
    
    # Main loop
    # NB: Exchanging _i and _j loop would lead to somewhat "order-independent" algorithm,
    #     allowing to choose the smallest length projection for the next step
    #    for _i in range(1,_n):
    #        _G[_i] = _G[_i] - sum( (_G[_i]*_G[_j])/_G[_j].norm()**2 * _G[_j] for _j in range(_i));
    for _i in range(1,_n):
        _mu_i  = [(_G[_i]*_G[_j])/_G[_j].norm()**2 for _j in range(_i)] + [1] + [0]*(_n-_i-1);
        _G[_i] = _G[_i] - sum(_mu_i[_j] * _G[_j] for _j in range(_i));
        _P[_i] = vector(_R, _mu_i);#[_mu_i[_j] for _j in range(_i)] + [_R(1)] + [_R(0)]*(_n-_i-1));
    
    # Orthonormalization (not by default)
    if (normalize == True):
        for _i in range(_n):
            _norm_i   = _G[_i].norm();
            _P[:,_i] *= _norm_i;
            _G[_i]   /= _norm_i;

    assert (_G.base_ring() == M.base_ring().fraction_field());
    assert (fp.fp_check_zero("M-PG", (M-_P*_G).coefficients(), target=_R.precision())); # **Warn** This assertion is not costless
    return _G, _P;



# ----------------------------------------------------------------------------------------------
# Lattice constants (reduced volume, Rankin-d, etc)

# Return Vol M
def lnvol(M):
    assert (M.nrows() <= M.ncols());
    _gM = gram_schmidt_ortho(M)[0];
    return sum([ ln(_gM[_k].norm()) for _k in range(M.nrows())]);
def vol(M):
    return exp(lnvol(M));


# Return |Vol M|^(1/n), where n is the matrix (the number of rows)
def lnvol_reduced(M):
    return (1/M.base_ring()(M.nrows()))*lnvol(M);
def vol_reduced(M):
    return exp(lnvol_reduced(M));


# Orthogonality defect
def lnrankin(M,d):
    assert (d <= M.nrows());
    _ln_rk_d   = sum([ln(M[_k].norm()) for _k in range(d)]);
    _ln_rk_d   = _ln_rk_d - lnvol(M)*M.base_ring()(d/M.nrows());
    return _ln_rk_d;
def rankin(M,d):
    return exp(lnvol(M));

def lnrankin_reduced(M,d):
    return (1/M.base_ring()(M.nrows()))*lnrankin(M,d);
def rankin_reduced(M, d):
    return exp(lnrankin_reduced(M,d));


# Minimum angle
def min_angle(M):
    _pi     = M.base_ring()(pi);
    _min_th = _pi/M.base_ring()(2);
    for _i in range(M.nrows()):
        _th_i = _pi/M.base_ring()(2);
        for _j in range(_i+1, M.nrows()):
            _theta  = arccos( (M[_i]*M[_j])/M[_i].norm()/M[_j].norm() );
            _th_ij  = min(_theta, _pi-_theta);
            _th_i   = min(_th_i, _th_ij);
        _min_th = min(_min_th, _th_i);

    return _min_th;

def mean_angle(M):
    _pi     = M.base_ring()(pi);
    _s      = 0;
    _N      = M.base_ring()(M.nrows()*(M.nrows()-1))/M.base_ring()(2);
    for _i in range(M.nrows()):
        for _j in range(_i+1, M.nrows()):
            _theta  = arccos( (M[_i]*M[_j])/M[_i].norm()/M[_j].norm() );
            _th_ij  = min(_theta, _pi-_theta);
            _s      = _s + _th_ij;
    return (_s/_N);



# ------------------------------------------------------------------------------------------
# Projection / full-rank isometry wizardry 

# Matrix fH of article §4. Pour §3, use with k=0
def get_minkH(r1, r2, k, b_prec=fp.BIT_PREC_DEFAULT):
    _dim   = r1+r2+k;
    _B_one = matrix([ vector([1 if (_i == _j) else (-1 if (_i == _j+1) else 0) for _j in range(_dim)])
                      for _i in range(1,_dim) ]);
    # To Minkowski space
    _B_c   = block_matrix(r2, [matrix([0,0]*i+[1/2,1/2]+[0,0]*(r2-i-1)) for i in range(r2)],
                          subdivide=False);
    _B_m   = block_diagonal_matrix([identity_matrix(r1), _B_c, identity_matrix(k)],
                                   subdivide=False);

    # Isometry matrix
    _mH    = matrix(RealField(b_prec), _B_one*_B_m);
    _OmH,_ = gram_schmidt_ortho(_mH, normalize=True);

    # Verify prec
    _OmH   = _OmH.transpose();
    assert (_OmH.base_ring().precision() >= b_prec);
    return _OmH;



# Projection on H = <1,...,1> (k times)
# ----------------------------------------------
# Returns the projection matrix on the orthogonal of <1,...,1> (k times)
def get_projection_H(n, k, b_prec=fp.BIT_PREC_DEFAULT): 
    _pH  = (-1/RealField(b_prec)(k))*ones_matrix(k) + identity_matrix(k);
    _pHi = block_diagonal_matrix(_pH, identity_matrix(n-k), subdivide=False);

    assert (_pHi.base_ring().precision() >= b_prec);
    return _pHi;



# Pruning
# ------------------------------------------
# Prune r2 coordinates + the last one
def get_ext_pruning(r1, r2):
    _mc   = matrix([[1],[0]]);
    _Mc   = block_diagonal_matrix([*[_mc]*r2], subdivide=False);
    _pr   = block_diagonal_matrix([identity_matrix(r1),_Mc], subdivide=False);
    _pr   = _pr.matrix_from_columns(range(_pr.ncols()-1));
    
    return _pr;



# ------------------------------------------------------------------------------------------
# FpLLL helpers
# The Sage version of fplll has some hard-to-track floating-point issues btw versions,
# so the best is probably to call directly fplll via system calls.
import subprocess

# Path details
__FPLLL_PATH = "/usr/local/bin/";

def fplll_get_path():
    return __FPLLL_PATH;

def fplll_change_path(new_path):
    global __FPLLL_PATH;
    __FPLLL_PATH = new_path;


# System machinery for external fplll calls
# Consistency is not verified [eg. asking for BKZ with a non empty vector is going to fail] 
# Final parsing is correct under the condition that matrices are output in a compatible Sage format (-of uk, -of bk, -of vk)
def __fplll_call(act, B, t=[], opts=""):
    _t_fname = tmp_filename();
    _t_file  = open(_t_fname, 'w');
    _t_file.write("[{0}]\n".format("\n".join([str(list(i)).replace(",","") for i in B.rows()])));
    if (len(t) != 0):
        _t_file.write("{0}\n".format(str(list(t)).replace(",", "")));
    _t_file.close(); # Automatically erased when Sage exits

    # Requires Python >= 3.7
    #print ("cmd: env -i {0}./fplll -a {1} {2} {3}".format(__FPLLL_PATH, act, opts, _t_fname));
    proc = subprocess.run(["env -i {0}./fplll -a {1} {2} {3}".format(__FPLLL_PATH, act, opts, _t_fname)],
                          shell=True, stdout=subprocess.PIPE, text=True);
    # Valid if output is with commas (Sage compliant), otherwise matrices will be torn apart.
    out  = proc.stdout.rstrip().replace(",\n", ",").split("\n");
    
    assert (out[0] != '[]'); # Empty output can hide a precision problem.
    return out;


def __fplll_call_cvp(B, t):
    # Weird bug? Output '[]' with -of c, and some vector without it. Doc says "-of c (default if -a cvp)".
    [_v_str] = __fplll_call("cvp", B, t=t, opts="");
    _v       = sage_eval("vector({0})".format(_v_str.replace(" ", ",")));
    return _v;


def __fplll_call_svp(B):
    [_s_str] = __fplll_call("svp", B, opts="-of s");
    _s       = sage_eval("vector({0})".format(_s_str.replace(" ", ",")));
    return _s;


# Official BKZ options
# -b: block_size -f: float_type (eg. mpfr) -p: precision
# -bkzmaxloops: int -bkzmaxtime: sec -bkzautoabort
# -s: strategy file -bkzghbound:... -bkzboundedlll -bkzdumgso filename
def __fplll_call_bkz(B, block_size=0, **kwargs):
    _bck_size = B.nrows() if (block_size == 0) else block_size;
    _opt_str  = "-b {0}".format(_bck_size);
    # Flexibility to add any arguments that are legitimate for fplll -a bkz or fplll -a hkz
    for _key, _value in kwargs.items():
        if _value == True:
            _opt_str += " -"+_key;
        else:
            _opt_str += " -"+_key+" "+str(_value);
            # Output format
    _opt_str += " -of bkuk";

    # BKZ call, _u is st. _u*B = _bkz
    _bkz_str, _u_str = __fplll_call("bkz", B, opts=_opt_str);
    _bkz     = sage_eval("matrix({0})".format(_bkz_str));
    _u       = sage_eval("matrix({0})".format(_u_str)); # _u * B = _bkz

    return _bkz, _u;



# -------------------------------------------------------------------------------------------
# Reduction : return always B, U st. B=U*M
def bkz_ZZ(MZ, block_size=0, **kwargs):
    return __fplll_call_bkz(MZ, block_size=block_size, **kwargs);

# Returns BKZ reduction of M
# Options: see __fplll_call_bkz() doc
def bkz(M, work_prec=0, block_size=0, **kwargs):
    # Wrap
    _MZ, _l = matvec_real_to_ZZ(M, work_prec=work_prec);
    # Reduce
    _bkz_Z, _uZ = __fplll_call_bkz(_MZ, block_size=block_size, **kwargs);
    # Unwrap
    _bkz    = matvec_ZZ_to_real(_bkz_Z, _l);

    return _bkz, _uZ;



# ------------------------------------------------------------------------------------------
# Solving SVP/CVP

# Exact SVP
def svp_exact_ZZ(MZ):
    return __fplll_call_svp(MZ);

def svp_exact(M, work_prec=0):
    # Wrap
    _MZ, _l = matvec_real_to_ZZ(M, work_prec=work_prec);
    # Call FPLLL SVP
    _svpZ = __fplll_call_svp(_MZ); 
    # Unwrap
    _svp    = matvec_ZZ_to_real(_svpZ, _l);

    return _svp;


# Babai
# NB: This follows [Gal, \S18.1, Alg.26], but significantly differs from eg. the code of [DPW19].
def cvp_babai_NP(M, t, G=0):
    if (G == 0):
        _G, _ = gram_schmidt_ortho(M, normalize=False);
    else:
        _G = G;
    
    _w = t;
    # Babai's Nearest Plane, see book [Gal, \S18.1, Alg.26]
    _n = M.nrows();
    _v = zero_vector(M.ncols());
    for _i in range(_n-1, -1, -1):
        _z = (_w*_G[_i])/_G[_i].norm()**2;
        _v = _v + round(_z)*M[_i]; # (Permutation needed here if GSO permutes rows)
        _w = _w - (_z - round(_z))*_G[_i] - round(_z)*M[_i]; # (and here)

    assert (_v.base_ring() ==  M.base_ring());
    return _v;


# Exact CVP
# The non-ZZ routine does take t back towards 0 using Babai NP before solving CVP(L,t),
def cvp_exact_ZZ(MZ, tZ):
    # Babai NP ?
    return __fplll_call_cvp(MZ, tZ);

def cvp_exact(M, t, work_prec=0):
    # Clean target
    _v = cvp_babai_NP(M, t);
    _t = t - _v;
    # Wrap
    _MZ, _l = matvec_real_to_ZZ(M,  work_prec=work_prec);
    _tZ, _  = matvec_real_to_ZZ(_t, work_prec=_l);
    # Call FPLLL CVP
    _cvpZ   = __fplll_call_cvp(_MZ, _tZ);
    # Unwrap (ZZ -> RR + Babai's NP)
    _cvp    = matvec_ZZ_to_real(_cvpZ, _l);
    _cvp    = _cvp + _v;
    
    return _cvp;



# ------------------------------------------------------------------------------------------
# Input/Output (to files) functions

# Format of the matrix output is defined as follow (to be compatible with Magma functions)
# [\n
# \t[<b_0>],\n
# ... (one line per b_i)
# \t[<b_n>],\n
# ]\n
# where each b_i is a comma separated list of elements ' c0, c1, ..., ck ',
# and each c0 is [-]d(ddd...).f(fff...)[Ee(eee...)]
def lattice_out_data(filename, L):
    _f_out = open(filename, "w");
    # Handle the case L is of type lattice instead of matrix ?
    # It was a problem in Magma, not sure it is in Sage
    # L = matrix(..., L);

    _L_str = ("[\n\t"
              + ",\n\t".join(map(str, [list(_row) for _row in L.rows()]))
              + ",\n]\n").replace('e','E');
    _f_out.write(_L_str);
    _f_out.close();

    return;


def lattice_read_data(filename, to_b_prec=0):
    _f_in = open(filename, "r");
    _M_ij = eval(preparse(_f_in.read())); # contains list of list of coeffs
    _f_in.close();

    # Get input precision
    _in_b_prec = max([_Mi[_j].parent().precision() for _Mi in _M_ij for _j in range(len(_Mi))]);
    # input precision might be way too high for both practical computations and usefulness
    # so a shrinking mechanism is provided (no control on _to_prec < _in_prec, though)
    _b_prec    = _in_b_prec if ((to_b_prec == 0) or (_in_b_prec < to_b_prec)) else to_b_prec;
    
    # Construct matrix
    _R = RealField(_b_prec);
    _M = matrix(_R, _M_ij);
    
    return _M;
