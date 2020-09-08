#!/usr/bin/env sage

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------

from sage.all import *
from lattice import *
from number_field import *
from twphs_algo import *

if len(sys.argv) != 4:
    print("Usage: {:s} <nf_tag> <bsz> <N>, nf_tag is :\n\tz<m> for Cyclotomic of conductor m\n\tn<p> for NTRU Prime fields of degree p\n    bsz and N designates the challenge file ../data_dir/<nf>/<nf>.chal_b<bsz>_n<N>".format(sys.argv[0]));
    sys.exit(2);

tag = sys.argv[1];
bsz = sage_eval(sys.argv[2]);
N   = sage_eval(sys.argv[3]);


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
print ("{}: solve Approx-idSVP".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# File names
data_dir    = "./data/"+tag+"/";
chal_name   = (data_dir + tag + ".chal_b{}_n{}").format(bsz, N); assert (os.path.exists(chal_name));
p_inf_name  = data_dir+tag+".inf"; assert (os.path.exists(p_inf_name));

# Typs : "tw", "opt0", "phs0", "opt", "phs"
typ_list    = ["tw", "opt0", "phs0", "opt", "phs"];
typ_opts    = ['TW', 'OPT',  'PHS',  'OPT', 'PHS'];
fb_names    = [ data_dir+tag+"_"+("tw" if (len(typ) == 4) else typ)+".fb" for typ in typ_list];
assert (all(os.path.exists(_f) for _f in fb_names));
su_names    = [ data_dir+tag+"_"+("tw" if (len(typ) == 4) else typ)+".su" for typ in typ_list];
assert (all(os.path.exists(_f) for _f in su_names));
cldl_names  = [ (data_dir+tag+ "_{}.cldl_b{}_n{}").format(("tw" if (len(typ) == 4) else typ), bsz, N) for typ in typ_list ];
assert (all(os.path.exists(_f) for _f in cldl_names));

lat_names   = [ data_dir+tag+"_"+typ+".lat" for typ in typ_list];
assert (all(os.path.exists(_f) for _f in lat_names));
bkz_names   = [ data_dir+tag+"_"+typ+".bkz" for typ in typ_list];
assert (all(os.path.exists(_f) for _f in bkz_names));


# --------------------------------------------------------------------------------------
# This part is easy, we work at (very) high precision just in case
HPREC = 5000; # We need good inf places, but we could lower the lattice precision.
SCALE = 2;


# --------------------------------------------------------------------------------------
# Load challenges
print ("Chals <- {}".format(chal_name), flush=True);
chal = chal_read_data(chal_name, K);

# --------------------------------------------------------------------------------------
# Load inf places
print ("Inf_places <- {} (prec:{})".format(p_inf_name, HPREC), flush=True);
p_inf = inf_places_read_data(p_inf_name, K);
# Trust is good, control is better
assert (abs(p_inf[0].codomain().precision() - SCALE*HPREC) < 10);

# --------------------------------------------------------------------------------------
# Load factor bases / Sunits / Lattices
print ("Load fbs: {}".format(fb_names), flush=True);
fb     = [ fb_read_data(_f_name, K) for _f_name in fb_names ];
print ("Load sunits: {}".format(su_names), flush=True);
u_su   = [ sum(sunits_read_data(_f_name, K), []) for _f_name in su_names ];
print ("Load cldls: {}".format(cldl_names), flush=True);
cldl   = [ cldl_read_data(_f_name, K)[0] for _f_name in cldl_names ];
print ("Load lattices: {}".format(lat_names), flush=True);
lsulat = [ lattice_read_data(_f_name) for _f_name in lat_names ];
print ("Load bkzlat: {}".format(bkz_names), flush=True);
lsubkz = [ lattice_read_data(_f_name) for _f_name in bkz_names ];
print ("[done] Precmp loaded.", flush=True);


# --------------------------------------------------------------------------------------
# Extend inf places to precision sufficient for cldl coefficients
#max_log = max([max([ max([max( RealField(1000)(log(_coef.numer().abs(),2)), RealField(1000)(log(_coef.denom().abs(),2)))  ] for _coef in _c.list()) for _c in _cldl]) for _cldl in cldl]);

# The max is 25128 for PHS z53
CLDL_PREC = 25500;
print ("Extending p_inf for clDL to prec={}...".format(CLDL_PREC), flush=True, end='');
t = cputime();
p_inf = extend_inf_places(K, p_inf, to_prec=CLDL_PREC);
t = cputime(t);
print ("\t[done] t={:.2f}".format(t), flush=True);


# --------------------------------------------------------------------------------------
# Pre compute fHcE
r1, r2 = K.signature();
print ("Precompute fHcE matrices...", flush=True, end='');
t = cputime();
fHcE_tab = [ get_twfHcE_matrix(r1, r2, fb[typ], typ_opts[typ], b_prec=lsubkz[typ].base_ring().precision()) for typ in range(len(typ_list)) ];
t = cputime(t);
print ("\t\t[done] t={:.2f}".format(t));


# --------------------------------------------------------------------------------------
# Output file
out_file = (data_dir + tag + ".sol_b{}_n{}").format(bsz, N);
f_out    = open(out_file, "w");
f_out.write("# nf:'{}' chal:'{}'\n".format(tag, chal_name));
f_out.write(("# exact"+"\t{}"*len(typ_list)+"\n").format(*typ_list));


# ----------------------------------------------------------------------------------------
# Tw-query
def test_protocol(_typ_idx, a, cldl_k, u_su_typ, p_inf, fb_typ, fhce_typ, lsulat_typ, lsubkz_typ, nb_iter):
    # guess beta
    beta_inf  = 0;
    beta_cur  = twphs_guess_beta(a, lsubkz_typ, K, fb_typ, typ_opts[_typ]); assert(beta_cur > 0);
    beta_sup  = 2*beta_cur;
    found_one = False;
    
    best_s = K(a.norm());
    for _i in range(nb_iter):
        print ("       beta={:5.2f}".format(float(beta_cur)), flush=True, end='');
        _t = cputime();
        t = twphs_get_target(cldl_k, a, p_inf, fb_typ, typ_opts[_typ], beta=beta_cur,
                             b_prec=lsubkz_typ.base_ring().precision(), _pcmp_fhce=fhce_typ);
        v = cvp_babai_NP(lsubkz_typ, t);
        s = twphs_build_solution(cldl_k, v, lsulat_typ, u_su_typ);
        _t = cputime(_t);
        print ("\t[{}]\tl2={:7.3f} t={:.2f}\t[alg_norm={}]".format(s in a, float(t2_norm(s)), _t, s.norm().factor()));

        if (s in a):
            found_one = True;
            beta_sup = beta_cur;
            beta_cur = (beta_sup+beta_inf)/2;
            best_s = s if (t2_norm(s) < t2_norm(best_s)) else best_s;
        else:
            beta_inf = beta_cur;
            beta_cur = (beta_sup+beta_inf)/2 if (found_one == True) else 2*beta_cur;
    
    return best_s;



# --------------------------------------------------------------------------------------
# For each challenge...
NB_ITER = 5; # For dichotomy on beta.
for _k in range(len(chal)):
    print ("Challenge #{}:".format(_k));
    _a = chal[_k];
    
    # Exact id-SVP (up to 47, and... for 53 ?)
    # ----------------------------------------
    _mink_prec = max(ceil((log( K.discriminant().abs(), 2) + log(_a.norm(), 2))+50), 500);
    print ("    exact id-svp using prec={}...".format(_mink_prec), flush=True, end='');
    t = cputime();
    if (K.degree() < 52):
        _s = idsvp_exact(_a, b_prec=_mink_prec);
    else:
        _s = idsvp_exact(_a, b_prec=_mink_prec, approx_bkz=True); # exact svp doesn't finish on z53
    t = cputime(t);
    print ("\t[done] t={:.2f} s={:5.2f}".format(t, float(t2_norm(_s))), flush=True);
    f_out.write("{:7.2f}".format(float(t2_norm(_s, b_prec=500)))); f_out.flush();

    # Running each Approx-SVP algorithm
    for _typ in range(len(typ_list)):
        print ("    method:'{}':".format(typ_list[_typ]), flush=True);
        # tester si cldl
        if (_k > len(cldl[_typ])):
            print ("        [next] No cldl", flush=True);
            f_out.write("\t0"); f_out.flush();
        # Test protocol: typ_index, chal, cldl, u+su, p_inf, fb, nb_iter
        _t  = cputime();
        _s  = test_protocol(_typ, _a, cldl[_typ][_k], u_su[_typ], p_inf, fb[_typ], fHcE_tab[_typ], lsulat[_typ], lsubkz[_typ], NB_ITER);
        _t  = cputime(_t);
        assert(_s in _a);
        f_out.write("\t{:7.2f}".format(float(t2_norm(_s, b_prec=500)))); f_out.flush();
        print ("    [Tot.] t={:.2f}".format(_t), flush=True);
        
    f_out.write("\n"); f_out.flush();

f_out.write("# --- END --- \n");
f_out.close();


exit;
