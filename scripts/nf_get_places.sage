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

if len(sys.argv) != 2:
    print("Usage: {:s} <nf_tag>, nf_tag is :\n\tz<m> for Cyclotomic of conductor m\n\tn<p> for NTRU Prime fields of degree p\n".format(sys.argv[0]));
    sys.exit(2);

tag = sys.argv[1];


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
print ("{}: get places".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# Output file names
data_dir    = "./data/"+tag+"/";
p_inf_name  = data_dir+tag+".inf";
fb_phs_name = data_dir+tag+"_phs.fb";
fb_opt_name = data_dir+tag+"_opt.fb";
fb_tw_name  = data_dir+tag+"_tw.fb";


# --------------------------------------------------------------------------------------
# Compute inf places at high precision
HPREC = 5000;
_t = cputime(); p_inf = get_inf_places(K, b_prec=HPREC); _t = cputime(_t);
print ("Inf_places -> {} prec={} t={:.2f}".format(p_inf_name, HPREC, _t), end='', flush=True);
_f_out = open(p_inf_name, "w");
nf_out_stream(_f_out, K);
inf_places_out_stream(_f_out, K, p_inf);
_f_out.write("# --- END ---\n");
_f_out.close();
print ("\t[done]", flush=True);


# --------------------------------------------------------------------------------------
# Compute factor bases for phs / opt / tw
fb_proc  = { "phs" : phs_get_fb,
             "opt" : opt_get_fb,
             "tw"  : tw_get_fb,
};
fb_names = { "phs" : fb_phs_name,
             "opt" : fb_opt_name,
             "tw"  : fb_tw_name,
};

for typ in ["phs", "opt", "tw"]:
    _t     = cputime(); _fb    = fb_proc.get(typ)(K); _t     = cputime(_t);
    print ("FB type {}->{} k={} t={:.2f}".format(typ, fb_names.get(typ), len(_fb), _t), end='', flush=True);
    _f_out = open(fb_names.get(typ), "w");
    nf_out_stream(_f_out, K);
    fb_out_stream(_f_out, K, _fb);
    _f_out.write("# --- END ---\n");
    _f_out.close();
    print ("\t[done]", flush=True);



exit;
