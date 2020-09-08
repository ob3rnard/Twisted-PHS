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
print ("{}: compute log S-unit lattices".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# File names
data_dir    = "./data/"+tag+"/";
p_inf_name  = data_dir+tag+".inf";
fb_phs_name = data_dir+tag+"_phs.fb";
fb_opt_name = data_dir+tag+"_opt.fb";
fb_tw_name  = data_dir+tag+"_tw.fb";
su_phs_name = data_dir+tag+"_phs.su";
su_opt_name = data_dir+tag+"_opt.su";
su_tw_name  = data_dir+tag+"_tw.su";
lat_phs_name= data_dir+tag+"_phs.lat";
lat_phs0_name= data_dir+tag+"_phs0.lat";
lat_opt_name= data_dir+tag+"_opt.lat";
lat_opt0_name= data_dir+tag+"_opt0.lat";
lat_tw_name = data_dir+tag+"_tw.lat";

# Bof
typ_link = { "phs" :
             {"fb_name" : fb_phs_name,
              "su_name" : su_phs_name,
              "lat_name": lat_phs_name,
              "fct_opt" : 'PHS'},
             "phs0" : # PHS using TW factor base
             {"fb_name" : fb_tw_name,
              "su_name" : su_tw_name,
              "lat_name": lat_phs0_name,
              "fct_opt" : 'PHS'},
             "opt" :
             {"fb_name" : fb_opt_name,
              "su_name" : su_opt_name,
              "lat_name": lat_opt_name,
              "fct_opt" : 'OPT'},
             "opt0" :
             {"fb_name" : fb_tw_name,
              "su_name" : su_tw_name,
              "lat_name": lat_opt0_name,
              "fct_opt" : 'OPT'},
             "tw" :
             {"fb_name" : fb_tw_name,
              "su_name" : su_tw_name,
              "lat_name": lat_tw_name,
              "fct_opt" : 'TW'},
};


# --------------------------------------------------------------------------------------
# This part is easy, we work at (very) high precision just in case
HPREC = 5000; # We need good inf places, but we could lower the lattice precision.
SCALE = 2;


# --------------------------------------------------------------------------------------
# Load inf places
print ("Inf_places <- {} (prec:{})".format(p_inf_name, HPREC), flush=True, end='');
_f_in = open(p_inf_name, "r");
_ = nf_in_stream(_f_in); assert (_ == K);
p_inf = inf_places_in_stream(_f_in, K);
# Trust is good, control is better
assert (abs(p_inf[0].codomain().precision() - SCALE*HPREC) < 10);
assert (len(p_inf) == get_nb_inf_places(K));
_f_in.close();
print ("\t[done]", flush=True);


# --------------------------------------------------------------------------------------
# For each of 'PHS', 'OPT', 'TW', load FB + SU, compute twlat and output
for typ in ["tw", "opt0", "phs0", "opt", "phs"]: # From smallest to biggest
    print ("'{}' method".format(typ), flush=True);
    # Load fb
    f_in = open(typ_link.get(typ).get("fb_name"), "r");
    tst  = nf_in_stream(f_in); assert (tst == K);
    fb   = fb_in_stream(f_in, K);
    f_in.close();
    print ("    #FB:{}".format(len(fb)), flush=True);
    # Load su (if it exists, otherwise continue)
    if (not os.path.exists(typ_link.get(typ).get("su_name"))):
        print ("    [next] S-Unit file '{}' does not exist.".format(typ_link.get(typ).get("su_name")), flush=True);
        continue;
    f_in  = open(typ_link.get(typ).get("su_name"), "r");
    tst   = nf_in_stream(f_in); assert (tst == K);
    u, su = sunits_in_stream(f_in, K);
    f_in.close();
    print ("    #u:{}\t#su:{}".format(len(u), len(su)), flush=True);
    # Compute lat
    print ("    Log S-unit lattice... ", end='', flush=True);
    t = cputime();
    L = twphs_get_matrix(u, su, p_inf, fb, typ_link.get(typ).get("fct_opt"), b_prec=HPREC);
    t = cputime(t);
    print ("\t[done] t={:.2f}".format(t), flush=True);
    # Output lat
    t = cputime();
    lattice_out_data(typ_link.get(typ).get("lat_name"), L);
    t = cputime(t);
    print ("    output t={:.2f}".format(t), flush=True);
    
    
exit;
