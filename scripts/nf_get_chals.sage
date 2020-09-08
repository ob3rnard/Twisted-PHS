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
import random


if len(sys.argv) != 4:
    print("Usage: {:s} <nf_tag> <bit_size> <nb>, nf_tag is :\n\tz<m> for Cyclotomic of conductor m\n\tn<p> for NTRU Prime fields of degree p\n".format(sys.argv[0]));
    sys.exit(2);

tag      = sys.argv[1];
chal_bsz = Integer(sys.argv[2]);
nb_chal  = Integer(sys.argv[3]);


# --------------------------------------------------------------------------------------
# Obtain number field
K     = nf_set_tag(tag);
print ("{}: get challenges (bit_size={},N={})".format(tag, chal_bsz, nb_chal), flush=True);
cyclo = True if tag[0] == "z" else False;


# --------------------------------------------------------------------------------------
# Output file names
data_dir    = "./data/"+tag+"/";
chal_name   = "{}.chal_b{}_n{}".format(data_dir+tag, chal_bsz, nb_chal);


# --------------------------------------------------------------------------------------
# Find random challenge of bsz
def rand_chal(K, bsz, is_cyclo):
    S = [];
    while (len(S) == 0):
        p = Integer(random.getrandbits(bsz));
        if (is_cyclo == True):
            p = K.next_split_prime(p);
            S = K.primes_above(p);
        else:
            p = next_prime(p); # Not ideal, but doesn't matter for getting challenges
            S = [pid for pid in K.primes_above(p) if pid.norm() == p];

    return random.choice(S); # Ramified primes should not occur statistically

print ("Generating challenges of prime norm...", flush=True, end='');
_t = cputime(); chals = [ rand_chal(K, chal_bsz, cyclo) for _i in range(nb_chal) ]; _t = cputime(_t);
print ("\t[done] t={:.2f}".format(_t), flush=True);


# --------------------------------------------------------------------------------------
# Write to file
_f_out = open(chal_name, "w");
print ("Output in {}".format(chal_name), flush=True, end='');
nf_out_stream(_f_out, K);
chal_out_stream(_f_out, K, chals, chal_bsz);
_f_out.write("# --- END ---\n");
_f_out.close();
print("\t[done]", flush=True);


exit;
