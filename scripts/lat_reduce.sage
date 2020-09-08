#!/usr/bin/env sage

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------

from sage.all import *
import fp
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
print ("{}: reduce log S-unit lattices".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# Output file names
data_dir  = "./data/"+tag+"/";
list_lats = [ data_dir+tag+"_"+_s+".lat" for _s in ["tw", "opt", "phs", "opt0", "phs0"]];
out_bkz   = [ data_dir+tag+"_"+_s+".bkz" for _s in ["tw", "opt", "phs", "opt0", "phs0"]];


# --------------------------------------------------------------------------------------
# This part is easy, we work at (very) high precision just in case
MAX_PREC = 500; # Beyond 500, fplll performance drops drastically (?)
MIN_PREC = 300; # It works fine everywhere.
VOL_SCALE = 2;  # Work at 2 times the volume log in base 2
MAX_LOOPS = 300;
BLOCK_SZ  = 40;
HKZ_MAX   = 60;

# --------------------------------------------------------------------------------------
# For each of 'PHS', 'OPT', 'TW', load FB + SU, compute twlat and output
for i in range(len(list_lats)):
    lat_file = list_lats[i];
    out_file = out_bkz[i];
    if (not os.path.exists(lat_file)):
        print ("[next] Lat file '{}' does not exist.".format(lat_file), flush=True);
        continue;
    print ("Lat '{}'".format(lat_file), flush=True);

    # Input lattice
    L = lattice_read_data(lat_file);
    print ("    dim={}".format(L.nrows()), flush=True);

    # Determine precision
    logv = lnvol(L)/ln(L.base_ring()(2));
    print ("    logvol={:.2f}".format(float(logv)));
    assert (logv < MAX_PREC); # Otherwise it will fail
    work_prec = ceil(logv*VOL_SCALE);
    work_prec = MIN_PREC if work_prec < MIN_PREC else MAX_PREC if work_prec > MAX_PREC else work_prec;

    # Determine block_size
    block_sz = L.nrows() if (L.nrows() <= HKZ_MAX) else BLOCK_SZ;

    # Launching BKZ
    print ("    bkz bk={} prec={} nloops={} ...".format(block_sz, work_prec, MAX_LOOPS), flush=True, end='');
    #t = cputime(subprocesses=True);
    t = walltime();
    B = bkz(L, work_prec=work_prec, block_size=block_sz, bkzmaxloops=MAX_LOOPS)[0];
    t = walltime(t);
    print ("\t[done] t={:.2f}".format(t), flush=True);
    fp.fp_check_zero("vol(L)=vol(BKZ(L))", [lnvol(B)-lnvol(L)], target=work_prec, sloppy=True);
    
    # Output
    lattice_out_data(out_file, B);
    print ("    out -> '{}'".format(out_file), flush=True);

exit;
