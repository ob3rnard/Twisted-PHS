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
print ("{}: eval Gram-Schmidt log norms of log S-unit lattices".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# Lattice file names: <data_dir>/<tag>/<tag>_<typ>.lat / <data_dir>/<tag>/<tag>_<typ>.bkz
data_dir = "./data/"+tag+"/";
typ_list = ["tw", "opt0", "phs0", "opt", "phs"]; # typ0 versions correspond to algo typ using twFB



# --------------------------------------------------------------------------------------
# Output one file per type
# Each file contain 5 columns: i ln(bi*) ln(ci*) ln(bi*)/ln(bi-1*) ln(ci*)/ln(ci*-1) 
# bi are GS vectors before BKZ, ci are GS vectors after BKZ
def print_data_line(f_out, data):
    RRout = RealField(106);
    #f_out.write(("{:d}"+"\t{}"*4+"\n").format(data[0], *[RRout(_d) for _d in data[1:5]]));
    f_out.write(("{:d}"+"\t{:14.12f}"*4+"\n").format(data[0], *[float(_d) for _d in data[1:5]]));
    f_out.flush();


for typ in typ_list:
    print ("'{}' method".format(typ), flush=True);

    # Reading lattices
    L_file = data_dir + tag + "_" + typ + ".lat";
    B_file = data_dir + tag + "_" + typ + ".bkz";
    if (not os.path.exists(L_file)):
        print ("    [next] Log S-unit lat file '{}' does not exist.".format(L_file), flush=True);
        continue;
    if (not os.path.exists(B_file)):
        print ("    [next] Bkz file '{}' does not exist.".format(B_file), flush=True);
        continue;
    print("    Reading data...", flush=True, end='');
    B = lattice_read_data(B_file);
    work_prec = B.base_ring().precision();
    L = lattice_read_data(L_file);
    print("\t[done]", flush=True);

    # Compute Gram-Schmidt
    print ("    Gram-Schmidt computation...", end='', flush=True);
    t = cputime();
    GSB, _ = gram_schmidt_ortho(B, normalize=False);
    GSL, _ = gram_schmidt_ortho(L, normalize=False);
    t = cputime(t);
    print ("\t[done] t={:.2f} (for 2)".format(t), flush=True);
    
    # Output file and header
    out_name = data_dir + tag + "_" + typ + ".gsn";
    f_out    = open(out_name, "w");
    f_out.write("# nf:'{}' in_nobkz='{}' in_bkz='{}'\n".format(tag, L_file, B_file));
    f_out.write("# i\tln(ai*)\t\tln(bi*)\t\tlnarat\t\tlnbrat\n");
    f_out.flush();

    # --------------------------------------------------------------------------------------
    # Read matrix: which prec ?
    # Row, Norm before BKZ, Norm after BKZ
    _prev_ln    = RealField(work_prec)(1);
    _prevBKZ_ln = RealField(work_prec)(1);

    for _k in range(GSB.nrows()):
        _gsl_ln = log(GSL[_k].norm());
        _gsb_ln = log(GSB[_k].norm());
        _gsl_rat = exp(_gsl_ln - _prev_ln);    assert(_gsl_rat.parent().precision() >= work_prec);
        _gsb_rat = exp(_gsb_ln - _prevBKZ_ln); assert(_gsb_rat.parent().precision() >= work_prec);
        # prepare next ratio
        _prev_ln = _gsl_ln;
        _prevBKZ_ln = _gsb_ln;

        # Output
        data_line = [_k, _gsl_ln, _gsb_ln, _gsl_rat, _gsb_rat];
        print_data_line(f_out, data_line);

    print ("    out in '{}'".format(out_name));
    f_out.close();

exit;
