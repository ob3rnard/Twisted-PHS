#!/usr/bin/env sage

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------
from sage.all import *
from number_field import *

bsz=100;
N=50;

# label = ARG1
#typs  = ["tw","opt0","phs0","opt","phs"]

data_dir = "./data/"

# Data files
cyclo_fields = [23,29,31,37,41,43,47,53];
ntru_fields  = [23,29,31,37];


def read_data(filename): # read by columns
    f_in = open(filename, "r");
    f_in.readline(); f_in.readline();

    _exact = [];
    _tw    = [];
    _opt   = [];
    _phs   = [];
    for _line in f_in:
        if (_line[0] == '#'):
            continue;
        _line.strip();
        _values = [sage_eval(_v) for _v in _line[:-1].split("\t")]; assert(len(_values) == 6);
        _exact += [_values[0]];
        _tw    += [_values[1]];
        _opt   += [_values[4]];
        _phs   += [_values[5]];

    assert (len(_exact) == len(_tw)) and (len(_exact) == len(_opt)) and (len(_exact) == len(_phs));
    return _exact, _tw, _opt, _phs;
    

# g = Graphics();
_moy_file = data_dir + "cyclo.solmoy"
_moy = open(_moy_file, "w");
_rat_file = data_dir + "cyclo.solrat";
_rat = open(_rat_file, "w");
for m in cyclo_fields:
    tag = "z{}".format(m);
    K = nf_set_tag(tag);
    sol_file = data_dir + tag + "/{}.sol_b{}_n{}".format(tag, bsz, N);
    _exact, _tw, _opt, _phs = read_data(sol_file);

    _tw_rat  = [ RealField(1000)(_tw[_i]/_exact[_i]) for _i in range(len(_exact)) ];
    _opt_rat = [ RealField(1000)(_opt[_i]/_exact[_i]) for _i in range(len(_exact)) ];
    _phs_rat  = [ RealField(1000)(_phs[_i]/_exact[_i]) for _i in range(len(_exact)) ];

    _tw_moy  = mean(_tw_rat);
    _opt_moy = mean(_opt_rat);
    _phs_moy = mean([_f for _f in _phs_rat]);

    # All ratios
    for _k in range(len(_exact)):
        _rat.write("{} {:7.3f} {:7.3f} {:7.3f}\n".format(K.degree(),float(_tw_rat[_k]),float(_opt_rat[_k]),float(_phs_rat[_k])));
    # Mean
    _moy.write("{} {} {} {}\n".format(K.degree(), *[RealField(100)(_f) for _f in [_tw_moy,_opt_moy,_phs_moy]]));

_moy.close();
_rat.close();


# ------------------------------
# Same for NTRU
_moy_file = data_dir + "ntru.solmoy";
_moy = open(_moy_file, "w");
_rat_file = data_dir + "ntru.solrat";
_rat = open(_rat_file, "w");
for m in ntru_fields:
    tag = "n{}".format(m);
    K = nf_set_tag(tag);
    sol_file = data_dir + tag + "/{}.sol_b{}_n{}".format(tag, bsz, N);
    _exact, _tw, _opt, _phs = read_data(sol_file);

    _tw_rat  = [ RealField(1000)(_tw[_i]/_exact[_i]) for _i in range(len(_exact)) ];
    _opt_rat = [ RealField(1000)(_opt[_i]/_exact[_i]) for _i in range(len(_exact)) ];
    _phs_rat  = [ RealField(1000)(_phs[_i]/_exact[_i]) for _i in range(len(_exact)) ];

    _tw_moy  = mean(_tw_rat);
    _opt_moy = mean(_opt_rat);
    _phs_moy = mean([_f for _f in _phs_rat]);

    # All ratios
    for _k in range(len(_exact)):
        _rat.write("{} {:7.3f} {:7.3f} {:7.3f}\n".format(K.degree(),float(_tw_rat[_k]),float(_opt_rat[_k]),float(_phs_rat[_k])));
    # Mean
    _moy.write("{} {} {} {}\n".format(K.degree(), *[RealField(100)(_f) for _f in [_tw_moy,_opt_moy,_phs_moy]]));

_moy.close();
_rat.close();

exit;
