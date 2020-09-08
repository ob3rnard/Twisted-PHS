# Floating Point arithmetic handling
from sage.all import *


# Default precision (arbitrary)
BIT_PREC_DEFAULT = 300;

# -------------------------------------------------------------------------------------------
# Test at zero
# 1. It is almost impossible to guarantee when working at precision N that the precision drift
#    will be contained at +20 (even if in most cases it is the case, some functions are not
#    numerically stable, e.g. matrix inversion is particularly difficult to this respect.
# 2. Hence, it is a problem to decide which is an acceptable drift, and what is not.
#    We choose to accept values under lN + O, and suggest l=__BIT_PREC_SLOPPY and O = __BIT_PREC_SHIFT for sloppy results.
#    The offset/shift is always tolerated, even in the "exact" case.
# 3. The provided interface allows to relax the requirements locally.
# 4. There is a non-modifiable maximum gap: if the result is True, it ensures val < 2^{__BIT_PREC_MAX_ACCEPT}.

__BIT_PREC_EXACT      = 1.0;
__BIT_PREC_SLOPPY     = 0.25; # Seems impossible to guarantee M*M^(-1) = 1 at prec(M)
__BIT_PREC_SHIFT      =  30; # 4 or 5 digits are generally sufficient (somewhere btw 10^4 and 10^5).
__BIT_PREC_MAX_ACCEPT = -32;


# Checks against zero
# If sloppy = False (default), __sloppyness is automatically set to __BIT_PREC_EXACT.
# Csqc: if you want to use another sloppyness factor, you must pass sloppy=True.
def fp_check_zero(description, lval, target=BIT_PREC_DEFAULT, sloppy=False,
                  __sloppyness=__BIT_PREC_SLOPPY, __shift=__BIT_PREC_SHIFT):
    __sloppyness = __BIT_PREC_EXACT if (sloppy == False) else __sloppyness;
    _threshold   = min(__BIT_PREC_MAX_ACCEPT, -target*__sloppyness+__shift);
    _l2_infty    = -Infinity if (len(lval) == 0) else max([log(abs(_val), 2) for _val in lval]);
    _is_zero     = True if (_l2_infty < _threshold) else False;

    # Comprehensive error message
    if (_is_zero == False):
        print ("[Err] log_2 infty-norm of '{}' is {:.2f} [thr:{:.2f}/exp:{:.2f}]".format(description,float(_l2_infty),float(_threshold),float(-target)));
    
    return _is_zero;
