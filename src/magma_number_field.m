SetClassGroupBounds("GRH");

/* -------------------------------------------------------------------------------------------------- */
/* ----- S-UNITS COMPUTATIONS ----- */

/* ----- S-Units wrpt. some FB ----- */
// Returns rk_u + #FB fundamental S-units wrpt. FB.
function get_S_units(K, FB)
    _r1, _r2   := Signature(K);
    _G_SU, _mU := SUnitGroup(FB);
    _SU        := [K| _mU(_G_SU.i): i in [1..NumberOfGenerators(_G_SU)]
                   | AbsoluteValue(Degree(K)-Length(_mU(_G_SU.i))) gt 2^(-100)];
    assert (#_SU eq (#FB + _r1 + _r2 - 1));
    return _SU;
end function;


/* ----- ClDL ----- */
/* WARN: This function assumes that I is coprime with FB */
/* Note: the S-unit are already huge, so scrambling is quickly not a good idea.
         If scramble := false, option S_fb is useless, otherwise you can give precomputed S-units wrpt. FB. */
function get_id_clDL(K, I, FB: scramble:=false, S_fb:=[])
    // For convenience, just assume I is prime and not in FB
    assert(IsPrime(I) and (I notin FB));
    r1, r2 := Signature(K);
    _I_fact := [ I ];
    _FBI := FB cat _I_fact;

    // Compute S-units for this extended factor base
    t := Cputime(); _SS := get_S_units(K, _FBI); t := Cputime(t);
    // printf "ClDL: S-units(fb+a)=%os\n", t;
    _SS:=_SS[r1+r2..#_SS]; assert(#_SS eq (#FB+1));

    // Construct clDL linear system from the relation matrix
    _sys_clDL:= Matrix([[ Valuation(_s, _pid): _pid in _I_fact ]: _s in _SS]);
    _vI      := Vector([Valuation(I, _pid): _pid in _I_fact]);
    
    // Solution Space
    _clDL, _clDL_ker := Solution(_sys_clDL, _vI);
    _clDL := Vector(_clDL); _clDL_ker := BasisMatrix(_clDL_ker); // (sic)

    /* Scrambling (random elt of kernel + random S-unit for FB) */
    if (scramble eq true) then
        if (#S_fb eq 0) then
            t := Cputime(); S_fb := get_S_units(K, FB); t := Cputime(t); // Include units but it's ok
            //printf "->scramble: S-units(fb)=%os\n", t;
        end if;
        _clDL +:= Vector([ Random([-2..2]): _k in [1..NumberOfRows(_clDL_ker)] ])*_clDL_ker;
        _sun   := PowerProduct(S_fb, [ Random([-2..2]): _k in [1..#S_fb] ]);
    else
        _sun := K!1;
    end if;

    /* Construct the answer eta */
    _eta := _sun * PowerProduct(_SS, Eltseq(_clDL));
    
    /* Trust does not exclude control, but:
       /!\ Asserting the result is principal and generator only differs by a unit is (very) costly */
    assert (_clDL*_sys_clDL eq _vI);
    //_P   := PowerProduct(_FBI, Eltseq(_clDL*_rel_mat)[1..#_FBI-#_I_fact] cat [0: _k in [1..#_I_fact]]);
    //_b,_w:= IsPrincipal(I*_P);
    // assert (_b eq true and Norm(_eta/_w) eq 1);
    
    return _eta;
end function;


/* -------------------------------------------------------------------------------------------------- */
/* INPUT/OUTPUT 
   Format SU, Chal, FB, NF :--> see Sage "number_field.py" */


/* ----------------------------------------- */
/* Number Field */
// in: 'z23' or 'n23' to resp. Cyclo(23) and NTRU(23)
function nf_set_tag(tag)
    Zx<x>:=PolynomialRing(Integers());
    typ  := tag[1];
    m    := StringToInteger(tag[2..#tag]);
    
    if (typ eq "z") then
        _K := CyclotomicField(m);
    elif (typ eq "n") then
        _K := NumberField(x^m-x-1);
    else
        printf "**Warn** nf tag %o unknown\n", tag;
    end if;
        
    return _K;
end function;

function __nf_out_stream(K, tag)
    nf_str := Sprintf("# number_field: tag=%o eq=%o", tag,
                      &cat[ss: ss in Split(Sprintf("%o",Eltseq(DefiningPolynomial(K))), " ")]);
    return nf_str;
end function;


/* ----------------------------------------- */
/* Prime ideals: for FB and Challenges */
function __fb_pid_in_str(line, K)
    g2 := Split(line, " \n"); assert(#g2 eq 2);
    p  := StringToInteger(g2[1]);
    a  := Evaluate(Polynomial([StringToInteger(_a): _a in Split(g2[2], "[,]")]), K.1);
    /* Didn't find a way not to compute OK.
       This is a problem for NTRU Prime fields, where the equation order IS generically maximal.
       We would like to tell Magma: "Please suppose it is maximal and don't bother". */
    pid := ideal< MaximalOrder(K)| p, a>;

    return pid;
end function;


function fb_in_stream(filename, K)
    F := Open(filename, "r");
    // ignore first line
    _ := Gets(F);
    // read second line to find k
    fb_head := Gets(F); assert(not IsEof(fb_head));
    exp     := "# fb_places: k=";
    assert(fb_head[1..#exp] eq exp);
    k       := StringToInteger(fb_head[#exp+1..#fb_head]);
    
    // read each line and convert
    FB := [];
    for _i in [1..k] do
        _line:= Gets(F); assert(not IsEof(_line));
        _pid := __fb_pid_in_str(_line, K);
        Append(~FB, _pid);
    end for;
    
    delete F;
    return FB;
end function;


/* ----------------------------------------- */
/* Sunits */

// NB: returns a string
function __su_out_stream(K, su)
    _su_str := &cat[_su_c: _su_c in Split(Sprintf("%o",Eltseq(Polynomial(Eltseq(su))))," ")];
    return _su_str;
end function;

procedure sunits_out_stream(filename, K, tag, S) // + Nf Head !
    F := Open(filename, "w");
    r1, r2 := Signature(K);
    nf_str := __nf_out_stream(K, tag);

    fprintf F, nf_str cat "\n";
    fprintf F, "# sunits(fb): nu=%o k=%o\n", r1+r2-1, #S-(r1+r2-1);
    for _i in [1..#S] do
        su_str := __su_out_stream(K, S[_i]);
        fprintf F, su_str cat "\n";
    end for;
    fprintf F, "# --- END ---\n";
    
    delete F;
end procedure;


/* ----------------------------------------- */
/* Challenges: in=prime ideal, out=element of K */
// out
function __cldl_out_stream(K, eta)
    return __su_out_stream(K, eta);
end function;

// in
function __chal_in_str(line, K)
    return __fb_pid_in_str(line, K);
end function;

function chal_in_stream(filename, K)
    F := Open(filename, "r");
    // ignore first line, which is "# nf:{} eq:[a,b,c,...]" for K = Q[x]/<x^d + ... + cx^2+bx+1.
    _ := Gets(F);
    // read second line to find k
    chal_head := Gets(F); assert(not IsEof(chal_head));
    exp := "# chals: k=";
    assert(chal_head[1..#exp] eq exp);
    kb := Split(chal_head[#exp+1..#chal_head], " "); assert(#kb eq 2);
    k  := StringToInteger(kb[1]);
    bsz:= StringToInteger(Split(kb[2],"=")[2]);
    
    // read each line and convert
    chals := [];
    for _i in [1..k] do
        _line:= Gets(F); assert(not IsEof(_line));
        _pid := __chal_in_str(_line, K);
        Append(~chals, _pid);
    end for;
    
    delete F;
    return chals, bsz;
end function;

