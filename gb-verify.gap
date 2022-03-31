#!
#! This GAP program verifies parameters of the GB codes constructed for the paper 
#! "Distance bounds for generalized bicycle code" 
#! by Renyu Wang and Leonid P. Pryadko
#! 

# adjust the path in the second argument if needed
SetPackagePath("QDistRnd","./QDistRnd" );
# package available from https://github.com/QEC-pages/QDistRnd/
LoadPackage("QDistRnd"); 

# In the following, each row has the structure: 
# [n0,d,[a0,a1,a2,...],[b0,b1,...]]
# where n0 is the circulant size, d0 is the GB code distance 
# components of the two vectors are degrees of binary polynomials 
# a(x)=x^a0+x^a1+... and b(x).

#! weight-four codes with a(x)=x^a0+x^a1, b(x)=1+x of length n=2*n0
parlst2:=[ [ 5,  3,  [ 0, 2 ],  [ 0, 1 ]],
           [ 11, 4,  [ 0, 3 ],  [ 0, 1 ]],
           [ 13, 5,  [ 0, 5 ],  [ 0, 1 ]],
           [ 19, 5,  [ 0, 4 ],  [ 0, 1 ]],
           [ 29, 7,  [ 0, 8 ],  [ 0, 1 ]],
           [ 37, 8,  [ 0, 8 ],  [ 0, 1 ]],
           [ 53, 9,  [ 0, 8 ],  [ 0, 1 ]],
           [ 59, 10, [ 0, 9 ],  [ 0, 1 ]],
           [ 61, 11, [ 0, 11 ], [ 0, 1 ]],
           [ 67, 11, [ 0, 12 ], [ 0, 1 ]],
           [ 83, 12, [ 0, 11 ], [ 0, 1 ]],
           [101, 13, [ 0, 12 ], [ 0, 1 ]],
           [107, 14, [ 0, 41 ], [ 0, 1 ]],
           [131, 15, [ 0, 20 ], [ 0, 1 ]],
           [139, 16, [ 0, 30 ], [ 0, 1 ]],
           [149, 17, [ 0, 34 ], [ 0, 1 ]],
           [163, 17, [ 0, 17 ], [ 0, 1 ]],
           [173, 18, [ 0, 51 ], [ 0, 1 ]],
           [179, 18, [ 0, 17 ], [ 0, 1 ]],
           [181, 19, [ 0, 19 ], [ 0, 1 ]],
           [197, 19, [ 0, 23 ], [ 0, 1 ]],
           [211, 20, [ 0, 20 ], [ 0, 1 ]],
           [227, 21, [ 0, 61 ], [ 0, 1 ]]
         ];

# List of GB codes of w=6 
# for circulant sizes n0 with primitive root 2  
parlst4:=[ [ 5,    3,    [ 0, 1, 2, 3 ],  [ 0, 1 ]],
           [ 11,   5,    [ 0, 1, 2, 5 ],  [ 0, 1 ]],
           [ 13,   5,    [ 0, 1, 2, 5 ],  [ 0, 1 ]],
           [ 19,   7,    [ 0, 2, 3, 7 ],  [ 0, 1 ]],
           [ 29,   9,    [ 0, 1, 3, 10 ], [ 0, 1 ]],
           [ 37,   11,   [ 0, 1, 4, 14 ], [ 0, 1 ]],
           [ 53,   13,   [ 0, 3, 4, 15 ], [ 0, 1 ]],
           [ 59,   15,   [ 0, 8, 11, 21 ], [ 0, 1 ]],
           [ 61,   15,   [ 0, 6, 9, 22 ], [ 0, 1 ]],
           [ 67,   16,   [ 0, 3, 24, 41 ], [ 0, 1 ]],
           [ 83,   18,   [ 0, 5, 13, 27 ], [ 0, 1 ]],
           [ 101,  21,   [ 0, 9, 15, 56 ],        [ 0, 1 ]],
           [ 107,  21,   [ 0, 11, 14, 29 ],       [ 0, 1 ]],
           [ 131,  23,   [ 0, 9, 17, 32 ],        [ 0, 1 ]], 
           [ 139,  25,   [ 0, 10, 39, 75 ],       [ 0, 1 ]], 
           [ 149, 25,    [ 0, 15, 20, 39 ],       [ 0, 1 ]], 
           [163, 27,     [ 0, 8, 25, 44 ],        [ 0, 1 ]], 
           [173, 28,     [ 0, 13, 22, 46 ],       [ 0, 1 ]], 
           [179, 29,     [ 0, 14, 51, 99 ],       [ 0, 1 ]], 
           [181, 29,     [ 0, 12, 52, 94 ],       [ 0, 1 ]], 
#           [ 187, 45,    [ 0, 12, 20, 58 ],       [0,1]], # incorrect distance
           [197, 30,     [ 0, 17, 30, 54 ],       [ 0, 1 ]], 
           [211, 31,     [ 0, 12, 34, 54 ],       [0, 1]],   
           [227, 33,     [ 0, 22, 27, 63 ],       [ 0, 1 ]]
         ];

## prime circulant size , w=6
parlst4B:=[ [5, 3,  [ 0, 1, 2, 3 ], [ 0, 1 ]],  
            [7, 3,  [ 0, 1, 2, 3 ], [ 0, 1 ]],  
            [11, 5,  [ 0, 1, 2, 5 ], [ 0, 1 ]], 
            [13, 5,  [ 0, 1, 2, 5 ], [ 0, 1 ]], 
            [17, 7,  [ 0, 1, 3, 9 ], [ 0, 1 ]], 
            [19, 7,  [ 0, 1, 2, 8 ], [ 0, 1 ]], 
            [23, 8,  [ 0, 1, 3, 9 ], [ 0, 1 ]], 
            [29, 9,  [ 0, 1, 3, 10 ], [ 0, 1 ]],
            [31, 10,  [ 0, 5, 11, 17 ], [ 0, 1 ]],
            [37, 11,  [ 0, 1, 4, 14 ], [ 0, 1 ]], 
            [41, 12,  [ 0, 4, 15, 27 ], [ 0, 1 ]],
            [43, 12,  [ 0, 3, 8, 17 ], [ 0, 1 ]], 
            [47, 13,  [ 0, 2, 5, 18 ], [ 0, 1 ]], 
            [53, 13,  [ 0, 2, 5, 17 ], [ 0, 1 ]], 
            [59, 15,  [ 0, 5, 22, 31 ], [ 0, 1 ]],
            [61, 15,  [ 0, 5, 13, 25 ], [ 0, 1 ]],
            [67, 16,  [ 0, 3, 24, 41 ], [ 0, 1 ]],
            [71, 17,  [ 0, 4, 11, 23 ], [ 0, 1 ]],
            [73, 17,  [ 0, 8, 20, 46 ], [ 0, 1 ]],
            [79, 17,  [ 0, 5, 11, 24 ], [ 0, 1 ]],
            [83, 18,  [ 0, 5, 13, 27 ], [ 0, 1 ]],
            [89, 19,  [ 0, 7, 16, 27 ], [ 0, 1 ]],
            [97, 20,  [ 0, 7, 15, 54 ], [ 0, 1 ]],
            [101, 21,  [ 0, 6, 15, 60 ], [ 0, 1 ]],
            [103, 21,  [ 0, 5, 24, 46 ], [ 0, 1 ]],
            [107, 21,  [ 0, 6, 16, 35 ], [ 0, 1 ]],
            [109, 21,  [ 0, 4, 15, 31 ], [ 0, 1 ]],
            [113, 21,  [ 0, 7, 18, 31 ], [ 0, 1 ]]
#            [127, 23,  [ 0, 3, 66, 86 ], [ 0, 1 ]] # incorrect distance
          ];

# this is bicycle w=8
parlst8:=[ [11, 6,     [ 0, 1, 2, 4, 5, 8 ],     [ 0, 1 ]],
           [13, 6,     [ 0, 2, 3, 6, 7, 9 ],     [ 0, 1 ]],
           [19, 7,     [ 0, 1, 2, 3, 4, 8 ],     [ 0, 1 ]],
           [29, 10,    [ 0, 1, 2, 7, 10, 13 ],   [ 0, 1 ]],
           [37, 12,    [ 0, 3, 4, 9, 13, 16 ],   [ 0, 1 ]],
           [53, 15,    [ 0, 2, 6, 12, 13, 19 ],  [ 0, 1 ]],
           [59, 17,    [ 0, 6, 7, 10, 23, 31 ],  [ 0, 1 ]],
#          [61, 17,    [ 0, 5, 10, 13, 14, 21 ], [ 0, 1 ]], # incorrect distance
           [67, 18,    [ 0, 4, 10, 16, 24, 33 ], [ 0, 1 ]], 
           [83, 21,    [ 0, 7, 10, 16, 18, 35 ], [ 0, 1 ]],   
           [101, 23,   [0, 1, 4, 16, 27, 36],    [ 0, 1 ]],
           [107, 25,   [ 0, 9, 13, 21, 27, 59 ], [ 0, 1 ]]
         ];

#131 28 ???
#139 
#149
#163
#173 33 #??? not converged 
#179
#181
#197
#211
#227

# construct m by n circulant matrix out of the coefficients of the polynomial,
# working with the field F
DoCirc:=function(poly,m,n,F)
    local v,perm,j,deg,mat;
    v:=CoefficientsOfUnivariatePolynomial(poly);
#    F:=DefaultField(v);
    deg:=Length(v);
    if (n>deg) then
        v:=Concatenation(v,ListWithIdenticalEntries(n-deg,0*One(F)));
    fi;
    mat:=[v];
    perm:=Concatenation([n],[1..n-1]);
    for j in [2..m] do
        v:=v{perm};
        Append(mat,[v]);
    od;
    return mat;
end;


parlst:=Concatenation(parlst2,parlst4,parlst4B,parlst8);


q:=2; F:=GF(q);  # field to be used
x:=Indeterminate(F,"x");;

docalc:=false;
if docalc then 
    for row in parlst do
        n0:=row[1];
        d0:=row[2];           
        v1:=row[3]; h1:=0;    
        for j in v1 do
            h1:=h1+One(F)*x^j;
        od;
        v2:=row[4]; h2:=0;    
        for j in v2 do
            h2:=h2+One(F)*x^j;
        od;
        m1:=DoCirc(h1,n0,n0,F);
        m2:=DoCirc(h2,n0,n0,F);
        G:=TransposedMatMutable(Concatenation(m2,m1));
        hh:=Gcd(h1,h2);
        h1x:=h1/hh;
        h2x:=h2/hh;
        #    Print(h1," ",h2," ",h1x, " ", h2x,"\n");
        
        m1:=QDR_DoCirc(h1x,n0,n0,F);
        m2:=QDR_DoCirc(h2x,n0,n0,F);
        H:=TransposedMatMutable(Concatenation(TransposedMat(m1),TransposedMat(-m2)));
        dZ:=DistRandCSS(G,H,10^5,d0-1,12 : field:=F,maxav:=500/n0);
        dX:=DistRandCSS(H,G,10^5,d0-1,12 : field:=F,maxav:=500/n0);
        if (dX<d0) or (dZ<d0) 
        then str:="##################";  
        else str:="";   
        fi;
            
        Print(" n0=",n0, " \t",v1," [0,1]\t n=",2*n0," dX=",dX," dZ=",dZ," d0=",d0," ",str,"\n");
        #    if (n0>20) then break; fi;    
    od;
else # just generate matrix files in ./codes/    
        for row in parlst do
        n0:=row[1];
        d0:=row[2];           
        v1:=row[3]; h1:=0;    
        for j in v1 do
            h1:=h1+One(F)*x^j;
        od;
        v2:=row[4]; h2:=0;    
        for j in v2 do
            h2:=h2+One(F)*x^j;
        od;
        m1:=DoCirc(h1,n0,n0,F);
        m2:=DoCirc(h2,n0,n0,F);
        G:=TransposedMatMutable(Concatenation(m2,m1));
        H:=TransposedMatMutable(Concatenation(TransposedMat(m1),TransposedMat(-m2)));
        w:=Length(v1)+Length(v2);        
        fname:=Concatenation("./codes/GB_",String(2*n0),"_w",String(w),"_X.mtx");
        comment:=Concatenation("GB code [[",String(2*n0),",2,",String(d0),
                               "]] w=",String(w)); 
        WriteMTXE(fname, 0,G, comment : field := F );           
        fname:=Concatenation("./codes/GB_",String(2*n0),"_w",String(w),"_Z.mtx");
        WriteMTXE(fname, 0,H, comment : field := F );
#            if (n0>20) then break; fi;    
                 
    od;

fi;
