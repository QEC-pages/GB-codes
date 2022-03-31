# GB-codes

This repository contains the quantum error-correcting codes constructed 
by Renyu Wang and Leonid P. Pryadko for the paper "Distance bounds for generalized bicycle codes".

## Files
- The zip arxive `gb-codes.zip` contains the actual code files in the MTX format.
- The GAP program `gb-verify.gap` contains the list of all codes (circulant length, distance, and two lists of degrees of the non-zero monomials) and a GAP program to verify the distance or write the files.
- The file `gb-verify_gap.out.gz` is the stdout generated by `gap.sh gb-verify.gap`.  


## Parameters of the codes
The following is the output of `zgrep dZ= gb-verify_gap.out.gz`
```
 n0=5   [ 0, 2 ] [0,1]   n=10 dX=3 dZ=3 d0=3
 n0=11  [ 0, 3 ] [0,1]   n=22 dX=4 dZ=5 d0=4
 n0=13  [ 0, 5 ] [0,1]   n=26 dX=6 dZ=5 d0=5
 n0=19  [ 0, 4 ] [0,1]   n=38 dX=5 dZ=6 d0=5
 n0=29  [ 0, 8 ] [0,1]   n=58 dX=9 dZ=7 d0=7
 n0=37  [ 0, 8 ] [0,1]   n=74 dX=9 dZ=8 d0=8
 n0=53  [ 0, 8 ] [0,1]   n=106 dX=9 dZ=10 d0=9
 n0=59  [ 0, 9 ] [0,1]   n=118 dX=10 dZ=11 d0=10
 n0=61  [ 0, 11 ] [0,1]  n=122 dX=12 dZ=11 d0=11
 n0=67  [ 0, 12 ] [0,1]  n=134 dX=13 dZ=11 d0=11
 n0=83  [ 0, 11 ] [0,1]  n=166 dX=12 dZ=13 d0=12
 n0=101         [ 0, 12 ] [0,1]  n=202 dX=13 dZ=13 d0=13
 n0=107         [ 0, 41 ] [0,1]  n=214 dX=14 dZ=15 d0=14
 n0=131         [ 0, 20 ] [0,1]  n=262 dX=15 dZ=16 d0=15
 n0=139         [ 0, 30 ] [0,1]  n=278 dX=17 dZ=16 d0=16
 n0=149         [ 0, 34 ] [0,1]  n=298 dX=17 dZ=17 d0=17
 n0=163         [ 0, 17 ] [0,1]  n=326 dX=18 dZ=17 d0=17
 n0=173         [ 0, 51 ] [0,1]  n=346 dX=18 dZ=19 d0=18
 n0=179         [ 0, 17 ] [0,1]  n=358 dX=18 dZ=19 d0=18
 n0=181         [ 0, 19 ] [0,1]  n=362 dX=20 dZ=19 d0=19
 n0=197         [ 0, 23 ] [0,1]  n=394 dX=20 dZ=19 d0=19
 n0=211         [ 0, 20 ] [0,1]  n=422 dX=21 dZ=20 d0=20
 n0=227         [ 0, 61 ] [0,1]  n=454 dX=22 dZ=21 d0=21
 
 n0=5   [ 0, 1, 2, 3 ] [0,1]     n=10 dX=3 dZ=3 d0=3
 n0=11  [ 0, 1, 2, 5 ] [0,1]     n=22 dX=5 dZ=5 d0=5
 n0=13  [ 0, 1, 2, 5 ] [0,1]     n=26 dX=5 dZ=6 d0=5
 n0=19  [ 0, 2, 3, 7 ] [0,1]     n=38 dX=7 dZ=7 d0=7
 n0=29  [ 0, 1, 3, 10 ] [0,1]    n=58 dX=9 dZ=9 d0=9
 n0=37  [ 0, 1, 4, 14 ] [0,1]    n=74 dX=12 dZ=11 d0=11
 n0=53  [ 0, 3, 4, 15 ] [0,1]    n=106 dX=15 dZ=13 d0=13
 n0=59  [ 0, 8, 11, 21 ] [0,1]   n=118 dX=15 dZ=15 d0=15
 n0=61  [ 0, 6, 9, 22 ] [0,1]    n=122 dX=16 dZ=15 d0=15
 n0=67  [ 0, 3, 24, 41 ] [0,1]   n=134 dX=17 dZ=16 d0=16
 n0=83  [ 0, 5, 13, 27 ] [0,1]   n=166 dX=18 dZ=19 d0=18
 n0=101         [ 0, 9, 15, 56 ] [0,1]   n=202 dX=21 dZ=21 d0=21
 n0=107         [ 0, 11, 14, 29 ] [0,1]  n=214 dX=21 dZ=21 d0=21
 n0=131         [ 0, 9, 17, 32 ] [0,1]   n=262 dX=23 dZ=23 d0=23
 n0=139         [ 0, 10, 39, 75 ] [0,1]  n=278 dX=25 dZ=25 d0=25
 n0=149         [ 0, 15, 20, 39 ] [0,1]  n=298 dX=25 dZ=25 d0=25
 n0=163         [ 0, 8, 25, 44 ] [0,1]   n=326 dX=28 dZ=27 d0=27
 n0=173         [ 0, 13, 22, 46 ] [0,1]  n=346 dX=28 dZ=29 d0=28
 n0=179         [ 0, 14, 51, 99 ] [0,1]  n=358 dX=29 dZ=29 d0=29
 n0=181         [ 0, 12, 52, 94 ] [0,1]  n=362 dX=29 dZ=29 d0=29
 n0=197         [ 0, 17, 30, 54 ] [0,1]  n=394 dX=30 dZ=31 d0=30
 n0=211         [ 0, 12, 34, 54 ] [0,1]  n=422 dX=33 dZ=31 d0=31
 n0=227         [ 0, 22, 27, 63 ] [0,1]  n=454 dX=35 dZ=33 d0=33
 
 n0=5   [ 0, 1, 2, 3 ] [0,1]     n=10 dX=3 dZ=3 d0=3
 n0=7   [ 0, 1, 2, 3 ] [0,1]     n=14 dX=3 dZ=4 d0=3
 n0=11  [ 0, 1, 2, 5 ] [0,1]     n=22 dX=5 dZ=5 d0=5
 n0=13  [ 0, 1, 2, 5 ] [0,1]     n=26 dX=5 dZ=6 d0=5
 n0=17  [ 0, 1, 3, 9 ] [0,1]     n=34 dX=8 dZ=7 d0=7
 n0=19  [ 0, 1, 2, 8 ] [0,1]     n=38 dX=8 dZ=7 d0=7
 n0=23  [ 0, 1, 3, 9 ] [0,1]     n=46 dX=8 dZ=9 d0=8
 n0=29  [ 0, 1, 3, 10 ] [0,1]    n=58 dX=9 dZ=9 d0=9
 n0=31  [ 0, 5, 11, 17 ] [0,1]   n=62 dX=10 dZ=11 d0=10
 n0=37  [ 0, 1, 4, 14 ] [0,1]    n=74 dX=12 dZ=11 d0=11
 n0=41  [ 0, 4, 15, 27 ] [0,1]   n=82 dX=13 dZ=12 d0=12
 n0=43  [ 0, 3, 8, 17 ] [0,1]    n=86 dX=13 dZ=12 d0=12
 n0=47  [ 0, 2, 5, 18 ] [0,1]    n=94 dX=14 dZ=13 d0=13
 n0=53  [ 0, 2, 5, 17 ] [0,1]    n=106 dX=13 dZ=13 d0=13
 n0=59  [ 0, 5, 22, 31 ] [0,1]   n=118 dX=15 dZ=15 d0=15
 n0=61  [ 0, 5, 13, 25 ] [0,1]   n=122 dX=16 dZ=15 d0=15
 n0=67  [ 0, 3, 24, 41 ] [0,1]   n=134 dX=17 dZ=16 d0=16
 n0=71  [ 0, 4, 11, 23 ] [0,1]   n=142 dX=17 dZ=17 d0=17
 n0=73  [ 0, 8, 20, 46 ] [0,1]   n=146 dX=17 dZ=17 d0=17
 n0=79  [ 0, 5, 11, 24 ] [0,1]   n=158 dX=17 dZ=17 d0=17
 n0=83  [ 0, 5, 13, 27 ] [0,1]   n=166 dX=18 dZ=19 d0=18
 n0=89  [ 0, 7, 16, 27 ] [0,1]   n=178 dX=19 dZ=19 d0=19
 n0=97  [ 0, 7, 15, 54 ] [0,1]   n=194 dX=21 dZ=20 d0=20
 n0=101         [ 0, 6, 15, 60 ] [0,1]   n=202 dX=22 dZ=21 d0=21
 n0=103         [ 0, 5, 24, 46 ] [0,1]   n=206 dX=22 dZ=21 d0=21
 n0=107         [ 0, 6, 16, 35 ] [0,1]   n=214 dX=22 dZ=21 d0=21
 n0=109         [ 0, 4, 15, 31 ] [0,1]   n=218 dX=21 dZ=21 d0=21
 n0=113         [ 0, 7, 18, 31 ] [0,1]   n=226 dX=21 dZ=21 d0=21
 
 n0=11  [ 0, 1, 2, 4, 5, 8 ] [0,1]       n=22 dX=7 dZ=6 d0=6
 n0=13  [ 0, 2, 3, 6, 7, 9 ] [0,1]       n=26 dX=6 dZ=7 d0=6
 n0=19  [ 0, 1, 2, 3, 4, 8 ] [0,1]       n=38 dX=7 dZ=7 d0=7
 n0=29  [ 0, 1, 2, 7, 10, 13 ] [0,1]     n=58 dX=10 dZ=11 d0=10
 n0=37  [ 0, 3, 4, 9, 13, 16 ] [0,1]     n=74 dX=12 dZ=13 d0=12
 n0=53  [ 0, 2, 6, 12, 13, 19 ] [0,1]    n=106 dX=15 dZ=15 d0=15
 n0=59  [ 0, 6, 7, 10, 23, 31 ] [0,1]    n=118 dX=18 dZ=17 d0=17
 n0=67  [ 0, 4, 10, 16, 24, 33 ] [0,1]   n=134 dX=18 dZ=19 d0=18
 n0=83  [ 0, 7, 10, 16, 18, 35 ] [0,1]   n=166 dX=21 dZ=21 d0=21
 n0=101         [ 0, 1, 4, 16, 27, 36 ] [0,1]    n=202 dX=23 dZ=23 d0=23
 n0=107         [ 0, 9, 13, 21, 27, 59 ] [0,1]   n=214 dX=26 dZ=25 d0=25
```

