$offlisting
*  
*  Equation counts
*      Total        E        G        L        N        X        C        B
*          1        1        0        0        0        0        0        0
*  
*  Variable counts
*                   x        b        i      s1s      s2s       sc       si
*      Total     cont   binary  integer     sos1     sos2    scont     sint
*        121        1      120        0        0        0        0        0
*  FX      0        0        0        0        0        0        0        0
*  
*  Nonzero counts
*      Total    const       NL      DLL
*        121        1      120        0
*
*  Solve m using MIQCP maximizing objvar;


Variables  objvar,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18
          ,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35
          ,b36,b37,b38,b39,b40,b41,b42,b43,b44,b45,b46,b47,b48,b49,b50,b51,b52
          ,b53,b54,b55,b56,b57,b58,b59,b60,b61,b62,b63,b64,b65,b66,b67,b68,b69
          ,b70,b71,b72,b73,b74,b75,b76,b77,b78,b79,b80,b81,b82,b83,b84,b85,b86
          ,b87,b88,b89,b90,b91,b92,b93,b94,b95,b96,b97,b98,b99,b100,b101,b102
          ,b103,b104,b105,b106,b107,b108,b109,b110,b111,b112,b113,b114,b115
          ,b116,b117,b118,b119,b120,b121;

Binary Variables  b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18
          ,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35
          ,b36,b37,b38,b39,b40,b41,b42,b43,b44,b45,b46,b47,b48,b49,b50,b51,b52
          ,b53,b54,b55,b56,b57,b58,b59,b60,b61,b62,b63,b64,b65,b66,b67,b68,b69
          ,b70,b71,b72,b73,b74,b75,b76,b77,b78,b79,b80,b81,b82,b83,b84,b85,b86
          ,b87,b88,b89,b90,b91,b92,b93,b94,b95,b96,b97,b98,b99,b100,b101,b102
          ,b103,b104,b105,b106,b107,b108,b109,b110,b111,b112,b113,b114,b115
          ,b116,b117,b118,b119,b120,b121;

Equations  e1;


e1.. 48*b2*b3 + 100*b2 + 6*b3 - 2*b2*b4 + 71*b4 - 84*b2*b8 - 55*b8 + 88*b2*b9
      + 100*b9 + 12*b2*b16 + 68*b16 + 82*b2*b17 - 66*b17 + 28*b2*b20 + 76*b20
      + 28*b2*b22 + 21*b22 - 74*b2*b29 - 34*b29 - 54*b2*b35 + 53*b35 + 94*b2*
     b36 + 34*b36 - 32*b2*b43 - 61*b43 + 88*b2*b48 + 38*b48 + 80*b2*b51 + 67*
     b51 + 90*b2*b55 + 60*b55 + 32*b2*b56 - 20*b56 + 98*b2*b58 + b58 - 54*b2*
     b62 - 6*b62 + 74*b2*b69 - 86*b69 - 68*b2*b70 - 71*b70 + 26*b2*b72 - 97*b72
      + 66*b2*b75 + 63*b75 - 72*b2*b80 + 19*b80 - 38*b2*b84 - 10*b84 + 90*b2*
     b85 + 25*b85 + 64*b2*b93 + 20*b93 + 38*b2*b101 + 85*b101 - 68*b2*b105 - 88
     *b105 + 20*b2*b106 - 25*b106 - 64*b2*b107 + 27*b107 - 4*b2*b110 - 41*b110
      - 6*b2*b119 + 31*b119 + 82*b2*b120 - 72*b120 - 32*b3*b4 + 48*b3*b5 + 39*
     b5 - 2*b3*b6 - 94*b6 - 98*b3*b7 + 40*b7 - 46*b3*b9 + 2*b3*b12 - 61*b12 - 
     94*b3*b13 - 13*b13 + 12*b3*b17 - 2*b3*b22 + 40*b3*b28 - 74*b28 + 72*b3*b31
      + 58*b31 + 48*b3*b33 + 30*b33 - 54*b3*b37 - 9*b37 - 10*b3*b38 + 21*b38 - 
     76*b3*b41 + 57*b41 - 6*b3*b45 - 100*b45 - 92*b3*b48 + 4*b3*b49 - 54*b49 - 
     10*b3*b52 - 70*b52 - 44*b3*b55 + 34*b3*b57 + 87*b57 - 96*b3*b58 - 60*b3*
     b59 + 92*b59 - 64*b3*b61 - 64*b61 + 76*b3*b62 - 62*b3*b70 + 100*b3*b73 - 
     20*b73 + 94*b3*b90 + 93*b90 - 28*b3*b93 - 96*b3*b104 - 3*b104 - 98*b3*b108
      + 27*b108 - 6*b3*b113 - 29*b113 + 64*b3*b118 + 24*b118 - 80*b3*b121 + 43*
     b121 + 100*b4*b13 - 54*b4*b15 - 41*b15 - 10*b4*b18 - 37*b18 - 66*b4*b19 + 
     15*b19 + 88*b4*b26 - 61*b26 - 86*b4*b30 - 67*b30 - 78*b4*b31 + 80*b4*b34
      + 73*b34 - 4*b4*b35 - 8*b4*b36 - 32*b4*b41 + 50*b4*b43 + 94*b4*b47 - 30*
     b47 + 24*b4*b48 - 48*b4*b49 + 48*b4*b50 + 41*b50 + 6*b4*b52 + 10*b4*b54 - 
     34*b54 - 100*b4*b57 - 56*b4*b59 + 6*b4*b61 - 22*b4*b62 - 40*b4*b63 + 43*
     b63 + 96*b4*b64 - 41*b64 + 64*b4*b72 - 38*b4*b74 - 62*b74 - 66*b4*b77 + 92
     *b77 - 36*b4*b84 - 8*b4*b85 - 66*b4*b92 + 12*b92 + 86*b4*b98 + 59*b98 - 88
     *b4*b104 - 16*b4*b105 + 28*b4*b107 + 72*b5*b7 - 22*b5*b13 + 88*b5*b14 + 78
     *b14 + 88*b5*b18 - 100*b5*b20 - 82*b5*b25 + 78*b25 + 78*b5*b34 - 56*b5*b36
      + 50*b5*b37 - 38*b5*b39 - 41*b39 + 36*b5*b40 - 84*b40 + 16*b5*b45 + 14*b5
     *b48 + 48*b5*b51 + 64*b5*b57 - 10*b5*b62 + 100*b5*b63 - 66*b5*b66 - 23*b66
      + 90*b5*b67 - 47*b67 + 16*b5*b68 - 37*b68 - 76*b5*b69 + 4*b5*b71 - 68*b71
      + 28*b5*b72 - 28*b5*b77 - 56*b5*b78 + 93*b78 + 40*b5*b84 + 40*b5*b87 + 41
     *b87 + 28*b5*b89 - 33*b89 - 54*b5*b96 + 96*b96 - 56*b5*b97 - 28*b97 - 100*
     b5*b99 - 51*b99 + 42*b5*b103 + 76*b103 + 66*b5*b109 - 44*b109 - 38*b5*b110
      + 64*b5*b114 + 24*b114 - 82*b5*b118 + 78*b5*b120 + 42*b5*b121 - 64*b6*b7
      + 44*b6*b9 + 74*b6*b10 + 96*b10 - 96*b6*b17 - 86*b6*b19 - 86*b6*b20 + 50*
     b6*b24 + 16*b24 + 20*b6*b32 - 77*b32 - 28*b6*b33 + 48*b6*b34 + 50*b6*b39
      + 40*b6*b40 + 52*b6*b42 + 39*b42 - 22*b6*b44 - 14*b44 - 92*b6*b47 - 56*b6
     *b52 + 18*b6*b56 + 90*b6*b60 - 21*b60 - 62*b6*b67 - 72*b6*b68 + 44*b6*b72
      - 80*b6*b73 - 64*b6*b74 + 48*b6*b76 + 28*b76 - 54*b6*b78 + 64*b6*b96 + 36
     *b6*b97 + 72*b6*b100 - 33*b100 - 46*b6*b101 - 2*b6*b102 - 9*b102 + 46*b6*
     b107 - 14*b6*b109 - 82*b6*b111 - 68*b111 - 40*b6*b113 - 82*b6*b116 + 55*
     b116 - 24*b6*b117 + 90*b117 - 74*b6*b120 + 14*b7*b8 + 28*b7*b11 - 5*b11 - 
     64*b7*b13 - 44*b7*b23 + 88*b23 - 90*b7*b24 + 98*b7*b25 + 60*b7*b27 - 8*b27
      - 28*b7*b30 - 2*b7*b31 - 86*b7*b32 - 82*b7*b36 + 48*b7*b39 - 76*b7*b40 - 
     30*b7*b46 + 42*b46 + 24*b7*b48 + 20*b7*b53 - 31*b53 - 96*b7*b62 - 82*b7*
     b71 - 40*b7*b74 - 16*b7*b75 + 86*b7*b77 + 96*b7*b78 + 60*b7*b80 - 66*b7*
     b87 + 70*b7*b88 - 64*b88 - 52*b7*b91 + 84*b91 - 76*b7*b92 + 8*b7*b94 - 84*
     b94 - 48*b7*b97 + 14*b7*b98 - 20*b7*b102 + 18*b7*b114 + 4*b7*b117 + 94*b7*
     b121 + 90*b8*b16 - 86*b8*b19 + 74*b8*b20 + 14*b8*b23 - 96*b8*b24 + 34*b8*
     b25 + 28*b8*b41 + 20*b8*b43 - 68*b8*b47 - 46*b8*b55 - 56*b8*b56 + 84*b8*
     b58 + 44*b8*b60 - 80*b8*b62 - 96*b8*b64 + 14*b8*b67 - 80*b8*b71 - 26*b8*
     b73 + 96*b8*b76 - 54*b8*b80 - 76*b8*b83 - 63*b83 - 88*b8*b84 - 22*b8*b90
      - 46*b8*b91 + 16*b8*b92 - 100*b8*b93 + 38*b8*b101 + 48*b8*b107 + 54*b8*
     b109 + 100*b8*b113 - 10*b9*b13 + 10*b9*b17 - 28*b9*b21 - 32*b21 - 68*b9*
     b22 + 28*b9*b24 + 28*b9*b34 + 18*b9*b38 - 34*b9*b39 + 16*b9*b40 - 32*b9*
     b43 + 62*b9*b53 + 68*b9*b61 - 6*b9*b62 + 48*b9*b63 + 14*b9*b66 + 14*b9*b69
      + 68*b9*b70 - 56*b9*b71 - 82*b9*b81 + 6*b81 - 80*b9*b90 - 56*b9*b101 - 8*
     b9*b102 + 24*b9*b106 - 14*b9*b108 - 10*b9*b113 - 14*b9*b114 - 70*b9*b116
      - 12*b9*b118 + 22*b9*b121 + 24*b10*b13 - 94*b10*b16 + 2*b10*b22 - 28*b10*
     b25 + 46*b10*b26 - 60*b10*b27 + 60*b10*b32 + 30*b10*b33 - 32*b10*b42 + 38*
     b10*b47 - 96*b10*b48 - 2*b10*b54 + 6*b10*b56 + 60*b10*b59 + 12*b10*b74 + 
     46*b10*b76 + 22*b10*b86 - 22*b86 - 18*b10*b89 + 22*b10*b90 - 2*b10*b92 + 
     32*b10*b93 + 48*b10*b102 + 78*b10*b103 + 62*b10*b108 + 10*b10*b114 - 28*
     b10*b116 + 20*b10*b118 + 54*b10*b119 + 90*b11*b13 + 80*b11*b20 + 66*b11*
     b22 - 32*b11*b24 + 18*b11*b27 + 44*b11*b30 + 36*b11*b31 - 96*b11*b33 - 44*
     b11*b37 - 26*b11*b38 + 78*b11*b44 + 4*b11*b47 - 26*b11*b51 + 46*b11*b52 + 
     90*b11*b55 - 24*b11*b57 + 66*b11*b58 + 94*b11*b59 - 74*b11*b66 - 4*b11*b67
      - 48*b11*b68 - 58*b11*b69 - 80*b11*b80 - 92*b11*b81 - 56*b11*b82 - 2*b82
      - 12*b11*b84 + 70*b11*b87 - 48*b11*b88 - 2*b11*b91 + 32*b11*b93 + 30*b11*
     b95 - 70*b95 + 96*b11*b96 + 8*b11*b99 + 12*b11*b100 - 78*b11*b101 - 86*b11
     *b103 - 76*b11*b110 - 28*b11*b111 - 86*b11*b114 + 94*b11*b116 + 36*b11*
     b121 + 58*b12*b14 - 96*b12*b21 + 40*b12*b26 + 32*b12*b36 - 38*b12*b40 + 84
     *b12*b42 + 92*b12*b47 - 32*b12*b49 + 94*b12*b50 - 10*b12*b51 + 44*b12*b52
      + 82*b12*b54 - 68*b12*b57 - 6*b12*b58 - 88*b12*b60 + 2*b12*b62 - 38*b12*
     b63 + 22*b12*b66 - 84*b12*b67 - 28*b12*b69 + 92*b12*b71 + 36*b12*b75 + 18*
     b12*b79 + 69*b79 + 18*b12*b80 - 76*b12*b85 + 68*b12*b88 - 98*b12*b90 + 22*
     b12*b94 - 12*b12*b104 - 2*b12*b105 + 32*b12*b107 + 100*b12*b111 - 100*b12*
     b113 + 68*b12*b118 + 36*b12*b120 - 62*b12*b121 + 2*b13*b26 - 76*b13*b32 + 
     98*b13*b35 + 66*b13*b38 + 94*b13*b41 + 74*b13*b43 + 84*b13*b45 - 56*b13*
     b53 + 14*b13*b61 - 54*b13*b65 + 28*b65 + 50*b13*b68 + 36*b13*b75 + 16*b13*
     b78 + 86*b13*b79 - 12*b13*b87 - 26*b13*b88 - 50*b13*b90 + 26*b13*b94 - 24*
     b13*b98 + 98*b13*b108 + 60*b13*b109 + 92*b13*b112 + 35*b112 - 10*b13*b113
      + 12*b13*b121 + 72*b14*b17 + 36*b14*b24 + 10*b14*b26 - 18*b14*b30 + 80*
     b14*b32 - 92*b14*b39 - 84*b14*b41 + 92*b14*b46 - 80*b14*b47 + 40*b14*b48
      - 14*b14*b52 - 34*b14*b56 + 22*b14*b60 + 70*b14*b68 + 28*b14*b71 - 32*b14
     *b72 + 78*b14*b73 - 96*b14*b74 + 84*b14*b75 + 100*b14*b82 - 42*b14*b85 + 
     90*b14*b90 - 52*b14*b91 - 96*b14*b92 + 70*b14*b93 + 4*b14*b94 + 28*b14*b98
      + 86*b14*b99 + 86*b14*b100 + 52*b14*b105 - 90*b14*b114 + 24*b14*b119 + 54
     *b15*b18 + 58*b15*b26 + 50*b15*b29 + 80*b15*b31 - 82*b15*b36 - 44*b15*b37
      - 36*b15*b40 - 90*b15*b43 - 4*b15*b44 - 70*b15*b51 - 26*b15*b56 - 60*b15*
     b58 + 2*b15*b61 - 4*b15*b62 + 22*b15*b66 - 82*b15*b68 - 6*b15*b72 + 34*b15
     *b74 + 54*b15*b75 + 16*b15*b78 + 56*b15*b79 + 82*b15*b80 + 98*b15*b82 + 30
     *b15*b83 + 60*b15*b93 - 74*b15*b95 + 68*b15*b97 + 26*b15*b104 - 50*b15*
     b105 - 34*b15*b110 + 52*b15*b115 + 94*b115 - 90*b15*b116 + 90*b15*b120 - 
     98*b16*b19 - 8*b16*b20 - 74*b16*b24 - 44*b16*b30 - 8*b16*b32 - 100*b16*b39
      + 32*b16*b46 + 28*b16*b48 - 2*b16*b51 + 86*b16*b53 + 8*b16*b60 + 12*b16*
     b70 + 76*b16*b72 + 40*b16*b75 + 54*b16*b82 - 50*b16*b84 + 10*b16*b86 + 6*
     b16*b88 + 40*b16*b91 + 34*b16*b93 + 2*b16*b103 + 36*b16*b104 + 82*b16*b109
      + 36*b16*b110 + 50*b16*b111 - 80*b16*b114 - 94*b16*b115 + 52*b16*b119 + 
     94*b17*b18 + 64*b17*b20 + 12*b17*b21 - 100*b17*b33 + 80*b17*b37 - 78*b17*
     b38 + 94*b17*b39 - 8*b17*b42 + 74*b17*b43 - 78*b17*b45 - 92*b17*b49 - 88*
     b17*b50 - 32*b17*b58 - 80*b17*b67 + 2*b17*b69 + 18*b17*b73 + 32*b17*b76 - 
     50*b17*b83 - 76*b17*b87 - 86*b17*b88 - 2*b17*b89 + 72*b17*b90 + 50*b17*b91
      + 90*b17*b95 + 36*b17*b96 + 94*b17*b97 + 58*b17*b121 - 48*b18*b23 - 70*
     b18*b25 - 32*b18*b27 + 14*b18*b34 + 70*b18*b36 - 50*b18*b37 - 4*b18*b38 - 
     34*b18*b40 + 20*b18*b42 - 56*b18*b43 - 74*b18*b45 + 58*b18*b50 - 46*b18*
     b55 + 86*b18*b62 - 82*b18*b63 - 62*b18*b74 + 4*b18*b76 + 68*b18*b77 - 70*
     b18*b80 - 90*b18*b89 - 16*b18*b91 + 58*b18*b92 + 14*b18*b93 - 22*b18*b95
      + 90*b18*b103 + 32*b18*b105 - 46*b18*b107 + 2*b18*b108 - 6*b19*b22 - 54*
     b19*b24 + 68*b19*b28 - 84*b19*b31 + 6*b19*b32 - 18*b19*b36 - 98*b19*b38 + 
     100*b19*b39 - 6*b19*b45 + 34*b19*b46 - 70*b19*b47 + 26*b19*b50 - 86*b19*
     b52 + 46*b19*b55 - 18*b19*b59 - 80*b19*b61 + 32*b19*b64 + 58*b19*b67 + 76*
     b19*b68 - 40*b19*b74 - 18*b19*b76 + 98*b19*b77 - 30*b19*b84 - 34*b19*b90
      + 10*b19*b96 - 32*b19*b97 - 22*b19*b98 + 48*b19*b104 + 82*b19*b110 + 84*
     b19*b114 + 40*b19*b115 - 6*b19*b117 + 94*b20*b21 - 72*b20*b22 + 12*b20*b25
      + 18*b20*b26 + 66*b20*b30 + 20*b20*b31 + 18*b20*b33 - 98*b20*b34 - 4*b20*
     b43 + 16*b20*b51 - 52*b20*b53 + 20*b20*b62 + 8*b20*b68 - 76*b20*b73 - 88*
     b20*b79 - 44*b20*b81 - 92*b20*b83 - 78*b20*b84 + 88*b20*b88 + 92*b20*b95
      + 88*b20*b99 - 2*b20*b105 + 64*b20*b108 - 30*b20*b110 + 8*b21*b23 + 86*
     b21*b31 - 72*b21*b33 - 62*b21*b34 - 96*b21*b40 - 88*b21*b45 - 40*b21*b47
      - 48*b21*b51 - 24*b21*b53 + 42*b21*b55 - 100*b21*b58 + 44*b21*b61 - 88*
     b21*b64 - 74*b21*b73 + 24*b21*b76 - 26*b21*b77 + 84*b21*b80 + 72*b21*b82
      + 36*b21*b84 - 78*b21*b91 - 28*b21*b97 + 38*b21*b99 + 84*b21*b102 - 38*
     b21*b103 + 56*b21*b105 + 62*b21*b106 + 94*b21*b107 + 92*b21*b109 + 76*b21*
     b112 - 44*b21*b114 + 94*b21*b116 + 30*b21*b118 + 6*b21*b120 + 92*b22*b23
      + 24*b22*b26 + 32*b22*b28 - 38*b22*b30 + 74*b22*b31 + 44*b22*b39 - 70*b22
     *b40 + 12*b22*b42 - 52*b22*b45 + 66*b22*b47 + 34*b22*b48 + 90*b22*b52 + 24
     *b22*b53 - 4*b22*b55 + 14*b22*b57 + 70*b22*b58 + 64*b22*b60 - 70*b22*b61
      - 72*b22*b66 - 72*b22*b79 + 98*b22*b80 - 32*b22*b82 - 68*b22*b89 + 60*b22
     *b94 + 10*b22*b95 + 94*b22*b96 - 98*b22*b99 + 36*b22*b102 - 92*b22*b110 + 
     70*b22*b112 - 10*b22*b120 + 20*b23*b24 + 34*b23*b30 + 62*b23*b42 + 82*b23*
     b43 + 56*b23*b44 + 34*b23*b50 - 100*b23*b51 + 72*b23*b52 + 66*b23*b53 - 78
     *b23*b54 + 40*b23*b57 + 100*b23*b64 + 8*b23*b65 + 88*b23*b67 + 80*b23*b68
      + 92*b23*b71 - 28*b23*b72 + 86*b23*b73 + 2*b23*b81 + 58*b23*b82 - 54*b23*
     b86 + 80*b23*b88 + 62*b23*b90 + 56*b23*b92 - 74*b23*b96 - 8*b23*b100 + 96*
     b23*b101 - 96*b23*b102 + 40*b23*b106 + 22*b23*b107 - 64*b23*b113 - 96*b23*
     b114 - 92*b23*b117 + 88*b24*b31 - 10*b24*b32 + 58*b24*b37 - 12*b24*b39 + 
     18*b24*b41 + 88*b24*b42 - 76*b24*b43 + 6*b24*b44 - 66*b24*b47 + 62*b24*b49
      - 14*b24*b50 - 68*b24*b52 + 28*b24*b53 - 28*b24*b54 + 90*b24*b55 + 8*b24*
     b57 + 42*b24*b63 - 60*b24*b65 + 40*b24*b68 + 100*b24*b69 - 76*b24*b70 - 10
     *b24*b71 - 44*b24*b72 - 12*b24*b73 - 82*b24*b76 + 98*b24*b79 + 52*b24*b82
      - 68*b24*b84 + 20*b24*b85 - 76*b24*b87 + 58*b24*b89 - 6*b24*b90 + 60*b24*
     b93 - 78*b24*b101 - 38*b24*b103 - 20*b24*b106 - 4*b24*b110 - 36*b24*b112
      - 32*b24*b115 + 72*b25*b28 - 96*b25*b41 + 20*b25*b43 + 54*b25*b48 - 90*
     b25*b49 - 46*b25*b52 - 98*b25*b55 - 48*b25*b63 - 22*b25*b68 + 36*b25*b69
      + 36*b25*b73 + 78*b25*b76 + 40*b25*b82 - 88*b25*b87 - 52*b25*b88 - 54*b25
     *b94 - 70*b25*b95 - 66*b25*b98 + 40*b25*b99 - 24*b25*b100 + 76*b25*b102 - 
     26*b25*b106 - 60*b25*b108 - 68*b25*b121 - 8*b26*b28 + 32*b26*b33 + 72*b26*
     b35 + 88*b26*b38 - 22*b26*b43 + 92*b26*b47 - 46*b26*b48 - 10*b26*b52 - 48*
     b26*b54 + 66*b26*b57 - 56*b26*b61 - 26*b26*b64 + 90*b26*b65 - 76*b26*b66
      - 68*b26*b68 + 68*b26*b69 - 52*b26*b78 - 36*b26*b80 + 24*b26*b82 + 2*b26*
     b87 + 12*b26*b90 - 52*b26*b91 + 72*b26*b93 + 16*b26*b96 - 52*b26*b98 - 44*
     b26*b99 + 64*b26*b100 + 36*b26*b106 - 68*b26*b109 + 62*b26*b110 + 30*b26*
     b111 + 76*b26*b112 - 26*b26*b114 - 52*b26*b115 + 46*b27*b28 - 60*b27*b29
      - 52*b27*b35 - 96*b27*b37 + 56*b27*b39 - 52*b27*b41 - 54*b27*b43 + 26*b27
     *b46 + 52*b27*b49 + 70*b27*b55 + 44*b27*b59 - 92*b27*b60 + 16*b27*b62 + 4*
     b27*b67 - 56*b27*b69 - 92*b27*b76 + 48*b27*b80 + 28*b27*b81 - 64*b27*b86
      - 46*b27*b93 - 44*b27*b94 - 18*b27*b102 - 72*b27*b103 + 72*b27*b106 - 100
     *b27*b108 - 94*b27*b110 - 82*b27*b112 - 100*b27*b113 + 4*b27*b114 - 20*b28
     *b29 + 42*b28*b38 + 56*b28*b39 + 68*b28*b41 - 90*b28*b49 - 64*b28*b51 + 4*
     b28*b52 - 86*b28*b53 - 34*b28*b54 + 20*b28*b56 + 48*b28*b57 + 22*b28*b58
      - 36*b28*b59 + 70*b28*b60 - 76*b28*b62 + 66*b28*b63 - 74*b28*b64 - 48*b28
     *b73 - 54*b28*b75 + 60*b28*b79 - 72*b28*b81 + 40*b28*b84 + 32*b28*b87 - 
     100*b28*b89 + 20*b28*b91 - 62*b28*b92 - 40*b28*b93 - 28*b28*b94 - 18*b28*
     b99 - 50*b28*b103 + 98*b28*b105 - 92*b28*b107 + 44*b28*b108 + 42*b28*b110
      + 42*b28*b111 - 20*b28*b113 - 52*b28*b115 - 44*b28*b118 - 56*b28*b121 + 
     22*b29*b32 + 42*b29*b33 - 16*b29*b35 + 34*b29*b36 - 42*b29*b40 - 26*b29*
     b42 + 42*b29*b43 + 60*b29*b48 + 98*b29*b50 - 70*b29*b51 - 54*b29*b52 + 94*
     b29*b59 - 38*b29*b60 - 6*b29*b61 - 16*b29*b70 - 34*b29*b71 + 50*b29*b72 + 
     46*b29*b75 + 82*b29*b76 - 84*b29*b78 - 52*b29*b79 + 38*b29*b92 - 2*b29*b99
      + 34*b29*b100 + 96*b29*b109 - 18*b29*b113 + 76*b29*b114 - 16*b29*b115 - 
     28*b30*b32 + 96*b30*b35 + 12*b30*b36 - 38*b30*b37 + 8*b30*b42 - 86*b30*b45
      - 50*b30*b50 + 80*b30*b51 + 72*b30*b57 + 20*b30*b59 + 74*b30*b60 - 32*b30
     *b67 - 70*b30*b72 - 36*b30*b73 - 18*b30*b77 + 28*b30*b78 - 44*b30*b79 - 74
     *b30*b80 + 64*b30*b81 - 48*b30*b84 - 86*b30*b85 + 50*b30*b97 - 78*b30*b98
      - 62*b30*b100 + 80*b30*b101 - 78*b30*b103 + 62*b30*b105 - 46*b30*b111 + 
     48*b30*b112 + 34*b31*b36 - 34*b31*b40 + 68*b31*b46 + 36*b31*b51 - 94*b31*
     b52 + 94*b31*b59 - 70*b31*b66 + 22*b31*b71 - 88*b31*b74 + 92*b31*b75 - 12*
     b31*b76 + 12*b31*b80 + 64*b31*b85 - 36*b31*b87 - 34*b31*b92 + 36*b31*b96
      + 8*b31*b98 + 44*b31*b101 - 40*b31*b103 + 86*b31*b109 + 34*b31*b117 - 30*
     b31*b119 - 48*b31*b121 - 68*b32*b34 - 10*b32*b38 - 82*b32*b41 + 54*b32*b42
      - 68*b32*b44 - 34*b32*b50 - 86*b32*b51 - 34*b32*b54 + 90*b32*b59 + 86*b32
     *b60 - 28*b32*b62 - 72*b32*b65 + 36*b32*b67 + 2*b32*b68 - 78*b32*b70 + 62*
     b32*b72 + 50*b32*b75 + 2*b32*b78 - 56*b32*b80 + 22*b32*b83 + 36*b32*b87 - 
     34*b32*b89 + 60*b32*b93 + 52*b32*b94 - 36*b32*b97 + 32*b32*b101 - 58*b32*
     b106 - 48*b32*b108 + 34*b32*b116 + 54*b32*b119 + 100*b32*b120 + 48*b33*b34
      - 82*b33*b39 + 70*b33*b42 + 88*b33*b43 - 40*b33*b47 - 42*b33*b50 - 48*b33
     *b52 + 36*b33*b69 + 70*b33*b76 + 76*b33*b77 + 2*b33*b81 - 12*b33*b83 - 78*
     b33*b87 + 30*b33*b91 + 64*b33*b99 - 60*b33*b100 - 98*b33*b109 + 62*b33*
     b114 + 88*b33*b115 + 24*b33*b118 - 96*b33*b121 - 58*b34*b36 - 74*b34*b43
      + 24*b34*b45 + 18*b34*b46 - 94*b34*b48 - 10*b34*b49 - 80*b34*b53 + 2*b34*
     b60 - 82*b34*b61 - 92*b34*b64 + 50*b34*b69 - 86*b34*b72 - 44*b34*b75 + 14*
     b34*b77 + 98*b34*b78 - 14*b34*b80 + 80*b34*b81 - 84*b34*b85 + 44*b34*b87
      + 72*b34*b93 + 46*b34*b97 + 48*b34*b99 - 30*b34*b103 + 100*b34*b111 + 38*
     b34*b115 + 98*b34*b117 + 90*b34*b121 - 62*b35*b36 + 32*b35*b45 - 12*b35*
     b46 - 56*b35*b47 - 58*b35*b50 - 32*b35*b56 + 2*b35*b69 + 52*b35*b73 + 8*
     b35*b77 + 42*b35*b79 + 20*b35*b82 - 40*b35*b98 + 92*b35*b99 + 58*b35*b101
      + 26*b35*b102 - 44*b35*b108 - 12*b35*b114 + 12*b35*b117 + 68*b36*b37 - 40
     *b36*b42 + 38*b36*b43 - 18*b36*b44 + 32*b36*b48 + 10*b36*b50 + 32*b36*b57
      - 62*b36*b58 - 56*b36*b59 + 68*b36*b60 - 38*b36*b61 - 76*b36*b66 + 20*b36
     *b68 - 34*b36*b69 - 16*b36*b74 - 40*b36*b79 - 14*b36*b84 - 6*b36*b89 + 32*
     b36*b94 + 88*b36*b100 + 72*b36*b103 - 100*b36*b105 + 52*b36*b114 + 34*b36*
     b119 - 70*b36*b120 + 14*b37*b39 - 74*b37*b42 + 88*b37*b44 - 72*b37*b45 + 
     58*b37*b48 - 6*b37*b49 - 94*b37*b50 - 18*b37*b51 - 4*b37*b55 + 14*b37*b56
      - 14*b37*b64 - 52*b37*b70 + 16*b37*b71 - 72*b37*b73 - 6*b37*b76 + 54*b37*
     b81 + 70*b37*b93 - 42*b37*b97 - 64*b37*b100 + 22*b37*b107 + 54*b37*b111 - 
     16*b37*b115 + 94*b38*b39 + 54*b38*b43 - 58*b38*b44 + 38*b38*b46 - 38*b38*
     b51 - 54*b38*b55 + 14*b38*b57 - 32*b38*b58 + 94*b38*b61 - 66*b38*b63 - 14*
     b38*b69 - 58*b38*b72 + 8*b38*b74 - 32*b38*b77 - 52*b38*b80 + 48*b38*b82 - 
     38*b38*b83 - 78*b38*b86 - 46*b38*b90 - 42*b38*b91 + 68*b38*b99 + 72*b38*
     b107 - 10*b38*b109 - 38*b38*b111 - 98*b38*b113 - 28*b38*b116 + 44*b38*b119
      - 74*b38*b120 - 40*b39*b41 - 12*b39*b44 - 68*b39*b47 + 50*b39*b52 + 2*b39
     *b53 + 38*b39*b64 - 10*b39*b73 + 56*b39*b78 - 42*b39*b79 - 26*b39*b83 + 58
     *b39*b85 + 50*b39*b92 - 50*b39*b96 + 48*b39*b103 - 4*b39*b109 + 50*b39*
     b110 + 20*b39*b112 + 72*b39*b120 - 86*b40*b41 + 6*b40*b47 + 22*b40*b54 + 
     26*b40*b59 + 68*b40*b60 - 42*b40*b62 + 40*b40*b63 + 52*b40*b66 - 24*b40*
     b70 - 78*b40*b72 + 82*b40*b75 + 76*b40*b84 - 32*b40*b86 + 44*b40*b87 - 20*
     b40*b88 + 10*b40*b89 - 60*b40*b91 + 60*b40*b93 + 38*b40*b95 + 80*b40*b97
      - 96*b40*b103 - 26*b40*b104 - 66*b40*b106 + 88*b40*b109 - 90*b41*b43 - 48
     *b41*b47 + 32*b41*b48 - 34*b41*b50 + 4*b41*b53 + 14*b41*b61 + 64*b41*b64
      + 12*b41*b76 + 8*b41*b81 + 100*b41*b85 - 80*b41*b86 - 64*b41*b90 - 74*b41
     *b93 + 92*b41*b94 + 96*b41*b96 - 52*b41*b103 - 44*b41*b106 - 18*b41*b107
      + 56*b41*b116 - 34*b41*b117 + 50*b41*b119 + 42*b42*b43 - 4*b42*b44 + 90*
     b42*b45 - 96*b42*b49 + 18*b42*b57 + 44*b42*b60 + 76*b42*b62 - 52*b42*b66
      + 52*b42*b67 - 80*b42*b68 + 46*b42*b69 + 66*b42*b72 + 84*b42*b73 - 94*b42
     *b76 - 56*b42*b78 + 68*b42*b81 - 46*b42*b82 + 60*b42*b84 - 20*b42*b85 - 38
     *b42*b86 - 66*b42*b89 - 32*b42*b90 - 90*b42*b95 - 26*b42*b99 - 20*b42*b100
      - 36*b42*b101 - 80*b42*b102 + 44*b42*b106 + 22*b42*b108 - 26*b42*b110 - 
     88*b42*b113 + 22*b42*b114 - 72*b43*b45 - 88*b43*b46 - 52*b43*b49 + 94*b43*
     b53 + 8*b43*b54 - 60*b43*b59 - 96*b43*b62 + 46*b43*b64 - 84*b43*b70 + 24*
     b43*b73 - 12*b43*b80 - 64*b43*b82 + 4*b43*b87 + 88*b43*b89 - 22*b43*b92 + 
     24*b43*b98 - 26*b43*b101 - 10*b43*b102 + 14*b43*b116 - 64*b43*b117 + 12*
     b44*b52 + 70*b44*b62 + 48*b44*b64 - 20*b44*b68 + 100*b44*b71 + 50*b44*b73
      - 38*b44*b74 - 90*b44*b76 - 98*b44*b80 + 90*b44*b82 - 28*b44*b85 - 44*b44
     *b90 + 46*b44*b91 - 60*b44*b102 - 56*b44*b103 - 38*b44*b104 + 78*b44*b113
      - 84*b44*b115 - 6*b44*b117 - 18*b44*b119 - 86*b44*b120 - 100*b45*b48 + 4*
     b45*b51 - 40*b45*b52 - 18*b45*b54 + 20*b45*b56 + 60*b45*b57 + 84*b45*b62
      + 20*b45*b73 + 70*b45*b74 + 96*b45*b76 + 58*b45*b81 - 96*b45*b84 - 20*b45
     *b86 + 78*b45*b87 + 42*b45*b90 + 58*b45*b92 + 92*b45*b96 + 2*b45*b102 + 72
     *b45*b103 - 34*b45*b108 + 28*b45*b110 + 72*b45*b112 - 38*b45*b113 - 52*b45
     *b114 - 34*b45*b121 + 2*b46*b52 + 36*b46*b55 - 62*b46*b58 + 66*b46*b59 + 
     78*b46*b60 - 98*b46*b61 + 54*b46*b68 - 70*b46*b70 + 54*b46*b71 - 2*b46*b72
      + 84*b46*b75 + 52*b46*b77 - 38*b46*b83 + 92*b46*b85 + 32*b46*b87 + 56*b46
     *b89 - 70*b46*b90 + 74*b46*b92 + 8*b46*b95 - 50*b46*b101 - 46*b46*b107 - 
     52*b46*b108 + 66*b46*b111 + 68*b46*b113 - 68*b46*b114 + 56*b46*b120 + 26*
     b47*b48 + 76*b47*b49 - 50*b47*b53 + 88*b47*b54 - 90*b47*b57 + 76*b47*b58
      - 92*b47*b60 + 56*b47*b62 - 98*b47*b63 + 80*b47*b70 + 20*b47*b77 + 50*b47
     *b81 - 90*b47*b82 + 78*b47*b85 - 78*b47*b86 - 68*b47*b87 - 76*b47*b91 - 44
     *b47*b94 - 100*b47*b99 + 82*b47*b101 + 50*b47*b107 - 38*b47*b112 + 24*b47*
     b115 + 34*b47*b116 - 90*b48*b51 - 84*b48*b57 + 36*b48*b61 + 44*b48*b64 - 
     14*b48*b79 - 90*b48*b87 + 62*b48*b89 - 42*b48*b92 + 76*b48*b96 - 70*b48*
     b98 - 80*b48*b101 - 88*b48*b105 - 14*b48*b107 - 84*b48*b109 - 74*b48*b113
      - 46*b48*b118 + 20*b48*b120 - 4*b49*b50 - 26*b49*b54 + 32*b49*b57 + 12*
     b49*b61 + 26*b49*b65 - 26*b49*b70 + 32*b49*b72 + 26*b49*b73 - 70*b49*b74
      + 34*b49*b78 + 100*b49*b81 + 30*b49*b85 + 30*b49*b89 + 20*b49*b94 + 34*
     b49*b96 + 30*b49*b97 + 94*b49*b101 - 30*b49*b102 - 88*b49*b107 - 28*b49*
     b110 + 64*b49*b113 - 48*b49*b114 - 78*b49*b119 + 40*b50*b51 - 78*b50*b53
      - 54*b50*b54 + 36*b50*b57 + 18*b50*b60 + 32*b50*b61 + 94*b50*b63 + 34*b50
     *b64 + 56*b50*b67 - 44*b50*b74 + 76*b50*b84 - 86*b50*b86 + 2*b50*b87 + 92*
     b50*b91 + 88*b50*b92 - 56*b50*b95 - 56*b50*b100 - 66*b50*b104 + 82*b50*
     b106 + 32*b50*b108 - 96*b50*b111 - 26*b50*b118 - 76*b50*b121 + 82*b51*b54
      - 90*b51*b55 - 2*b51*b63 - 34*b51*b66 + 32*b51*b72 + 96*b51*b74 - 52*b51*
     b77 - 32*b51*b78 - 88*b51*b87 + 100*b51*b88 - 20*b51*b89 - 82*b51*b93 - 56
     *b51*b97 - 38*b51*b101 + 6*b51*b113 - 36*b51*b115 + 16*b52*b62 - 62*b52*
     b63 + 90*b52*b65 + 90*b52*b67 + 64*b52*b71 - 62*b52*b83 - 12*b52*b85 - 6*
     b52*b86 - 24*b52*b87 + 26*b52*b101 - 2*b52*b102 + 52*b52*b112 + 42*b52*
     b116 + 38*b53*b55 + 26*b53*b61 - 82*b53*b62 + 88*b53*b67 + 36*b53*b73 - 72
     *b53*b78 - 64*b53*b82 + 88*b53*b85 + 34*b53*b87 + 74*b53*b88 - 96*b53*b93
      + 44*b53*b95 - 34*b53*b96 - 24*b53*b97 - 42*b53*b105 + 36*b53*b106 - 30*
     b53*b108 - 52*b53*b111 + 44*b53*b112 + 60*b53*b114 + 86*b53*b118 + 70*b54*
     b57 + 84*b54*b58 - 98*b54*b60 - 16*b54*b72 - 2*b54*b74 - 76*b54*b76 + 36*
     b54*b78 - 30*b54*b82 + 100*b54*b84 - 90*b54*b85 + 60*b54*b88 - 28*b54*b93
      - 84*b54*b97 - 40*b54*b99 + 34*b54*b101 - 96*b54*b103 + 54*b54*b105 + 84*
     b54*b108 - 86*b54*b109 - 14*b54*b114 + 12*b54*b120 - 2*b55*b57 - 76*b55*
     b60 + 2*b55*b61 - 14*b55*b68 + 82*b55*b69 - 54*b55*b75 - 42*b55*b76 + 32*
     b55*b79 + 78*b55*b83 - 4*b55*b89 - 62*b55*b93 - 42*b55*b102 - 4*b55*b103
      + 6*b55*b105 - 100*b55*b109 - 4*b55*b111 - 50*b55*b117 + 2*b56*b62 - 98*
     b56*b65 - 6*b56*b69 - 72*b56*b70 + 44*b56*b71 - 86*b56*b76 + 60*b56*b79 + 
     28*b56*b83 - 46*b56*b88 + 100*b56*b92 - 74*b56*b93 - 60*b56*b94 - 16*b56*
     b95 + 48*b56*b104 - 92*b56*b105 + 84*b56*b106 - 32*b56*b107 + 78*b56*b110
      + 24*b56*b114 - 92*b56*b115 - 6*b56*b117 - 48*b56*b121 + 94*b57*b59 - 50*
     b57*b60 - 88*b57*b63 + 22*b57*b69 + 8*b57*b74 + 62*b57*b76 - 60*b57*b77 + 
     82*b57*b79 + 66*b57*b88 - 10*b57*b90 - 8*b57*b93 + 16*b57*b98 + 80*b57*b99
      - 54*b57*b101 + 72*b57*b103 - 98*b57*b105 - 8*b57*b115 - 100*b57*b117 + 
     26*b57*b121 - 78*b58*b71 + 24*b58*b75 - 66*b58*b76 + 32*b58*b77 - 26*b58*
     b80 - 14*b58*b82 - 94*b58*b86 - 64*b58*b87 + 18*b58*b92 - 60*b58*b98 + 32*
     b58*b100 - 62*b58*b102 + 16*b58*b105 + 82*b58*b108 - 50*b58*b110 - 6*b58*
     b112 - 22*b58*b113 + 6*b58*b117 + 36*b58*b118 - 34*b58*b119 - 82*b59*b66
      + 58*b59*b67 + 34*b59*b68 + 92*b59*b74 - 28*b59*b77 - 52*b59*b78 + 84*b59
     *b79 - 88*b59*b83 + 2*b59*b84 - 18*b59*b85 - 4*b59*b87 - 82*b59*b88 - 64*
     b59*b99 + 2*b59*b102 + 20*b59*b107 + 14*b59*b108 + 36*b59*b109 - 64*b59*
     b116 - 10*b59*b117 - 38*b59*b120 + 24*b60*b61 - 82*b60*b68 - 12*b60*b70 - 
     62*b60*b72 + 92*b60*b78 - 90*b60*b82 + 80*b60*b83 - 64*b60*b87 + 34*b60*
     b88 - 18*b60*b90 + 96*b60*b92 + 60*b60*b97 - 70*b60*b101 + 56*b60*b102 + 
     82*b60*b110 + 72*b60*b111 + 78*b60*b117 - 54*b60*b118 - 94*b60*b119 - 56*
     b60*b120 + 78*b61*b62 - 78*b61*b63 + 20*b61*b70 + 64*b61*b71 + 96*b61*b72
      - 14*b61*b73 - 56*b61*b80 + 2*b61*b81 + 74*b61*b91 + 80*b61*b92 - 8*b61*
     b95 + 54*b61*b100 + 36*b61*b101 - 42*b61*b102 - 24*b61*b110 + 24*b62*b65
      - 66*b62*b68 + 42*b62*b70 + 90*b62*b74 + 72*b62*b75 - 88*b62*b76 - 52*b62
     *b81 + 62*b62*b83 - 98*b62*b86 + 48*b62*b88 + 14*b62*b92 - 76*b62*b93 - 42
     *b62*b97 - 98*b62*b101 - 58*b62*b103 - 40*b62*b104 - 76*b62*b106 - 20*b62*
     b109 - 94*b62*b110 - 20*b62*b113 + 70*b62*b114 + 80*b62*b115 + 98*b63*b66
      - 100*b63*b68 + 8*b63*b69 + 68*b63*b70 + 14*b63*b72 + 32*b63*b75 - 54*b63
     *b77 - 50*b63*b85 - 6*b63*b87 + 34*b63*b90 - 44*b63*b92 + 2*b63*b95 - 58*
     b63*b102 - 16*b63*b106 - 86*b63*b108 - 14*b63*b110 + 80*b63*b115 + 12*b63*
     b116 + 94*b63*b118 - 96*b63*b120 - 44*b64*b71 - 66*b64*b73 + 92*b64*b75 - 
     66*b64*b77 - 84*b64*b82 + 74*b64*b89 - 100*b64*b90 - 60*b64*b92 + 66*b64*
     b101 - 22*b64*b102 - 64*b64*b104 + 84*b64*b108 + 82*b64*b113 + 48*b64*b115
      + 70*b64*b118 - 88*b64*b120 + 18*b65*b69 + 76*b65*b74 - 32*b65*b77 - 18*
     b65*b85 - 22*b65*b92 - 92*b65*b93 - 94*b65*b95 + 68*b65*b98 - 34*b65*b104
      + 92*b65*b109 + 68*b65*b113 + 28*b66*b71 + 62*b66*b77 + 18*b66*b80 + 20*
     b66*b82 + 34*b66*b83 - 12*b66*b87 + 22*b66*b88 - 18*b66*b90 - 50*b66*b94
      - 16*b66*b96 + 64*b66*b97 + 36*b66*b105 + 2*b66*b109 - 4*b66*b112 + 94*
     b66*b115 - 58*b66*b121 - 34*b67*b70 + 62*b67*b80 - 12*b67*b81 - 96*b67*b87
      - 82*b67*b92 - 30*b67*b97 - 10*b67*b102 + 26*b67*b114 + 50*b68*b69 - 8*
     b68*b71 - 60*b68*b74 + 66*b68*b87 - 68*b68*b91 - 20*b68*b94 - 92*b68*b97
      + 24*b68*b99 - 26*b68*b102 + 80*b68*b104 - 18*b68*b105 + 66*b68*b106 - 90
     *b68*b108 - 60*b68*b110 - 72*b68*b115 - 4*b68*b120 + 100*b69*b75 + 38*b69*
     b76 + 64*b69*b79 + 76*b69*b84 + 28*b69*b85 + 34*b69*b90 + 2*b69*b92 + 40*
     b69*b104 - 94*b69*b111 - 12*b69*b112 - 14*b69*b116 + 98*b69*b119 + 66*b69*
     b120 + 78*b69*b121 + 14*b70*b74 - 12*b70*b75 - 2*b70*b77 - 14*b70*b81 + 4*
     b70*b82 - 18*b70*b83 - 90*b70*b88 - 64*b70*b89 - 70*b70*b90 - 58*b70*b91
      + 26*b70*b95 - 88*b70*b104 - 90*b70*b106 + 16*b70*b107 + 98*b70*b113 - 60
     *b70*b114 - 34*b70*b115 - 84*b70*b117 - 90*b70*b120 + 40*b71*b76 - 6*b71*
     b78 - 60*b71*b79 - 8*b71*b81 - 8*b71*b82 - 10*b71*b84 - 82*b71*b87 - 68*
     b71*b93 + 46*b71*b94 + 42*b71*b95 + 4*b71*b103 - 32*b71*b109 - 8*b71*b111
      - 12*b71*b113 + 28*b71*b118 + 52*b71*b120 - 16*b71*b121 + 70*b72*b76 - 70
     *b72*b91 - 36*b72*b95 + 48*b72*b99 - 78*b72*b110 - 72*b72*b114 - 80*b72*
     b121 - 52*b73*b77 + 2*b73*b78 - 32*b73*b81 + 84*b73*b82 + 8*b73*b89 + 64*
     b73*b92 - 46*b73*b95 + 84*b73*b99 - 82*b73*b101 - 30*b73*b102 - 48*b73*
     b104 + 96*b73*b105 - 100*b73*b113 - 68*b73*b114 + 16*b73*b120 + 58*b74*b77
      - 46*b74*b78 + 62*b74*b84 + 18*b74*b85 - 56*b74*b90 - 30*b74*b91 - 60*b74
     *b93 - 88*b74*b97 - 28*b74*b99 - 68*b74*b100 - 28*b74*b109 + 10*b74*b110
      - 84*b74*b115 - 66*b74*b116 + 60*b75*b77 + 34*b75*b79 - 96*b75*b81 + 4*
     b75*b93 + 100*b75*b94 + 50*b75*b98 + 48*b75*b101 - 60*b75*b102 + 30*b75*
     b103 - 32*b75*b111 + 6*b75*b114 - 82*b75*b115 - 4*b76*b80 + 42*b76*b81 + 
     18*b76*b86 + 62*b76*b92 + 24*b76*b95 + 72*b76*b99 - 52*b76*b102 + 20*b76*
     b103 + 14*b76*b108 + 8*b76*b109 - 14*b76*b111 - 84*b76*b113 + 14*b76*b114
      + 4*b76*b115 + 22*b76*b116 + 14*b76*b121 - 36*b77*b84 - 4*b77*b87 + 12*
     b77*b88 - 68*b77*b90 - 8*b77*b91 + 30*b77*b92 - 34*b77*b95 + 80*b77*b101
      + 52*b77*b102 + 46*b77*b107 - 84*b77*b113 + 56*b77*b114 + 10*b77*b121 - 
     66*b78*b83 + 86*b78*b85 + 60*b78*b89 + 44*b78*b92 - 76*b78*b93 + 90*b78*
     b96 + 4*b78*b112 - 32*b78*b114 + 48*b78*b115 - 50*b78*b116 + 82*b78*b119
      - 8*b78*b121 - 6*b79*b80 + 82*b79*b83 + 26*b79*b87 - 90*b79*b94 - 36*b79*
     b102 - 36*b79*b104 - 48*b79*b107 - 42*b79*b108 - 14*b79*b110 - 12*b79*b117
      - 26*b79*b119 - 60*b80*b81 - 48*b80*b89 - 66*b80*b91 + 28*b80*b97 + 80*
     b80*b104 + 88*b80*b112 - 40*b80*b113 + 10*b80*b114 - 90*b80*b120 + 44*b81*
     b86 + 4*b81*b88 - 20*b81*b91 - 30*b81*b94 - 62*b81*b97 - 48*b81*b99 + 60*
     b81*b112 + 86*b81*b118 - 96*b82*b84 - 92*b82*b85 - 42*b82*b93 + 8*b82*b99
      + 34*b82*b103 - 52*b82*b105 - 24*b82*b107 + 30*b82*b113 + 40*b82*b116 - 
     16*b82*b120 + 92*b83*b86 + 64*b83*b90 - 78*b83*b91 - 32*b83*b93 + 74*b83*
     b95 + 16*b83*b96 + 88*b83*b101 + 68*b83*b102 + 44*b83*b103 + 38*b83*b104
      + 74*b83*b105 - 24*b83*b108 + 70*b83*b109 - 48*b83*b110 + 16*b83*b111 + 
     86*b83*b114 - 2*b83*b116 + 74*b83*b117 - 8*b84*b87 - 50*b84*b88 - 22*b84*
     b92 - 42*b84*b94 - 88*b84*b95 + 52*b84*b98 + 22*b84*b99 - 40*b84*b105 + 38
     *b84*b106 + 94*b84*b115 + 74*b84*b116 + 22*b84*b119 - 24*b85*b89 - 84*b85*
     b95 + 6*b85*b98 + 36*b85*b102 - 32*b85*b105 + 96*b85*b109 - 78*b85*b118 - 
     74*b86*b88 + 76*b86*b91 + 74*b86*b98 + 70*b86*b99 - 96*b86*b100 + 8*b86*
     b101 + 56*b86*b103 + 96*b86*b105 - 78*b86*b112 - 90*b86*b115 + 82*b86*b117
      + 66*b86*b119 + 10*b87*b88 + 90*b87*b92 - 38*b87*b96 + 66*b87*b97 + 88*
     b87*b99 - 64*b87*b109 - 4*b87*b115 + 98*b87*b119 - 86*b88*b91 + 36*b88*b96
      + 22*b88*b98 + 86*b88*b99 + 24*b88*b100 + 32*b88*b102 - 80*b88*b105 + 2*
     b88*b110 - 28*b88*b116 - 16*b88*b121 + 38*b89*b101 - 34*b89*b102 - 70*b89*
     b107 - 78*b89*b109 - 74*b89*b110 + 78*b89*b114 - 54*b89*b115 + 42*b89*b116
      + 16*b90*b93 + 4*b90*b94 + 24*b90*b98 + 96*b90*b102 - 48*b90*b103 + 100*
     b90*b105 + 52*b90*b112 - 66*b90*b113 - 26*b90*b115 + 90*b90*b120 + 54*b91*
     b95 + 6*b91*b96 + 48*b91*b97 + 80*b91*b99 + 96*b91*b100 + 18*b91*b119 + 96
     *b91*b121 - 92*b92*b95 + 18*b92*b98 - 80*b92*b104 - 54*b92*b106 + 48*b92*
     b109 + 26*b92*b114 - 96*b92*b116 - 22*b92*b118 + 26*b92*b120 - 40*b93*b95
      - 98*b93*b98 - 20*b93*b99 + 28*b93*b103 - 80*b93*b106 + 6*b93*b110 + 92*
     b93*b113 + 90*b93*b115 - 50*b93*b118 - 36*b93*b121 + 30*b94*b98 + 26*b94*
     b100 + 62*b94*b101 - 44*b94*b102 - 28*b94*b109 - 70*b94*b110 + 80*b94*b114
      - 38*b94*b118 + 44*b94*b120 - 42*b95*b98 - 28*b95*b102 + 92*b95*b103 - 90
     *b95*b106 - 10*b95*b108 + 6*b95*b111 + 60*b95*b113 - 32*b95*b119 - 12*b95*
     b121 - 20*b96*b97 + 58*b96*b100 + 94*b96*b101 + 82*b96*b109 + 50*b96*b110
      + 90*b96*b112 - 26*b96*b119 + 70*b96*b121 + 40*b97*b98 - 26*b97*b101 + 6*
     b97*b104 + 6*b97*b105 + 26*b97*b109 + 48*b97*b110 - 2*b97*b115 - 26*b97*
     b117 + 100*b97*b119 + 72*b98*b103 - 54*b98*b105 + 16*b98*b109 + 54*b98*
     b110 + 4*b98*b112 - 52*b98*b116 + 66*b98*b117 + 64*b99*b100 - 22*b99*b101
      - 56*b99*b103 + 42*b99*b104 - 20*b99*b106 + 76*b99*b107 - 94*b99*b110 + 
     64*b99*b117 - 86*b99*b119 - 46*b100*b104 - 90*b100*b110 + 36*b100*b114 - 
     48*b100*b115 + 14*b100*b119 + 72*b101*b105 - 16*b101*b107 + 26*b101*b108
      + 84*b101*b109 + 4*b101*b110 + 28*b101*b112 + 24*b101*b116 + 40*b101*b119
      - 38*b101*b121 - 28*b102*b103 + 92*b102*b107 - 92*b102*b108 + 8*b102*b110
      + 36*b102*b114 + 12*b102*b115 - 72*b102*b117 - 14*b103*b107 - 16*b103*
     b108 + 40*b103*b113 - 54*b103*b115 + 28*b103*b118 + 4*b103*b120 + 88*b104*
     b105 + 92*b104*b110 - 12*b104*b113 - 94*b104*b115 - 74*b104*b120 + 66*b105
     *b110 - 54*b105*b111 - 38*b105*b113 - 48*b105*b114 + 98*b105*b116 + 24*
     b105*b120 + 10*b106*b108 + 96*b106*b109 - 28*b106*b111 + 12*b106*b114 + 84
     *b106*b115 + 24*b106*b118 + 42*b107*b110 - 58*b107*b111 - 98*b107*b116 - 
     12*b107*b117 + 44*b107*b118 - 42*b107*b119 + 22*b108*b110 + 68*b108*b112
      + 84*b108*b121 + 2*b109*b110 - 18*b109*b117 + 72*b109*b119 - 54*b109*b120
      - 48*b110*b111 - 16*b110*b114 - 28*b111*b112 + 2*b111*b114 - 18*b112*b113
      - 74*b112*b114 + 94*b112*b115 - 18*b112*b116 - 24*b112*b117 - 24*b112*
     b118 + 70*b112*b121 + 60*b113*b115 + 78*b113*b117 + 22*b114*b116 - 74*b114
     *b118 - 30*b115*b116 + 84*b115*b120 - 94*b115*b121 + 98*b116*b117 - 94*
     b116*b118 + 44*b117*b118 + 50*b117*b121 - 14*b118*b120 - 78*b119*b121
      - objvar =E= 0;

Model qplib_5881 / all /;

qplib_5881.limrow=0; qplib_5881.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

qplib_5881.tolproj = 0.0;
$if not set MIQCP $set MIQCP MIQCP

*
*$onEcho>convert.opt
*dumpgdx 5881.gdx
*GDXQuadratic 1
*$offEcho
*
*option miqcp=convert;
*
*m.optfile=1;
*m.reslim=5;
*
*Solve qplib_5881 using %MIQCP% maximizing objvar;

$batInclude qubo_solve qplib_5881 miqcp max objvar 1

display objvar.l;