Minimize
  obj: 0
Subject To
 
\ Row assignment constraints
+ x0x0 + x0x1 + x0x2 + x0x3 + x0x4 = 1
+ x1x0 + x1x1 + x1x2 + x1x3 + x1x4 = 1
+ x2x0 + x2x1 + x2x2 + x2x3 + x2x4 = 1
+ x3x0 + x3x1 + x3x2 + x3x3 + x3x4 = 1
+ x4x0 + x4x1 + x4x2 + x4x3 + x4x4 = 1
 
\ Column assignment constraints
+ x0x0 + x1x0 + x2x0 + x3x0 + x4x0 = 1
+ x0x1 + x1x1 + x2x1 + x3x1 + x4x1 = 1
+ x0x2 + x1x2 + x2x2 + x3x2 + x4x2 = 1
+ x0x3 + x1x3 + x2x3 + x3x3 + x4x3 = 1
+ x0x4 + x1x4 + x2x4 + x3x4 + x4x4 = 1
 
\multiple objectives
 +  6 x0x0 +  1 x0x1 + 20 x0x2 +  2 x0x3 +  3 x0x4 +  2 x1x0 +  6 x1x1 +  9 x1x2 + 10 x1x3 + 18 x1x4 +  1 x2x0 +  6 x2x1 + 20 x2x2 +  5 x2x3 +  9 x2x4 +  6 x3x0 +  8 x3x1 +  6 x3x2 +  9 x3x3 +  6 x3x4 +  7 x4x0 + 10 x4x1 + 10 x4x2 +  6 x4x3 +  2 x4x4 < 1
 + 17 x0x0 + 20 x0x1 +  8 x0x2 +  8 x0x3 + 20 x0x4 + 10 x1x0 + 13 x1x1 +  1 x1x2 + 10 x1x3 + 15 x1x4 +  4 x2x0 + 11 x2x1 +  1 x2x2 + 13 x2x3 +  1 x2x4 + 19 x3x0 + 13 x3x1 +  7 x3x2 + 18 x3x3 + 17 x3x4 + 15 x4x0 +  3 x4x1 +  5 x4x2 +  1 x4x3 + 11 x4x4 < 2
 + 10 x0x0 +  7 x0x1 +  1 x0x2 + 19 x0x3 + 12 x0x4 +  2 x1x0 + 15 x1x1 + 12 x1x2 + 10 x1x3 +  3 x1x4 + 11 x2x0 + 20 x2x1 + 16 x2x2 + 12 x2x3 +  9 x2x4 + 10 x3x0 + 15 x3x1 + 20 x3x2 + 11 x3x3 +  7 x3x4 +  1 x4x0 +  9 x4x1 + 20 x4x2 +  7 x4x3 +  6 x4x4 < 3

BINARY
x0x0
x0x1
x0x2
x0x3
x0x4
x1x0
x1x1
x1x2
x1x3
x1x4
x2x0
x2x1
x2x2
x2x3
x2x4
x3x0
x3x1
x3x2
x3x3
x3x4
x4x0
x4x1
x4x2
x4x3
x4x4
End
