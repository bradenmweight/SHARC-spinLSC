%chk=geometry.chk
%nprocshared=1
%mem=1GB

#p opt b3lyp/sto-3g

TitleMe

0 1
 C                 -0.48500713   -2.13775299   -0.53864796
 C                  0.87261856   -2.08708570   -0.57654145
 C                  1.32994259   -0.60431529   -0.59501175
 C                  0.23694003    0.20299875   -0.56780773
 H                 -1.08435361   -3.02391320   -0.51855879
 H                  1.52371298   -2.93606029   -0.59131369
 H                  2.34597031   -0.27007532   -0.62452628
 H                  0.24091793    1.27298284   -0.57208562
 C                 -1.00067356   -0.69670908   -0.52994143
 C                 -2.33927040   -0.28340946   -0.49438727
 H                 -3.12338180   -1.01105257   -0.46978563
 H                 -2.55579335    0.66590835   -0.49207448






