The Lagrangian equation for x1 is:
m1*diff(x1(t), t, t) - (c1*(2*diff(x0(t), t) - 2*diff(x1(t), t)))/2 + (c2*(2*diff(x1(t), t) - 2*diff(x2(t), t)))/2 - (k1*(2*x0(t) - 2*x1(t)))/2 + (k2*(2*x1(t) - 2*x2(t)))/2
symbolic function inputs: t
 
The Lagrangian equation for x2 is:
m2*diff(x2(t), t, t) - (c2*(2*diff(x1(t), t) - 2*diff(x2(t), t)))/2 + (c3*(2*diff(x2(t), t) - 2*diff(x3(t), t)))/2 - (k2*(2*x1(t) - 2*x2(t)))/2 + (k3*(2*x2(t) - 2*x3(t)))/2
symbolic function inputs: t
 
The Lagrangian equation for x3 is:
(c4*(2*diff(x3(t), t) - 2*L*sin(theta3(t))*diff(theta3(t), t)))/2 - (m3*(2*Lg*cos(theta3(t))*diff(theta3(t), t)^2 - 2*diff(x3(t), t, t) + 2*Lg*sin(theta3(t))*diff(theta3(t), t, t)))/2 - (c3*(2*diff(x2(t), t) - 2*diff(x3(t), t)))/2 - (k3*(2*x2(t) - 2*x3(t)))/2 + (k4*(2*x3(t) - 2*L*cos(theta3_hat) + 2*L*cos(theta3(t))))/2
symbolic function inputs: t
 
The Lagrangian equation for y3 is:
(k5*(2*y3(t) + 2*L*sin(theta3(t)) - 2*L*sin(theta3_hat)))/2 + (m3*(2*Lg*cos(theta3(t))*diff(theta3(t), t, t) - 2*Lg*sin(theta3(t))*diff(theta3(t), t)^2 + 2*diff(y3(t), t, t)))/2 + (c5*(2*diff(y3(t), t) + 2*L*cos(theta3(t))*diff(theta3(t), t)))/2
symbolic function inputs: t
 
The Lagrangian equation for theta3 is:
ct3*diff(theta3(t), t) - (kt3*(2*theta3_hat - 2*theta3(t)))/2 + J3g*diff(theta3(t), t, t) + L*c5*cos(theta3(t))*(diff(y3(t), t) + L*cos(theta3(t))*diff(theta3(t), t)) - L*k4*sin(theta3(t))*(x3(t) - L*cos(theta3_hat) + L*cos(theta3(t))) + L*k5*cos(theta3(t))*(y3(t) + L*sin(theta3(t)) - L*sin(theta3_hat)) + Lg*m3*cos(theta3(t))*(Lg*cos(theta3(t))*diff(theta3(t), t, t) - Lg*sin(theta3(t))*diff(theta3(t), t)^2 + diff(y3(t), t, t)) - L*c4*sin(theta3(t))*(diff(x3(t), t) - L*sin(theta3(t))*diff(theta3(t), t)) + Lg*m3*sin(theta3(t))*(Lg*cos(theta3(t))*diff(theta3(t), t)^2 - diff(x3(t), t, t) + Lg*sin(theta3(t))*diff(theta3(t), t, t))
symbolic function inputs: t
 
 
M(t) =
 
[m1,  0,                     0,                    0,                                                         0]
[ 0, m2,                     0,                    0,                                                         0]
[ 0,  0,                    m3,                    0,                                     -Lg*m3*sin(theta3(t))]
[ 0,  0,                     0,                   m3,                                      Lg*m3*cos(theta3(t))]
[ 0,  0, -Lg*m3*sin(theta3(t)), Lg*m3*cos(theta3(t)), m3*Lg^2*cos(theta3(t))^2 + m3*Lg^2*sin(theta3(t))^2 + J3g]
 
 
C(t) =
 
[c1 + c2,     -c2,                    0,                   0,                                                       0]
[    -c2, c2 + c3,                  -c3,                   0,                                                       0]
[      0,     -c3,              c3 + c4,                   0,                                    -L*c4*sin(theta3(t))]
[      0,       0,                    0,                  c5,                                     L*c5*cos(theta3(t))]
[      0,       0, -L*c4*sin(theta3(t)), L*c5*cos(theta3(t)), c5*L^2*cos(theta3(t))^2 + c4*L^2*sin(theta3(t))^2 + ct3]
 
 
K(t) =
 
[k1 + k2,     -k2,                    0,                   0,                                                                                                                                                                                                 0]
[    -k2, k2 + k3,                  -k3,                   0,                                                                                                                                                                                                 0]
[      0,     -k3,              k3 + k4,                   0,                                                                                                                                                                              -L*k4*sin(theta3(t))]
[      0,       0,                    0,                  k5,                                                                                                                                                                               L*k5*cos(theta3(t))]
[      0,       0, -L*k4*sin(theta3(t)), L*k5*cos(theta3(t)), kt3 + L^2*k5*cos(theta3(t))^2 + L^2*k4*sin(theta3(t))^2 - L*k4*cos(theta3(t))*(x3(t) - L*cos(theta3_hat) + L*cos(theta3(t))) - L*k5*sin(theta3(t))*(y3(t) + L*sin(theta3(t)) - L*sin(theta3_hat))]
 
 
M_at_theta3(t) =
 
[m1,  0,      0,  0,             0]
[ 0, m2,      0,  0,             0]
[ 0,  0,     m3,  0,        -Lg*m3]
[ 0,  0,      0, m3,             0]
[ 0,  0, -Lg*m3,  0, m3*Lg^2 + J3g]
 
 
C_at_theta3(t) =
 
[c1 + c2,     -c2,       0,  0,            0]
[    -c2, c2 + c3,     -c3,  0,            0]
[      0,     -c3, c3 + c4,  0,        -L*c4]
[      0,       0,       0, c5,            0]
[      0,       0,   -L*c4,  0, c4*L^2 + ct3]
 
 
K_at_theta3(t) =
 
[k1 + k2,     -k2,       0,  0,            0]
[    -k2, k2 + k3,     -k3,  0,            0]
[      0,     -k3, k3 + k4,  0,        -L*k4]
[      0,       0,       0, k5,            0]
[      0,       0,   -L*k4,  0, k4*L^2 + kt3]
 
 
Q(t) =
 
(c1 + c2)*diff(x0(t), t) + x0(t)*(k1 + k2)
            - c2*diff(x0(t), t) - k2*x0(t)
                                         0
                                         0
                                         0
 
 
M_substituted(t) =
 
[1/2,     0,           0,     0,                0]
[  0, 27/20,           0,     0,                0]
[  0,     0,       21/10,     0,      -9471/25000]
[  0,     0,           0, 21/10,                0]
[  0,     0, -9471/25000,     0, 5271421/62500000]
 
 
K_substituted(t) =
 
[219410,  -23610,         0,     0,             0]
[-23610,  508210,   -484600,     0,             0]
[     0, -484600,   1000000,     0,     -845256/5]
[     0,       0,         0, 40250,             0]
[     0,       0, -845256/5,     0, 69313617/1250]
 
 
C_substituted(t) =
 
[ 773/4,  -804/5,          0,        0,              0]
[-804/5,  3772/5,    -2968/5,        0,              0]
[     0, -2968/5,     3641/5,        0,     -27593/625]
[     0,       0,          0, 5001/100,              0]
[     0,       0, -27593/625,        0, 6002127/312500]
 

M_num =

    0.5000         0         0         0         0
         0    1.3500         0         0         0
         0         0    2.1000         0   -0.3788
         0         0         0    2.1000         0
         0         0   -0.3788         0    0.0843


K_num =

   1.0e+06 *

    0.2194   -0.0236         0         0         0
   -0.0236    0.5082   -0.4846         0         0
         0   -0.4846    1.0000         0   -0.1691
         0         0         0    0.0403         0
         0         0   -0.1691         0    0.0555


C_num =

  193.2500 -160.8000         0         0         0
 -160.8000  754.4000 -593.6000         0         0
         0 -593.6000  728.2000         0  -44.1488
         0         0         0   50.0100         0
         0         0  -44.1488         0   19.2068


eigenmodes_undamped =

    0.0790   -0.0000    1.4078   -0.1082    0.0101
    0.7159   -0.0000   -0.0106    0.3542   -0.3205
    0.7252   -0.0000   -0.0667   -0.2088    1.3915
    0.0000    0.6901    0.0000    0.0000         0
    2.1934   -0.0000   -0.0097    2.1562    7.2828


omega_squared =

   1.0e+06 *

    0.0109         0         0         0         0
         0    0.0192         0         0         0
         0         0    0.4392         0         0
         0         0         0    0.5934         0
         0         0         0         0    1.9356


eigen_freq_undamped =

   1.0e+03 *

    0.1043         0         0         0         0
         0    0.1384         0         0         0
         0         0    0.6627         0         0
         0         0         0    0.7703         0
         0         0         0         0    1.3913


eigenmodes_damped =

  Columns 1 through 6

   0.0125 + 0.0092i   0.0125 - 0.0092i  -0.3323 - 0.7694i  -0.3323 + 0.7694i  -0.0195 + 0.1071i  -0.0195 - 0.1071i
  -0.0593 - 0.0096i  -0.0593 + 0.0096i  -0.0311 + 0.0739i  -0.0311 - 0.0739i  -0.0876 + 0.0146i  -0.0876 - 0.0146i
   0.2133 + 0.0134i   0.2133 - 0.0134i  -0.0971 + 0.1371i  -0.0971 - 0.1371i   0.0123 - 0.0818i   0.0123 + 0.0818i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.9749 - 0.0036i   0.9749 + 0.0036i   0.0789 + 0.5067i   0.0789 - 0.5067i  -0.9730 - 0.1633i  -0.9730 + 0.1633i

  Columns 7 through 10

   0.0140 + 0.0285i   0.0140 - 0.0285i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.2694 + 0.1181i   0.2694 - 0.1181i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.2730 + 0.1198i   0.2730 - 0.1198i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i  -0.0863 - 0.9963i  -0.0863 + 0.9963i
   0.8364 + 0.3521i   0.8364 - 0.3521i   0.0000 + 0.0000i   0.0000 + 0.0000i


eigen_freq_damped =

   1.0e+03 *

  -1.0902 + 0.7818i
  -1.0902 - 0.7818i
  -0.2029 + 0.6218i
  -0.2029 - 0.6218i
  -0.1493 + 0.7706i
  -0.1493 - 0.7706i
  -0.0466 + 0.0969i
  -0.0466 - 0.0969i
  -0.0119 + 0.1379i
  -0.0119 - 0.1379i

