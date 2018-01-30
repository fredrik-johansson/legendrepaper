from pylab import *

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

data1 = """10 2   6.6e-06     6.8e-06    7e-06    1.92e-05
12 2   8.8e-06     8.3e-06    9.2e-06    5.8e-05
14 2   1.11e-05     1.07e-05    1.16e-05    6.6e-05
16 2   1.38e-05     1.27e-05    1.43e-05    8.2e-05
20 2   1.98e-05     1.93e-05    2.05e-05    9.8e-05
24 2   2.71e-05     2.61e-05    2.79e-05    0.000121
28 2   3.5e-05     3.5e-05    3.6e-05    0.000144
34 3   5.1e-05     6.6e-05    4.8e-05    0.00017
40 4   6.5e-05     0.000112    6.7e-05    0.000218
48 4   9.2e-05     0.000175    9.3e-05    0.000265
58 5   0.000129     0.000287    0.000132    0.00033
70 7   0.000184     0.00065    0.000186    0.00045
84 8   0.000266     0.00084    0.000262    0.00054
100 10   0.00036     0.00127    0.00036    0.00064
120 12   0.00051     0.00183    0.00052    0.00083
144 14   0.00075     0.0031    0.00075    0.00096
172 17   0.00103     0.0043    0.00103    0.00129
206 20   0.00152     0.0061    0.00147    0.00157
248 24   0.00212     0.009    0.0021    0.00207
298 29   0.003     0.0154    0.003    0.00269
358 35   0.0043     0.0216    0.0041    0.0035
430 43   0.0065     0.032    0.0046    0.0046
516 51   0.0104     0.061    0.0057    0.0058
620 62   0.0147     0.084    0.0085    0.0085
744 74   0.0225     0.129    0.0106    0.0106
892 89   0.031     0.203    0.0138    0.0139
1070 107   0.056     0.399    0.0209      329
1284 128   0.084     0.616    0.0271      436
1540 154   0.123     1.158    0.035      351
1848 184   0.194     1.748    0.052      367
2218 221   0.283     3.746    0.068      719
2662 266   0.472     4.155    0.1      767
3194 319   0.773     6.392    0.136      625
3832 383   1.248     9.75    0.2      677
4598 459   2.221     18.765    0.295      1643
5518 551   3.887     24.523    0.436      1390
6622 662   6.127     36.069    0.602      1167
7946 794   10.704     54.753    0.878      1278
9536 953   18.868     86.118    1.334      3566
11444 1144   33.661     127.578    2.21      2568
13732 1373   62.475     203.157    3.51      2250
16478 1647   114.694     nan    5.184      9223372036854775807
19774 1977   207.868     nan    8.278      9223372036854775807
23728 2372   382.774     nan    13.162      9223372036854775807
28474 2847   741.873     nan    20.423      9223372036854775807
34168 3416   1350.64     nan    33.51      9223372036854775807
41002 4100   2507.92     nan    51.401      9223372036854775807"""

data2 = """10 10   6.8e-06     7.6e-06    7.1e-06    1.96e-05
12 12   8.9e-06     9.4e-06    9.3e-06    2.8e-05
14 14   1.14e-05     1.25e-05    1.19e-05    3.6e-05
16 16   1.41e-05     1.57e-05    1.46e-05    5.9e-05
20 20   2.07e-05     2.81e-05    2.12e-05    9.1e-05
24 24   2.8e-05     3.8e-05    2.84e-05    0.000132
28 28   3.7e-05     5.2e-05    3.7e-05    0.000179
34 34   5.3e-05     0.000106    5.2e-05    0.000241
40 40   6.7e-05     0.000186    7e-05    0.0003
48 48   9.8e-05     0.000254    9.9e-05    0.00042
58 58   0.000153     0.00041    0.000154    0.0006
70 70   0.000222     0.00092    0.000225    0.00084
84 84   0.0003     0.00121    0.00031    0.00116
100 100   0.00043     0.00184    0.00043    0.00158
120 120   0.00079     0.00246    0.00081    0.00227
144 144   0.00115     0.0041    0.00115    0.0032
172 172   0.00165     0.0058    0.00163    0.0051
206 206   0.00273     0.0088    0.00275    0.0063
248 248   0.0042     0.013    0.0043    0.0096
298 298   0.0064     0.022    0.0062    0.0127
358 358   0.0102     0.032    0.0126    0.0187
430 430   0.0174     0.05    0.0284    0.0306
516 516   0.0282     0.093    0.042    0.051
620 620   0.048     0.141    0.063    0.062
744 744   0.078     0.213    0.095    0.095
892 892   0.148     0.337    0.146    0.147
1070 1070   0.266     0.853    0.227      1249
1284 1284   0.481     1.167    0.343      1546
1540 1540   0.895     1.694    0.529      1691
1848 1848   1.615     2.591    0.852      1969
2218 2218   3.261     4.347    1.475      2670
2662 2662   5.9     6.469    2.229      3116
3194 3194   10.985     9.866    3.497      3455
3832 3832   19.657     15.129    5.728      4058
4598 4598   38.164     23.361    9.336      5736
5518 5518   69.814     35.924    17.137      6302
6622 6622   132.328     54.212    24.994      7080
7946 7946   256.422     85.587    41.479      8394
9536 9536   472.762     134.089    69.394      12096
11444 11444   875.029     197.038    117.501      12805"""

data3 = """10 100   8.8e-06     1.19e-05    9.1e-06    2.29e-05
12 120   1.39e-05     1.74e-05    1.43e-05    3.6e-05
14 140   1.79e-05     2.51e-05    1.84e-05    4.7e-05
16 160   2.27e-05     3e-05    2.34e-05    5.7e-05
20 200   3.5e-05     5e-05    3.5e-05    8.7e-05
24 240   5e-05     6.6e-05    5e-05    0.000119
28 280   7e-05     9.2e-05    7.1e-05    0.000158
34 340   0.00011     0.000155    0.000112    0.000225
40 400   0.000162     0.00034    0.000164    0.00031
48 480   0.000247     0.00044    0.000252    0.00044
58 580   0.00042     0.00071    0.00043    0.00066
70 700   0.00075     0.0017    0.00074    0.00099
84 840   0.00125     0.00249    0.00126    0.00149
100 1000   0.0021     0.004    0.00222    0.00209
120 1200   0.0039     0.0065    0.0031    0.0031
144 1440   0.0071     0.013    0.0047    0.0047
172 1720   0.013     0.019    0.0072    0.0073
206 2060   0.0238     0.031    0.0111    0.0111
248 2480   0.044     0.05    0.0179    0.0185
298 2980   0.081     0.092    0.0279    0.028
358 3580   0.156     0.183    0.045    0.045
430 4300   0.292     0.273    0.074    0.074
516 5160   0.546     0.512    0.123    0.126
620 6200   1.191     0.697    0.216    0.207
744 7440   1.948     1.152    0.349    0.35
892 8920   3.866     1.949    0.6    0.603
1070 10700   7.485     3.529    1.091      10875
1284 12840   14.636     5.372    1.815      13099
1540 15400   26.27     7.918    3.32      15548
1848 18480   49.705     11.253    5.554      18597
2218 22180   90.779     19.282    9.386      22628
2662 26620   166.583     29.134    16.132      27071
3194 31940   307.592     46.119    27.362      32197
3832 38320   591.934     69.316    48.613      38542
4598 45980   1073.92     108.226    82.342      47114
5518 55180   1895.4     145.912    138.899      55961
6622 66220   3348.37     224.232    251.924      66674
7946 79460   6214.14     341.384    399.46      79904
9536 95360   11379.3     533.752    797.335      97917"""

data4 = """10 64   8.1e-06     9.2e-06    8.5e-06    2.18e-05
12 64   1.07e-05     1.16e-05    1.13e-05    3.1e-05
14 64   1.36e-05     1.6e-05    1.41e-05    4e-05
16 64   1.69e-05     1.89e-05    1.75e-05    4.9e-05
20 64   2.41e-05     3e-05    2.48e-05    7.7e-05
24 64   3.2e-05     5.4e-05    3.3e-05    0.000103
28 64   4.2e-05     7.2e-05    4.3e-05    0.000134
34 64   5.9e-05     0.000115    6e-05    0.000182
40 64   7.9e-05     0.000199    8.1e-05    0.000246
48 64   0.000108     0.000281    0.000109    0.00033
58 64   0.000154     0.00042    0.000156    0.00064
70 64   0.000216     0.00093    0.000221    0.00081
84 64   0.0003     0.00117    0.0003    0.00099
100 64   0.00042     0.00174    0.00042    0.00118
120 64   0.0006     0.00233    0.0006    0.00147
144 64   0.00085     0.0038    0.00086    0.00173
172 64   0.0012     0.0052    0.00121    0.00211
206 64   0.00169     0.0075    0.00173    0.00246
248 64   0.00245     0.0105    0.00245    0.003
298 64   0.0034     0.0172    0.0034    0.0037
358 64   0.005     0.024    0.0049    0.0046
430 64   0.0071     0.035    0.0056    0.0056
516 64   0.0102     0.063    0.0065    0.0065
620 64   0.0147     0.083    0.0085    0.0086
744 64   0.0213     0.126    0.01    0.0102
892 64   0.03     0.188    0.0122    0.0122
1070 64   0.043     0.357    0.0142    0.0144
1284 64   0.062     0.514    0.0166    0.0167
1540 64   0.09     0.938    0.0199    0.02
1848 64   0.13     1.342    0.0251    0.0253
2218 64   0.188     2.284    0.03    0.03
2662 64   0.271     3.13    0.037    0.037
3194 64   0.385     4.73    0.042    0.043
3832 64   0.561     6.662    0.055    0.055
4598 64   0.798     11.09    0.065    0.066
5518 64   1.191     15.76    0.078    0.079
6622 64   1.668     22.473    0.094    0.095
7946 64   2.404     33.128    0.11    0.111
9536 64   3.452     53.909    0.129    0.131
11444 64   4.976     77.029    0.166    0.169
13732 64   7.211     112.118    0.196    0.198
16478 64   10.97     190.867    0.236    0.24
19774 64   15.61     nan    0.301    0.301
23728 64   21.985     nan    0.353    0.351
28474 64   30.643     nan    0.42    0.429
34168 64   46.433     nan    0.521    0.523
41002 64   68.251     nan    0.608    0.612
49202 64   97.388     nan    0.711    0.706
59042 64   135.958     nan    0.897    0.907
70850 64   195.674     nan    1.07    1.078"""

ns = []
ts_rec1 = []
ts_rec2 = []
ts_rec3 = []
ts_rec4 = []
ts_multi1 = []
ts_multi2 = []
ts_multi3 = []
ts_multi4 = []
ts_auto1 = []
ts_auto2 = []
ts_auto3 = []
ts_auto4 = []
ts_norec1 = []
ts_norec2 = []
ts_norec3 = []
ts_norec4 = []

for line in data1.splitlines():
    n, prec, tr, tm, ta, tnr = map(float, line.split())
    ts_rec1.append(tr)
    ts_multi1.append(tm)
    ts_auto1.append(ta)
    if n < 1000:
        ts_norec1.append(tnr)

for line in data2.splitlines():
    n, prec, tr, tm, ta, tnr = map(float, line.split())
    ts_rec2.append(tr)
    ts_multi2.append(tm)
    ts_auto2.append(ta)
    if n < 1000:
        ts_norec2.append(tnr)

for line in data3.splitlines():
    n, prec, tr, tm, ta, tnr = map(float, line.split())
    ts_rec3.append(tr)
    ts_multi3.append(tm)
    ts_auto3.append(ta)
    if n < 1000:
        ts_norec3.append(tnr)

for line in data4.splitlines():
    n, prec, tr, tm, ta, tnr = map(float, line.split())
    ns.append(n)
    ts_rec4.append(tr)
    ts_multi4.append(tm)
    ts_auto4.append(ta)
    if n < 1000:
        ts_norec4.append(tnr)

while len(ts_rec1) < len(ns): ts_rec1.append(float("nan"))
while len(ts_multi1) < len(ns): ts_multi1.append(float("nan"))
while len(ts_auto1) < len(ns): ts_auto1.append(float("nan"))
while len(ts_norec1) < len(ns): ts_norec1.append(float("nan"))

while len(ts_rec2) < len(ns): ts_rec2.append(float("nan"))
while len(ts_multi2) < len(ns): ts_multi2.append(float("nan"))
while len(ts_auto2) < len(ns): ts_auto2.append(float("nan"))
while len(ts_norec2) < len(ns): ts_norec2.append(float("nan"))

while len(ts_rec3) < len(ns): ts_rec3.append(float("nan"))
while len(ts_multi3) < len(ns): ts_multi3.append(float("nan"))
while len(ts_auto3) < len(ns): ts_auto3.append(float("nan"))
while len(ts_norec3) < len(ns): ts_norec3.append(float("nan"))

while len(ts_rec4) < len(ns): ts_rec4.append(float("nan"))
while len(ts_multi4) < len(ns): ts_multi4.append(float("nan"))
while len(ts_auto4) < len(ns): ts_auto4.append(float("nan"))
while len(ts_norec4) < len(ns): ts_norec4.append(float("nan"))

ns = array(ns)
ts_rec1 = array(ts_rec1)
ts_rec2 = array(ts_rec2)
ts_rec3 = array(ts_rec3)
ts_rec4 = array(ts_rec4)
ts_multi1 = array(ts_multi1)
ts_multi2 = array(ts_multi2)
ts_multi3 = array(ts_multi3)
ts_multi4 = array(ts_multi4)
ts_auto1 = array(ts_auto1)
ts_auto2 = array(ts_auto2)
ts_auto3 = array(ts_auto3)
ts_auto4 = array(ts_auto4)
ts_norec1 = array(ts_norec1)
ts_norec2 = array(ts_norec2)
ts_norec3 = array(ts_norec3)
ts_norec4 = array(ts_norec4)

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(8, 4)

subplot(2,2,1)
loglog(ns, ts_multi4 / ts_rec4, linewidth=2.0, color="red", linestyle="--", label="Fast multipoint evaluation")
loglog(ns, ts_auto4 / ts_rec4, linewidth=2.0, color="blue", label="Main algorithm")
loglog(ns, ts_norec4 / ts_rec4, linewidth=2.0, color="purple", linestyle=":", label="Main algorithm without recurrence")
ylim([0.01,10.0])
xlim([1e1,1e5])
grid(True, color="gray")
ylabel("Relative time")
title("$p = 64$")

plt.legend(bbox_to_anchor=(1.05, -2.6), loc='lower center', ncol=1)

subplot(2,2,2)
loglog(ns, ts_multi1 / ts_rec1, linewidth=2.0, color="red", linestyle="--", label="Fast multipoint evaluation")
loglog(ns, ts_auto1 / ts_rec1, linewidth=2.0, color="blue", label="Main algorithm")
loglog(ns, ts_norec1 / ts_rec1, linewidth=2.0, color="purple", linestyle=":", label="Main algorithm without recurrence")
ylim([0.01,10.0])
xlim([1e1,1e5])
grid(True, color="gray")
title("$p = n/10$")

subplot(2,2,3)
loglog(ns, ts_multi2 / ts_rec2, linewidth=2.0, color="red", linestyle="--", label="Fast multipoint evaluation")
loglog(ns, ts_auto2 / ts_rec2, linewidth=2.0, color="blue", label="Main algorithm")
loglog(ns, ts_norec2 / ts_rec2, linewidth=2.0, color="purple", linestyle=":", label="Main algorithm without recurrence")
ylim([0.01,10.0])
xlim([1e1,1e5])
grid(True, color="gray")
ylabel("Relative time")
title("$p = n$")
xlabel("$n$")

subplot(2,2,4)
loglog(ns, ts_multi3 / ts_rec3, linewidth=2.0, color="red", linestyle="--", label="Fast multipoint evaluation")
loglog(ns, ts_auto3 / ts_rec3, linewidth=2.0, color="blue", label="Main algorithm")
loglog(ns, ts_norec3 / ts_rec3, linewidth=2.0, color="purple", linestyle=":", label="Main algorithm without recurrence")
ylim([0.01,10.0])
xlim([1e1,1e5])
grid(True, color="gray")
xlabel("$n$")
title("$p = 10n$")




fig.subplots_adjust(hspace=0.4)
fig.subplots_adjust(wspace=0.2)

savefig("benchplot.pdf", bbox_inches="tight", dpi=200)
savefig("benchplot.eps", bbox_inches="tight", dpi=200)


