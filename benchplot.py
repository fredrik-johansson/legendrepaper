from pylab import *

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

data1 = """10 2   5.8e-06     6.8e-06    6.2e-06    1.3e-05
12 2   7.7e-06     8.2e-06    8.1e-06    5.6e-05
14 2   9.6e-06     1.06e-05    1.01e-05    6.3e-05
16 2   1.16e-05     1.26e-05    1.22e-05    7.8e-05
20 2   1.71e-05     1.9e-05    1.78e-05    9.6e-05
24 2   2.35e-05     2.55e-05    2.45e-05    0.000117
28 2   3e-05     3.4e-05    3.1e-05    0.000142
34 3   4.5e-05     6.5e-05    4.1e-05    0.000169
40 4   5.6e-05     0.00011    5.7e-05    0.000218
48 4   7.7e-05     0.000173    7.9e-05    0.000259
58 5   0.000109     0.000275    0.000111    0.00033
70 7   0.000155     0.0007    0.000161    0.00043
84 8   0.000217     0.00083    0.00022    0.00052
100 10   0.0003     0.00127    0.0003    0.00063
120 12   0.00042     0.00182    0.00043    0.0008
144 14   0.0006     0.0031    0.0006    0.00099
172 17   0.00085     0.0044    0.00086    0.00132
206 20   0.00121     0.0063    0.00123    0.00161
248 24   0.00176     0.0092    0.00177    0.00208
298 29   0.00254     0.0157    0.0025    0.00211
358 35   0.0037     0.0223    0.0034    0.00273
430 43   0.0052     0.033    0.0035    0.0035
516 51   0.0081     0.065    0.0047    0.0048
620 62   0.0121     0.085    0.007    0.007
744 74   0.0176     0.158    0.0088    0.0088
892 89   0.0259     0.239    0.0116    0.0116
1070 107   0.047     0.453    0.0158    0.0158
1284 128   0.071     0.696    0.023    0.0231
1540 154   0.101     1.31    0.03    0.03
1848 184   0.163     1.946    0.043    0.043
2218 221   0.244     3.063    0.058    0.058
2662 266   0.398     4.246    0.082    0.082
3194 319   0.643     6.703    0.119    0.119
3832 383   1.049     10.061    0.171    0.171
4598 459   1.674     16.117    0.247    0.248
5518 551   2.939     24.15    0.358    0.359
6622 662   4.874     35.207    0.537    0.538
7946 794   8.551     52.394    0.797    0.797
9536 953   15.959     85.679    1.242    1.239
11444 1144   29.542     126.174    1.931    1.925
13732 1373   56.005     191.697    3.051    3.046
16478 1647   103.252     nan    4.751    4.715
19774 1977   187.652     nan    7.655    7.647
23728 2372   355.529     nan    12.259    12.233
28474 2847   668.099     nan    19.971    20
34168 3416   1243.2     nan    31.925    31.874
41002 4100   2386.57     nan    51.593    52.516"""

data2 = """12 12   8.1e-06     9.2e-06    8.5e-06    2.29e-05
14 14   1.02e-05     1.27e-05    1.06e-05    2.99e-05
16 16   1.3e-05     1.54e-05    1.36e-05    5.4e-05
20 20   1.86e-05     2.8e-05    1.92e-05    8.7e-05
24 24   2.49e-05     3.8e-05    2.58e-05    0.000128
28 28   3.1e-05     5e-05    3.2e-05    0.00016
34 34   4.5e-05     0.000103    4.6e-05    0.000233
40 40   5.8e-05     0.000178    6e-05    0.000297
48 48   7.9e-05     0.000242    8.1e-05    0.00041
58 58   0.000124     0.0004    0.000127    0.00058
70 70   0.000184     0.00094    0.000187    0.00083
84 84   0.000256     0.00124    0.000261    0.00114
100 100   0.00036     0.0019    0.00037    0.00156
120 120   0.00065     0.00249    0.00066    0.00222
144 144   0.00093     0.0042    0.00095    0.0032
172 172   0.00129     0.0058    0.0013    0.0043
206 206   0.00213     0.0086    0.00221    0.0061
248 248   0.0033     0.0127    0.0034    0.0093
298 298   0.005     0.0223    0.0051    0.0117
358 358   0.0081     0.032    0.0083    0.0173
430 430   0.0122     0.049    0.0126    0.0258
516 516   0.0231     0.091    0.0229    0.039
620 620   0.035     0.137    0.035    0.059
744 744   0.065     0.207    0.063    0.09
892 892   0.124     0.33    0.14    0.14
1070 1070   0.226     0.944    0.22    0.22
1284 1284   0.42     1.281    0.335    0.336
1540 1540   0.799     1.87    0.524    0.522
1848 1848   1.456     3.867    0.835    0.834
2218 2218   2.815     4.605    1.378    1.379
2662 2662   5.294     7.035    2.223    2.271
3194 3194   9.868     10.328    3.507    3.534
3832 3832   18.306     15.51    5.758    5.744
4598 4598   34.75     23.865    9.305    9.384
5518 5518   66.447     37.164    15.172    15.042
6622 6622   122.528     54.146    24.723    24.625
7946 7946   239.458     83.98    40.644    40.448
9536 9536   443.605     130.391    67.667    67.469
11444 11444   824.06     188.683    118.243    117.129"""

data3 = """10 100   7.9e-06     1.13e-05    8.3e-06    1.63e-05
12 120   1.06e-05     1.66e-05    1.11e-05    2.93e-05
14 140   1.55e-05     2.35e-05    1.6e-05    3.7e-05
16 160   2.02e-05     3e-05    2.09e-05    4.6e-05
20 200   3.1e-05     5.1e-05    3.2e-05    7.2e-05
24 240   4.3e-05     6.8e-05    4.4e-05    9.8e-05
28 280   5.9e-05     9.1e-05    6e-05    0.000132
34 340   9.3e-05     0.000153    9.5e-05    0.000187
40 400   0.000134     0.00033    0.000135    0.000258
48 480   0.000209     0.00044    0.000213    0.00036
58 580   0.00035     0.00071    0.00035    0.00053
70 700   0.0006     0.00169    0.0006    0.00081
84 840   0.00104     0.00249    0.00105    0.00122
100 1000   0.00179     0.0042    0.00177    0.00177
120 1200   0.0032     0.0066    0.00273    0.00272
144 1440   0.0062     0.0131    0.0041    0.0041
172 1720   0.0112     0.0235    0.0063    0.0063
206 2060   0.021     0.036    0.0098    0.0098
248 2480   0.038     0.058    0.016    0.016
298 2980   0.073     0.106    0.0252    0.0252
358 3580   0.141     0.208    0.041    0.04
430 4300   0.27     0.31    0.067    0.067
516 5160   0.52     0.543    0.114    0.114
620 6200   0.973     0.741    0.19    0.191
744 7440   1.852     1.197    0.328    0.328
892 8920   3.433     1.854    0.566    0.567
1070 10700   6.5     3.194    0.96    0.96
1284 12840   12.518     4.671    1.683    1.647
1540 15400   26.466     6.906    2.8    2.806
1848 18480   44.515     10.387    4.854    4.855
2218 22180   82.641     17.183    8.424    8.392
2662 26620   153.746     25.696    14.612    14.552
3194 31940   277.671     40.791    24.336    24.241
3832 38320   528.184     61.862    42.009    41.48
4598 45980   944.251     93.854    77.26    71.299
5518 55180   1751.25     146.992    121.816    122.187
6622 66220   3287.19     225.087    251.29    249.521
7946 79460   6090.74     344.73    434.731    439.657
9536 95360   11161.1     541.502    735.766    736.776"""


# n, prec, recurrence, multipoint, auto, norec
data4 = """10 64   7.1e-06     9.5e-06    7.5e-06    1.67e-05
12 64   9.2e-06     1.18e-05    9.7e-06    2.49e-05
14 64   1.15e-05     1.63e-05    1.21e-05    3.2e-05
16 64   1.4e-05     1.9e-05    1.46e-05    3.9e-05
20 64   2.07e-05     3e-05    2.16e-05    6.1e-05
24 64   2.86e-05     5.4e-05    2.96e-05    8.3e-05
28 64   3.6e-05     7.2e-05    3.7e-05    0.000109
34 64   5e-05     0.000115    5.2e-05    0.00015
40 64   6.6e-05     0.000197    6.8e-05    0.000203
48 64   9.1e-05     0.000282    9.3e-05    0.000277
58 64   0.000127     0.00041    0.000129    0.00061
70 64   0.000179     0.00093    0.000182    0.0008
84 64   0.000248     0.00116    0.000252    0.00098
100 64   0.00034     0.00175    0.00035    0.00119
120 64   0.00048     0.00231    0.00048    0.00146
144 64   0.00068     0.0038    0.00069    0.00173
172 64   0.00097     0.0052    0.00097    0.00213
206 64   0.00137     0.0073    0.00143    0.00251
248 64   0.00197     0.0104    0.00201    0.0031
298 64   0.00283     0.0171    0.00283    0.0031
358 64   0.0041     0.0238    0.0041    0.0038
430 64   0.0058     0.035    0.0046    0.0046
516 64   0.0084     0.063    0.0057    0.0057
620 64   0.0121     0.083    0.0071    0.0072
744 64   0.0174     0.123    0.0082    0.0082
892 64   0.025     0.186    0.0097    0.0098
1070 64   0.036     0.354    0.0117    0.0117
1284 64   0.051     0.497    0.014    0.014
1540 64   0.074     0.922    0.0166    0.0167
1848 64   0.106     1.306    0.0195    0.0196
2218 64   0.154     2.173    0.0233    0.0234
2662 64   0.223     3.017    0.0281    0.0284
3194 64   0.324     4.624    0.033    0.034
3832 64   0.467     6.577    0.04    0.04
4598 64   0.675     10.711    0.048    0.049
5518 64   0.971     15.132    0.058    0.058
6622 64   1.398     22.012    0.07    0.07
7946 64   2.02     32.606    0.083    0.084
9536 64   2.918     52.19    0.101    0.102
11444 64   4.218     74.917    0.122    0.124
13732 64   6.096     108.946    0.147    0.148
16478 64   8.895     1e+300    0.177    0.179
19774 64   12.839     1e+300    0.212    0.214
23728 64   18.543     1e+300    0.255    0.259
28474 64   26.788     1e+300    0.308    0.31
34168 64   38.583     1e+300    0.368    0.372
41002 64   55.653     1e+300    0.443    0.447
49202 64   80.224     1e+300    0.536    0.543
59042 64   115.184     1e+300    0.648    0.655
70850 64   166.467     1e+300    0.784    0.788"""

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


