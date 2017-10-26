from pylab import *

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

xs = {}
ys = {}

for line in open("timingsw.txt").readlines():
    if line.startswith("1"):
        try:
            n, prec, x, t = map(float, line.split())
        except:
            break
        if (n, prec) not in xs:
            xs[n, prec] = []
            ys[n, prec] = []
        xs[n, prec].append(x)
        ys[n, prec].append(t)

for np in ys:
    ym = median(ys[np])
    for i in range(len(ys[np])):
        ys[np][i] /= ym
    if np[0] == 100000:
        print "h", np[1], max(ys[np]), ys[np].index(max(ys[np]))

def ff(X):
    return [x for x in X]

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(8, 4)

nind = 1
for nn in [100,1000,10000,100000]:
    subplot(2,2,nind)
    nind += 1
    th = 0.5
    for np in sorted(xs.keys()):
        if np[0] == nn:
            plot(ff(xs[np]), ys[np], label="%s" % np[1], linewidth=th, color="black")
    ylim([0,3])
    title("$n = %s$" % nn)
    if nn in [100, 10000]:
        ylabel("Relative time")
    if nn in [10000,100000]:
        xlabel("$k / n$")

fig.subplots_adjust(hspace=0.5)
fig.subplots_adjust(wspace=0.2)

savefig("timeplot.pdf", bbox_inches="tight", dpi=200)
savefig("timeplot.eps", bbox_inches="tight", dpi=200)


