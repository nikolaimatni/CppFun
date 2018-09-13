import matplotlib.pyplot as plt

filepath = 'xlog.tr'  

lines = [line.rstrip('\n') for line in open(filepath)]

xvalues = []
xtemp = []

cnt = 0
for line in lines:
    xtemp.append(float(line.rstrip(';')))
    if ";" in line:
        cnt+=1
        xvalues.append(xtemp)
        xtemp = [];

t = range(cnt)
plt.xlabel("Time")
plt.ylabel("State")

for i in range(len(xvalues[1])):
    plt.plot(t,[pt[i] for pt in xvalues],label = 'x(%s)'%(i+1))
plt.legend()


plt.figure()
filepath = 'ulog.tr'

lines = [line.rstrip('\n') for line in open(filepath)]

uvalues = []
utemp = []

cnt = 0
for line in lines:
    utemp.append(float(line.rstrip(';')))
    if ";" in line:
        cnt+=1
        uvalues.append(utemp)
        utemp = [];

t = range(cnt)
plt.xlabel("Time")
plt.ylabel("Input")

for i in range(len(uvalues[1])):
    plt.plot(t,[pt[i] for pt in uvalues],label = 'u(%s)'%(i+1))
plt.legend()
plt.show()
