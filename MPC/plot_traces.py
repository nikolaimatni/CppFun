import matplotlib.pyplot as plt

def plotFile(filepath,label):
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
        plt.plot(t,[pt[i] for pt in xvalues],label = '{}({})'.format(label,i+1))

    plt.legend()
       

plt.figure()
plt.subplot(2,2,1)
plotFile('lti_xlog.tr','x')
plt.subplot(2,2,2)
plotFile('mpc_xlog.tr','x')

plt.subplot(2,2,3)
plotFile('lti_ulog.tr','u')
plt.subplot(2,2,4)
plotFile('mpc_ulog.tr','u')
plt.show()

