#save data in the XRTplot
def save_data(plot):
    dump = []
    x = plot.xaxis.binCenters
    y = plot.yaxis.binCenters
    Ixy = plot.total2D
    dump.append([x, y, Ixy])
    fileName = '{0}'.format(plot.title)
    pickleName = '{0}.pickle'.format(fileName)
    with open(pickleName, 'wb') as f:
        pickle.dump(dump, f, protocol=2)
    print(plot.title," Data save :Done")
