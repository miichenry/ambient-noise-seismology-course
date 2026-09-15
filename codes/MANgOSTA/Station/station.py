
import numpy
import utm
import csv

class Station:

    # Stations coordinates
    staname=[]      # This is taken from stations coordinates
    stalat=[]
    stalon=[]
    stael=[]

    X0=[]
    Y0=[]
    xsta=[]
    ysta=[]

    subnets=[]
    matsub=[]

    def __init__(self, file):

        """
        # Read stations coordinates
        f = open(file, 'r')
        for row in f:
            col=row.split()
            if len(col)>=3:
                self.staname.append(col[0])
                self.stalat.append(float(col[1]))
                self.stalon.append(float(col[2]))
        f.close()
        """

        self.staname = []
        self.stalat = []
        self.stalon = []
        self.stael = []

        f=open(file,'r')
        reader = csv.reader(f, delimiter=',')
        i=0
        for row in reader:
            if len(row)==0: continue
            if i==0:
                subnets=row[4:]
                nsub=len(subnets)
            else:
                self.staname.append(row[0])
                self.stalat.append(float(row[1]))
                self.stalon.append(float(row[2]))
                self.stael.append(float(row[3]))
            i=i+1
        f.close()

        self.stalat = numpy.array(self.stalat)
        self.stalon = numpy.array(self.stalon)
        self.stael = numpy.array(self.stael)
        alat = numpy.mean(self.stalat)
        alon = numpy.mean(self.stalon)

        X0, Y0, Z0, L0 = utm.from_latlon(alat, alon)
        nsta=len(self.staname)
        self.xsta = numpy.empty(nsta)
        self.ysta = numpy.empty(nsta)
        for i in range(nsta):
            x, y, z, l = utm.from_latlon(self.stalat[i], self.stalon[i], Z0)
            self.xsta[i] = x - X0
            self.ysta[i] = y - Y0

        # Subnmet matrix
        self.matsub = numpy.empty((nsta,nsub), dtype=bool)
        f = open(file, 'r')
        reader = csv.reader(f, delimiter=',')
        i=0
        for row in reader:
            if i==0:
                pass
            else:
                for j in range(nsub):
                    if row[j]=='0':
                        self.matsub[i - 1][j] = False
                    else:
                        self.matsub[i - 1][j] = True

            i = i + 1
        f.close()

        self.distfile('pair_distance.dat')

    def getll(self, sta):
        i=self.staname.index(sta)
        return self.stalat[i], self.stalon[i]

    def getxy(self, sta):
        i = self.staname.index(sta)
        return self.xsta[i], self.ysta[i]

    def getdist(self, sta1, sta2):

        x1, y1 = self.getxy(sta1)
        x2, y2 = self.getxy(sta2)

        return numpy.sqrt( (x1-x2)**2 + (y1-y2)**2 )

    def distfile(self, filename):

        print("NSTA=",len(self.staname))

        f=open(filename,'w')
        for i in range(len(self.staname)):
            for j in range(len(self.staname)):
                dist=self.getdist(self.staname[i], self.staname[j])
                if i!=j:
                    f.write('%s_%s %f\n' % (self.staname[i], self.staname[j], dist))
        f.close()
