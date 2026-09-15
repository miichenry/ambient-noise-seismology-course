
class Parser:

    ifile = []
    counter = 0

    def __init__(self, file):

        nlines = sum(1 for line in open(file))

        ###############################
        # Preprocessing input file

        self.ifile = []
        rem = False
        f = open(file, 'r')
        for i in range(nlines):
            row = f.readline().split()
            nf = len(row)

            if len(row) == 0:
                continue

            if row[0] == "#<" and rem == False:
                rem = True
                continue

            if row[0] == "#>" and rem == True:
                rem = False
                continue

            if rem == True:
                continue

            if row[0][0] == '#':
                continue

            self.ifile.append(row)

        f.close()

        ###############################
        # Variable substitution

        var={}
        for k in range(len(self.ifile)):

            # Replace
            for j in range(len(self.ifile[k])):

                flagv=False
                for l in range(len(self.ifile[k][j])):
                    if self.ifile[k][j][l]=="$":
                        if self.ifile[k][j][l+1]=="{":
                            flagv=True
                            vstart=l+2

                    if self.ifile[k][j][l]=="}" and flagv==True:
                        vend=l
                        varname=self.ifile[k][j][vstart:vend]
                        flagv=False
                        if varname in var.keys():
                            self.ifile[k][j]=self.ifile[k][j].replace(self.ifile[k][j][vstart-2:vend+1],var[varname])

            # Upgrade
            if self.ifile[k][0].upper()=="VAR":
                name = self.ifile[k][1]
                val = self.ifile[k][2]
                var[name]=val

        ###############################
        # Delete VAR
        k=0
        while k<len(self.ifile):
            if self.ifile[k][0].upper() == "VAR":
                del self.ifile[k]
            else:
                k=k+1

        self.counter=0

    def get(self, rewind=False):

        if rewind: self.counter = self.counter - 1

        row=self.ifile[self.counter]
        self.counter = self.counter + 1

        return row

    def finished(self):
        return self.counter>=len(self.ifile)

    def getNrows(self):

        return len(self.ifile)
