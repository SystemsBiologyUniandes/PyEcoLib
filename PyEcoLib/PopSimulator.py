import platform
import scipy
import numpy as np
import math
from scipy import integrate
from scipy import optimize as opt
from scipy.stats import gamma
from colorama import init, Fore, Back, Style

from cell import Cell

class PopSimulator:
    def __init__(self, ncells, gr, sb, steps, CV2div = 0, CV2gr = 0, lamb=1, V0array=None, nu=1):
        """
        :param ncells: int
        :param gr: float
        :param sb: float
        :param steps: float
        :param CV2div: float
        :param CV2gr: float
        :param lamb: float
        :param V0array: list
        """

        #self.__title()
        self.__check_errors(ncells, gr, sb, steps, CV2div, CV2gr, lamb, nu)
        self.n = ncells # Number of cells to study
        self.smplt = 0 # Sampling time
        self.gr = gr # Growth  rate
        self.total_steps = steps # Division steps
        self.sb = sb #Initial size
        self.l = lamb
        self.nu=nu
        if lamb ==1:
            self.K = self.total_steps *self.gr/(self.sb)
        else:
            self.K = self.total_steps*self.getk()
        self.CV2div = CV2div
        self.CV2gr = CV2gr



        self.output = "" # string to export data in dynamic simulation
        self.output_size = "" # string to export data in divison strategy

        self.num_steps = 0 # Initial steps
        self.V = self.sb # Cell size
        self.time = 0 # Simulation time

        self.cells = [] # Array of cells
        if hasattr(V0array, "__len__"):
            self.V0arr = V0array
        else:
            self.V0arr = []
        self.initialize_cells(V0array=self.V0arr) #Initialize cells
        self.DivFile=''

    def __title(self):
        """
        Initial title with the name of the project
        :return: None
        """

        if platform.system() == "Windows":
            print(" ___    __     __   _______    ______    _____    __       ___   _____")
            print("|  _ \  \  \  |  | |  _____|  /  ____|  / ___ \  |  |     |   | |  __ \\")
            print("| | \ |  \  \ |  | | |       | /       | /   \ | |  |     |___| | |  \ |")
            print("| |_/ /   \  \|  | | |___    | |       | |   | | |  |      ___  | |__/ /")
            print("|  __/     \__   | |  ___|   | |       | |   | | |  |     |   | |  __  \\")
            print("| |           /  / | |       | |       | |   | | |  |     |   | | |  \  |")
            print("| |       ___/  /  | |_____  | \_____  | \___/ | |  |___  |   | | |__/  |")
            print("|_|      |_____/   |_______|  \______|  \_____/  |______| |___| |______/")
        else:
            print("\x1b[1,32m"+" ___    __     __   _______    ______    _____    __       ___   _____"+'\033[0m')
            print("\x1b[1,32m"+"|  _ \  \  \  |  | |  _____|  /  ____|  / ___ \  |  |     |   | |  __ \\"+'\033[0m')
            print("\x1b[1,32m"+"| | \ |  \  \ |  | | |       | /       | /   \ | |  |     |___| | |  \ |"+'\033[0m')
            print("\x1b[1,32m"+"| |_/ /   \  \|  | | |___    | |       | |   | | |  |      ___  | |__/ /"+'\033[0m')
            print("\x1b[1,32m"+"|  __/     \__   | |  ___|   | |       | |   | | |  |     |   | |  __  \\"+'\033[0m')
            print("\x1b[1,32m"+"| |           /  / | |       | |       | |   | | |  |     |   | | |  \  |"+'\033[0m')
            print("\x1b[1,32m"+"| |       ___/  /  | |_____  | \_____  | \___/ | |  |___  |   | | |__/  |"+'\033[0m')
            print("\x1b[1,32m"+"|_|      |_____/   |_______|  \______|  \_____/  |______| |___| |______/"+'\033[0m')


    def __check_errors(self, ncells, gr, sb, steps, CV2div, CV2gr, lamb,nu):
        """
         it generate an error if some param does not comply with the established
         :param ncells: int
         :param gr: float
         :param sb: float
         :param steps: int
         :param CV2div: float
         :param CV2gr: float
         :param lamb: float
         :return: None
        """
        if not(nu in [1,2]):
            raise NameError('nu must be in [1,2]')
        elif ncells <= 0:
            raise NameError('ncells must be positive')
        elif gr < 0:
            raise NameError('gr must be positive')
        elif sb < 0:
            raise NameError('sb must be positive or zero')
        elif steps < 0:
            raise NameError('steps must be positive or zero')
        elif CV2div < 0:
            raise NameError('CV2div must be positive or zero')
        elif CV2gr < 0:
            raise NameError('CV2gr must be positive or zero')
        elif lamb < 0.5 or lamb > 2:
            raise NameError('lamb must be higher than 0.5 and less than 2')


    def newgr(self,CV2):
        """
        Give a new growth rate
        :param CV2: float
        :return: float
        """

        if CV2 ==0:
            return 1.
        else:
            return np.random.gamma(shape=1/CV2,scale=CV2)

    def newdivpar(self,CV2):
        """
        *
        :param CV2: float
        :return: None
        """
        if CV2 ==0:
            return 0.5
        else:
            beta = 0.5*((1/CV2)-1)
            return np.random.beta(a=beta,b=beta)




    def getsb(self,k):
        """
        *
        :param k: float
        :return: None
        """
        def root(tt):
            return self.multimean(tt,k)-2*tt
        def meansb():
            return opt.bisect(root,0.00001,100000)
        sb = meansb()
        return sb


    def multimean(self,s,k):
        """
        *
        :param s: float
        :param k: float
        :return: None
        """

        sb=s
        def moment(sd):
            return self.rhomulti(sb,sd,k)*sd
        v=integrate.quad(moment, sb, np.inf)[0]
        return v

    def rhomulti(self,sb,sd,k):
        """
        *
        :param sb: float
        :param sd: float
        :param k: float
        :return: None
        """

        n=self.total_steps
        lamb=self.l
        gr=self.gr
        c=n*k/gr
        x=c*((sd**lamb-sb**lamb)/lamb)
        return gamma.pdf(x, n)*c*sd**(lamb-1)

    def opti(self,k):
        """
        *
        :param k: float
        :return: float
        """

        return self.getsb(k)-self.sb
    def getk(self):
        """
        return k when it cannot be calculate with the equation gr/sb
        :return: float
        """

        return opt.bisect(self.opti,0.001,1.5)


    def initialize_cells(self, V0array):
        """
         Give the initial params to the cells
         :param V0array: list
         :return: None
        """
        self.cells=[]
        if len(V0array)!=0:
            idx = 0
            for v in V0array:
                gr = self.newgr(self.CV2gr)
                divpar = self.newdivpar(self.CV2div)
                cell = Cell(idx, v, num_steps=self.total_steps, gr=gr, divpar=divpar, k = gr)
                cell.nextt = self.nextt(v,cell.rv,cell)
                self.cells.append(cell)
                idx += 1
        else:
            for i in range(self.n):
                gr = self.newgr(self.CV2gr)
                divpar = self.newdivpar(self.CV2div)
                cell = Cell(i, self.sb, num_steps=self.total_steps, gr = gr, divpar = divpar, k = gr)
                cell.nextt = self.nextt(self.sb,cell.rv,cell)
                self.cells.append(cell)

    def open_file(self, FileName="./DataSz.csv",DivEventsFile=None):
        
        """
        Here open the file to write the .csv outputs
        :param nameCRM: string
        :return: None
        """
        if hasattr(DivEventsFile, "__len__"):
            self.DivFile = open(DivEventsFile, "w")
            output="Sample,Cell,Mother,MotherSize,BirthTime,Sb,GrowthRate,DivPar\n"
            for m in range(len(self.cells)):
                output+=str(m)+','+str(m)+','+str(m)+','+str(np.nan)+',0.000,'+str(np.round(self.cells[m].V,8))\
                +','+str(np.round(self.cells[m].gr*self.gr,8))+','+str(self.cells[m].dp)+'\n'
            self.DivFile.write(output)
        self.output = ""
        self.file = open(FileName, "w")
        self.output += "Time,Sample,Cell,Size,DivSteps\n"
        self.file.write(self.output)
        kk=0
        for cell in self.cells:
            self.output = ""
            self.output += "0.00,"+str(kk)+","+str(kk)+","
            self.output += str(np.round(cell.get_size(), 4) )+",0\n"
            self.file.write(self.output)
            kk+=1


    def nextt (self,s0,r,cell):
        """
        *
        :param s0: float
        :param r: float
        :param cell: Cell
        :return: None
        """
        mu= (self.gr*cell.gr)
        k= self.K*cell.k
        l=self.l
        return (1/(l*mu))*np.log(1-(l*mu*np.log(r)/(k*s0**l))) 

    def simulate(self,tmax):
        
        """
        This function do all operations
        :param tmax: int
        :return: None
        """
        if self.nu==2:
            raise NameError('This function was designed for nu=1.')
        else:
            for cell in self.cells:
                t=0
                while t<tmax:
                    tt = cell.nextt
                    if ((t+tt) <= tmax):
                        cell.num_steps += 1
                        Vn=cell.V*np.exp(self.gr*cell.gr*tt)
                        if cell.num_steps >= cell.total_steps:
                            dp = self.newdivpar(self.CV2div)
                            gr = self.newgr(self.CV2gr)
                            cell.division(Vn,dp,gr,k=gr)
                        else:
                            cell.change(Vn)
                            cell.rv=np.random.rand()
                        cell.nextt = self.nextt(cell.V,cell.rv,cell)
                    else:
                        Vn = cell.V*np.exp(self.gr*cell.gr*(tmax-t))
                        cell.change(Vn)
                        cell.nextt = cell.nextt - (tmax-t)
                    t += tt

    def szdyn(self, tmax, sample_time, FileName = "./DataSz.csv", DivEventsFile=None):
        self.initialize_cells(self.V0arr) #Initialize cells
        if hasattr(DivEventsFile, "__len__"):
            self.open_file(FileName = FileName, DivEventsFile=DivEventsFile)
        else:
            self.open_file(FileName = FileName)
            """
            *
            :param tmax: int
            :param sample_time: int
            :param nameCRM: string
            :return: None
            """
        if self.nu==2:
            if tmax>9*np.log(2)/self.gr:
                raise NameError('Please select tmax<9*doublingtime')
            elif self.n>=2000:
                raise NameError('Please select ncells<=2000')
        self.smplt = sample_time
        self.time = 0
        nextarray=np.empty(len(self.cells))
        for m in range(len(self.cells)):
            nextarray[m]=self.cells[m].nextt
            self.cells[m].idx=m
        tref=self.smplt
        while self.time<tmax:
            nexttm=np.min(nextarray)
            if nexttm>tref-self.time:
                tt=tref-self.time
                self.time+=tt
                self.output = ""
                for ll in range(len(self.cells)):
                    cell=self.cells[ll]
                    cell.nextt+=-tt
                    nextarray[ll]+=-tt
                    g=cell.gr*self.gr
                    cell.V=cell.V*np.exp(g*tt)
                    "Time,Sample,Cell,Size,DivSteps\n"
                    self.output += str(self.time)+","+str(cell.popidx)+","+str(cell.idx)+","\
                    +str(np.round(cell.V,8))+","+str(cell.num_steps)+"\n"
                self.file.write(self.output)
                tref+=self.smplt
            else:
                m = np.argmin(nextarray)
                cell = self.cells[m]
                cell.num_steps+=1
                for ll in range(len(self.cells)):
                    self.cells[ll].nextt+=-nexttm
                    nextarray[ll]+=-nexttm
                    g=self.cells[ll].gr*self.gr
                    self.cells[ll].V=self.cells[ll].V*np.exp(g*nexttm)
                if cell.num_steps>=self.total_steps:
                    momsz=cell.V
                    sz=cell.V*cell.dp
                    sz2=cell.V*(1-cell.dp)
                    cell.V=sz
                    cell.num_steps=0
                    cell.dp = self.newdivpar(self.CV2div)
                    cell.gr = self.newgr(self.CV2gr)
                    gr=cell.gr
                    dp=cell.dp
                    cell.nextt = self.nextt(sz,np.random.rand(),cell)
                    nextarray[m]=cell.nextt
                    if self.nu==2:
                        gr2 = self.newgr(self.CV2gr)
                        dp2 = self.newdivpar(self.CV2div)
                        idx=len(self.cells)
                        cell2 = Cell(idx, V0=sz2, num_steps=self.total_steps, gr=gr2, divpar=dp2, k = gr2)
                        cell2.popidx=cell.popidx
                        cell2.nextt = self.nextt(sz2,np.random.rand(),cell2)
                        nextarray=np.concatenate((nextarray, [cell2.nextt]), axis=None)
                        self.cells.append(cell2)
                    if hasattr(DivEventsFile, "__len__"):
                        self.DivFile.write(str(cell.popidx)+","+str(int(cell.idx))+","+str(int(cell.idx))+","+str(momsz)\
                                                               +","+str(nexttm+self.time)+","+str(np.round(sz,8))\
                                           +","+str(np.round(gr*self.gr,8))+","+str(dp)+"\n")
                        if self.nu==2:
                            self.DivFile.write(str(cell.popidx)+","+str(int(cell2.idx))+","+str(int(cell.idx))+","+str(momsz)\
                                                               +","+str(nexttm+self.time)+","+str(np.round(sz2,8))\
                                           +","+str(np.round(gr2*self.gr,8))+","+str(dp2)+"\n")
                else:
                    cell.rv=np.random.rand()
                    cell.nextt = self.nextt(cell.V, cell.rv, cell)
                    nextarray[m]=cell.nextt
                self.time+=nexttm
        self.file.close()
        if hasattr(DivEventsFile, "__len__"):
            self.DivFile.close()


    def divstrat(self, tmax, sample_time, nameDSM = "./dataDSM.csv"):
        
        """
        *
        :param tmax: int
        :param sample_time: int
        :param nameDSM: string
        :return: None
        """
        if self.nu==2:
            raise NameError('This function was designed for nu==1.')
        else:
            self.initialize_cells(self.V0arr) #Initialize cells
            self.file_size = open(nameDSM, "w")
            self.file_size.write("S_b,S_d,time\n")
            self.smplt = sample_time
            self.time = 0
            self.open_file()
            self.time = 0
            divarray = np.array([])
            tgt = (tmax/10)
            cnt = 0
            for i in range(len(self.cells)):
                divarray = np.concatenate((divarray,[self.get_ndiv(i)]),axis=0)
            while self.time<tmax:
                self.simulate(self.smplt)
                cnt2 = 0
                self.time += self.smplt
                line = ""
                for cell in self.cells:
                    if self.get_ndiv(i) > divarray[cnt2]:
                        line+=str(self.truncate(cell.Vb, 4))+","+str(self.truncate(cell.Vd, 4))+","+str(self.truncate(self.time, 4))+"\n "
                        divarray[cnt2] = self.get_ndiv(i)
                    cnt2+=1
                self.file_size.write(line)
                cnt +=self.smplt
                if cnt >= tgt:
                    print(str(np.int(100*self.time/tmax))+"%")
                    cnt = 0
    
            self.file_size.close()
    
    def du(self,u,sb,t,dt):
        """
        *
        :param u: array
        :param sb: float
        :param t: int
        :param dt: float
        :return: array
        """
        mu=self.gr
        lamb=self.l
        k=self.K
        v=np.zeros_like(u)
        s=sb*np.exp(mu*t)
        for l in range(len(u)):
            if l==0:
                v[0]=(-k*(s**lamb)*u[0])*dt
            elif l==len(u)-1:
                v[len(u)-1]=(k*(s**lamb)*u[len(u)-2])*dt
            elif l==len(u)-2:
                v[len(u)-2]=(-k*(s**lamb)*u[len(u)-2]+k*(s**lamb)*u[len(u)-3])*dt
            else:
                v[l]=(-k*(s**lamb)*u[l]+k*(s**lamb)*u[l-1])*dt
        return v


    def SdStat(self,sb):
        """
        *
        :param sb: float
        :return: float, float
        """
        mu=self.gr
        tmax=5/self.gr
        dt=0.001/self.gr
        u=np.zeros(self.total_steps+1)
        t=0
        count=10
        plim=[]
        tarrayfsp=[]
        u[0]=1
        while t<tmax:
            u+=self.du(u,sb,t,dt)
            t+=dt
            count+=1
            if count>9:
                plim.append(u[-1])
                tarrayfsp.append(t)
                count=0
        tt=np.array(tarrayfsp)
        h=tt[1]-tt[0]
        rhot=np.diff(plim)/h
        trho=0.5*(tt[1:] + tt[:-1])
        sarray=sb*np.exp(mu*tt)
        ds=np.diff(sarray)
        ss=0.5*(sarray[1:] + sarray[:-1])
        rhos=rhot=np.diff(plim)/ds
        mn=np.trapz(rhos*ss,x=ss)
        var=np.trapz(rhos*(ss)**2,x=ss)
        CV2=(var-mn**2)/(mn-sb)**2
        return mn-sb,CV2


    def szdynFSP(self, tmax, CV2sz = 0, nameFSP = "./dataFSP.csv"):
        """
        *
        :param tmax: int
        :param CV2sz: float
        :param nameFSP: string
        :return: None
        """
        if self.nu==2:
            raise NameError('This function was designed for nu==1.')
        else:
            file = open(nameFSP, "w")
            output = "time,Meansize,VarSize\n"
            nsteps=self.total_steps
            gr=self.gr
            k=self.K
            lamb=self.l
            tmax=tmax
            ndivs=int(1.5*tmax*self.gr/np.log(2))
            dt=0.0001*np.log(2)/self.gr
            if CV2sz==0:
                s0arr=[self.V]
            else:
                s0arr = np.linspace(gamma.ppf(0.001,a=1/CV2sz,scale=self.V*CV2sz),
                    gamma.ppf(0.999, a=1/CV2sz,scale=self.V*CV2sz), 30)
                dx=(s0arr[1]-s0arr[0])
                wgs=[]
                for l in s0arr:
                    wgs.append((gamma.cdf(l+dx/2,a=1/CV2sz,scale=self.V*CV2sz)-gamma.cdf(l-dx/2,a=1/CV2sz,scale=self.V*CV2sz))/dx)

            allp=np.zeros([ndivs,len(s0arr),1000])
            obj=0
            countv0=0
            for v0 in s0arr:
                if obj%3==2:
                    print(str(np.int(100*obj/30))+"%")
                obj+=1
                t=0
                steps=int(np.floor(tmax/dt))
                u=np.zeros([ndivs,nsteps])#(DIVS,STEPS)
                u[0]=np.zeros(nsteps)
                u[0][0]=1#P_00
                time=[]#time array
                count=int(np.floor(tmax/(dt*1000)))-1
                count2=0
                for l in range(steps):
                    utemp=u
                    for n in range(len(utemp)):#n=divs,
                        for m in range(len(utemp[n])):#m=steps
                            arg=lamb*(gr*t-n*np.log(2))
                            if (m==0):#m=steps
                                if(n==0):#n=divs
                                    dun=-k*v0**lamb*np.exp(lamb*gr*t)*(utemp[0][0])
                                    u[n][m]+=dun*dt
                                else:
                                    dun=k*v0**lamb*np.exp(arg)*(2**lamb*utemp[n-1][len(utemp[n])-1]-utemp[n][0])
                                    u[n][m]+=dun*dt
                            elif(m==len(utemp[n])-1):
                                if(n==len(utemp)-1):
                                    dun=k*v0**lamb*np.exp(arg)*(utemp[n][len(utemp[n])-2])
                                    u[n][m]+=dun*dt
                                else:
                                    dun=k*v0**lamb*np.exp(arg)*(utemp[n][m-1]-utemp[n][m])
                                    u[n][m]+=dun*dt
                            else:
                                dun=k*v0**lamb*np.exp(arg)*(utemp[n][m-1]-utemp[n][m])
                                u[n][m]+=dun*dt
                    t+=dt
                    count=count+1
                    if count==int(np.floor(tmax/(dt*1000))):
                        time.append(t)
                        mean=0
                        for ii in range(len(allp)):
                            allp[ii][countv0][count2]=np.sum(u[ii])
                            count=0
                    count2+=1
                countv0=countv0+1
            if CV2sz==0:
                fullmeansz=[]
                fullvarsz=[]
                fulltime=[]
                t=0
                dt=tmax/1000
                for ll in range(len(allp[0][0])):
                    ms=0
                    for ctv0 in range(len(s0arr)):
                        tempms=0
                        for ii in range(ndivs):
                            arg=gr*t-np.log(2)*ii
                            tempms+=np.exp(arg)*allp[ii][ctv0][ll]
                        ms+=s0arr[ctv0]*tempms
                    fullmeansz.append(ms)
                    mvar=0
                    for ctv0 in range(len(s0arr)):
                        tempms=0
                        for ii in range(ndivs):
                            arg=gr*t-np.log(2)*ii
                            tempms+=(ms-s0arr[ctv0]*np.exp(arg))**2*allp[ii][ctv0][ll]
                        mvar+=tempms
                    fullvarsz.append(mvar)
                    fulltime.append(t)
                    t+=dt
            else:
                fullmeansz=[]
                fullvarsz=[]
                fulltime=[]
                t=0
                dt=tmax/1000
                for ll in range(len(allp[0][0])):
                    ms=0
                    for ctv0 in range(len(s0arr)):
                        tempms=0
                        for ii in range(ndivs):
                            arg=gr*t-np.log(2)*ii
                            tempms+=np.exp(arg)*allp[ii][ctv0][ll]
                        ms+=s0arr[ctv0]*tempms*wgs[ctv0]*dx
                    fullmeansz.append(ms)
                    mvar=0
                    for ctv0 in range(len(s0arr)):
                        tempms=0
                        for ii in range(ndivs):
                            arg=gr*t-np.log(2)*ii
                            tempms+=(ms-s0arr[ctv0]*np.exp(arg))**2*allp[ii][ctv0][ll]
                        mvar+=tempms*wgs[ctv0]*dx
                    fullvarsz.append(mvar)
                    fulltime.append(t)
                    t+=dt
            for m in range(len(fullmeansz)):
                output += str(fulltime[m])+","+str(fullmeansz[m])+","+str(fullvarsz[m])+"\n"
            file.write(output)
  
    
    
    
    def get_sz(self, n, cells=[]):
        """
        Give the size of a cell
        :param n: int
        :param cells: list
        :return: float
        """

        if len(cells) > 0:
            return cells[n].V
        else:
            return self.cells[n].V
    def get_ndiv(self, n, cells=[]):
        if len(cells) > 0:
            return cells[n].ndiv
        else:
            return self.cells[n].ndiv
    def get_gr(self, n, cells=[]):
        """
        Give the growth rate of a given index cell
        :param n: int
        :param cells: list
        :return: float
        """
        if len(cells) > 0:
            return cells[n].gr
        else:
            return self.cells[n].gr

    def get_dp(self, n, cells=[]):
        """
        *
        :param n: int
        :param cells: array
        :return: float
        """
        if len(cells) > 0:
            return cells[n].dp
        else:
            return self.cells[n].dp

    def get_next_t(self, n, cells=[]):
        """
        Get the next time
        :param n: int
        :param cells: array
        :return: int
        """
        if len(cells) > 0:
            return cells[n].nextt
        else:
            return self.cells[n].nextt


    def truncate(self, num, ciphers):
        """
        This functions return a number with the n number of ciphers
        :param num: float
        :param ciphers: int
        :return: float
        """
        pos = pow(10.0, ciphers)
        return math.trunc(pos * num)/pos


    def __str__(self):
        out = "Initial Params: {\n   tmax: "+str(self.total_time)+", \n   sample time: "+str(self.smplt)+", \n   ncells: "+str(self.n)+", \n   dt: "+str(self.dt)+", \n   alpha: "+str(self.alpha)+", \n   k: "+str(self.K)+"\n}"
        for cell in self.cells:
            out+= str(cell)+"\n"
        return out
