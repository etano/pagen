import sys
import os
from numpy import *
import h5py as h5

def cds(group,name,data):
  a = zeros(1,float)
  try:
    if not data.dtype == a.dtype:
      a = array(data)
    else:
      a = data
  except:
    a = array([data])
  ds = group.create_dataset(name,a.shape,a.dtype)
  ds[:] = data

def check(s, sCheck):
  if(s != sCheck):
    print "MISMATCH: expected",sCheck,"got",s
    return(False)
  return(True);


class LongRangeParser:
  def __init__(self,baseName,out_file_name,k_cut):
     self.file_name = baseName
     self.a = h5.File(out_file_name)
     self.k_cut = k_cut
  def ReadYK(self):
     vals = loadtxt(self.file_name+"yk")
     self.b = self.a.create_group("long_range")
     cds(self.b,'type','yk')
     cds(self.b,'k_cut',self.k_cut)
     cds(self.b,'k_points',vals[:,0])
     cds(self.b,'n_k',len(vals[:,0]))
     cds(self.b,'u_k',vals[:,1])
     potgen_inFile=open(self.file_name+"dm")
     a=potgen_inFile.readlines()
     mass_1=float(a[1].split()[2])
     mass_2=float(a[2].split()[2])
     cds(self.b,'mass_1',mass_1)
     cds(self.b,'mass_2',mass_2)
     n_d=int(a[5].split()[5])
     box=zeros(n_d,float)
     cds(self.b,'n_d',n_d)
     box[0]=float(a[5].split()[-1])
     for dim in range(n_d):
       #box[dim]=float(a[5].split()[8+dim])
       box[dim]=box[0] # HACK ONLY CUBIC BOXES
     cds(self.b,'box',box)
     print box
     print mass_1,mass_2
  def Done(self):
    self.a.close()

class Squarer2HDFParser:
  def __init__(self,filename,UseVimage):
    self.f = open(filename,'r')
    self.basename = filename[:-2]
    out_file_name = self.basename + 'h5'
    self.a = h5.File(out_file_name,'w')
    self.numFits = 0
    self.SpitPotential = True
    self.Ucounter = 0
    self.dUcounter = 0
    self.samplingcounter = 0
    self.spec1 = ''
    self.spec2 = ''
    self.UseVimage = UseVimage

  def Done(self):
    print "Done."
    self.f.close()
    self.a.close()

  def next(self):
    w = ''
    c = self.f.read(1)
    # check for EOF
    if(c == ''):
      print "ENCOUNTERED EOF"
      sys.exit()
      return w
    empty = True
    while(empty):
      while(c!=' ' and c!='\t' and c!='\n' and c!=''):
        w += c
        c = self.f.read(1)
        empty = False
      if(empty):
        c = self.f.read(1)
    #print "parsed",w
    return w

  def find(self, target):
    s = self.next()
    while(s != target):
      s = self.next()

  def ProcessSquarerInfo(self):
    ### collect squarer run info
    self.b = self.a.create_group("squarer")
    g = self.next()
    check(g,'UNITS')
    self.c = self.b.create_group("units")
    cds(self.c,'temp',self.next())
    cds(self.c,'length',self.next())
    g = self.next()
    check(g,'TYPE')
    self.c = self.b.create_group("type_0")
    self.spec1 = self.next()
    cds(self.c,'species',self.spec1)
    cds(self.c,'lambda',float(self.next()))
    g = self.next()
    check(g,'TYPE')
    self.c = self.b.create_group("type_1")
    self.spec2 = self.next()
    cds(self.c,'species',self.spec2)
    cds(self.c,'lambda',float(self.next()))

    #### get important stats for remainder of read
    self.find("SQUARER")
    self.next()
    self.next()
    self.next()
    self.numFits = int(self.next())
    cds(self.b,'n_fits',self.numFits)
    self.find("POT")
    if self.UseVimage:
      print self.next(), self.next()
      self.k_cut=float(self.next())
      print "KCUT IS",self.k_cut
      cds(self.b,'k_cut',self.k_cut)
      self.find("VIMAGE")
      self.vimage=float(self.next())
      print "VIMAGE IS",self.vimage
      cds(self.b,'v_image',self.vimage)
  ### end SquarerInfo ###

  ### get potential
  def ProcessPotential(self):
    self.find("RANK")
    print "RANK"
    self.next()
    size = int(self.next())
    self.find("BEGIN")
    self.next()
    self.next()
    u = zeros(size) + 0.
    for i in range(0,size):
      u[i] = float(self.next())
    self.b = self.a.create_group("potential")
    cds(self.b,'data',u)
    # dump potential to ASCII
    if(self.SpitPotential):
      pout = open(self.basename + 'potential.dat', 'w')
      for i in range(0,len(u)):
        pout.write(str(u[i]) + '\n')
      pout.close()
  ### end ProcessPotential ###

  def ProcessU(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    n_u_kj = int(self.next())
    n_tau = int(self.next())

    self.find("GRID")
    self.next()
    gridtype = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    n_max = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)

    # load array from .dm file
    u_kj = zeros([numPts, n_u_kj, n_tau]) + 0.
    for cT in range(0,n_tau):
      for cU in range(0,n_u_kj):
        for cG in range(0,numPts):
          u_kj[cG, cU, cT] = float(self.next())

    taus = zeros(n_tau) + 0.
    tau0 = loTau/2
    for t in range(0,n_tau):
      tau0 *= 2.
      taus[t] = tau0

    print 'n_max:', n_max, 'derv:', derv
    section_title = 'u_kj_' + str(self.Ucounter)
    self.Ucounter += 1
    self.b = self.a.create_group(section_title)
    self.c = self.b.create_group("grid")
    cds(self.c,"n_grid_points",numPts)
    cds(self.c,"type",gridtype)
    cds(self.c,"start",start)
    cds(self.c,"end",end)
    cds(self.b,"n_u_kj",n_u_kj)
    cds(self.b,"n_tau",n_tau)
    cds(self.b,"n_max",n_max)
    cds(self.b,"derv",derv)
    cds(self.b,"Rank",rnk)
    cds(self.b,"taus",taus)
    cds(self.b,"data",u_kj)
  ### end ProcessU ###

  def ProcessdU_dbeta(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    n_u_kj = int(self.next())
    n_tau = int(self.next())

    self.find("GRID")
    self.next()
    gridtype = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    n_max = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)

    # load array from .dm file
    u_kj = zeros([numPts, n_u_kj, n_tau]) + 0.
    for cT in range(0,n_tau):
      for cU in range(0,n_u_kj):
        for cG in range(0,numPts):
          u_kj[cG, cU, cT] = float(self.next())

    taus = zeros(n_tau) + 0.
    tau0 = loTau/2
    for t in range(0,n_tau):
      tau0 *= 2.
      taus[t] = tau0

    print 'n_max:', n_max, 'derv:', derv
    section_title = 'du_kj_dbeta_' + str(self.dUcounter)
    self.dUcounter += 1
    self.b = self.a.create_group(section_title)
    self.c = self.b.create_group("grid")
    cds(self.c,"num_grid_points",numPts)
    cds(self.c,"type",gridtype)
    cds(self.c,"start",start)
    cds(self.c,"end",end)
    cds(self.b,"n_u_kj",n_u_kj)
    cds(self.b,"n_tau",n_tau)
    cds(self.b,"n_max",n_max)
    cds(self.b,"derv",derv)
    cds(self.b,"rank",rnk)
    cds(self.b,"taus",taus)
    cds(self.b,"data",u_kj)
  ### end ProcessdU_dbeta ###

  def ProcessSampling(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    n_u_kj = int(self.next())
    n_tau = int(self.next())

    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    derv = int(self.next())
    print 'derv:', derv

    # load array from .dm file
    u_kj = zeros([numPts, n_u_kj, n_tau]) + 0.
    for cT in range(0,n_tau):
      for cU in range(0,n_u_kj):
        for cG in range(0,numPts):
          u_kj[cG, cU, cT] = float(self.next())

    taus = zeros(n_tau) + 0.
    tau0 = loTau/2
    for t in range(0,n_tau):
      tau0 *= 2.
      taus[t] = tau0

    section_title = 'sampling_' + str(self.samplingcounter)
    self.samplingcounter += 1
    self.b = self.a.create_group(section_title)
    cds(self.b,"n_u_kj",n_u_kj)
    cds(self.b,"n_tau",n_tau)
    cds(self.b,"derv",derv)
    cds(self.b,"taus",taus)
    cds(self.b,"data",u_kj)
  ### end ProcessSampling ###

### end class Squarer2HDFParser ###

def Parse(dmfile):
  basename = dmfile[:-2]
  infilename=basename+"yk"
  UseVimage = False
  if os.path.exists(infilename):
    UseVimage = True

  print 'Squarer Parsing'
  sq = Squarer2HDFParser(dmfile,UseVimage)
  print "Process Squarer Info"
  sq.ProcessSquarerInfo()
  print "Process Potential"
  sq.ProcessPotential()
  print "Process U"
  for sec in range(0,sq.numFits + 1):
    sq.ProcessU()
  print "Process dU/dbeta"
  for sec in range(0,sq.numFits + 1):
    sq.ProcessdU_dbeta()
  print "Species are",sq.spec1,sq.spec2
  print "Process Sampling",
  if(sq.spec1 == sq.spec2):
    sq.ProcessSampling()
    sq.ProcessSampling()
  else:
    print "Skipping sampling table for different species"
  sq.Done()

  out_file_name =basename + 'h5'
  infilename=basename+"yk"
  if os.path.exists(infilename):
    print "About to process the long range file",infilename
    #should check to see if it's consistent somehow?
    longRangeParse=LongRangeParser(basename,out_file_name,sq.k_cut)
    longRangeParse.ReadYK()
    longRangeParse.Done()

def usage():
  print "Usage:  %s dmfile.dm" % os.path.basename(sys.argv[0])

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()
    sys.exit(2)

  Parse(argv[1])

if __name__ == "__main__":
  sys.exit(main())
