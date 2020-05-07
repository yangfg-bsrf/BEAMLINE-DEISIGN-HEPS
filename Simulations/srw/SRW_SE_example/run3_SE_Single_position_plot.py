#%%
try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

if not srwl_uti_proc_is_master(): exit()
import matplotlib.pyplot as plt
####################################################
# LIGHT SOURCE
#%%
part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.2
part_beam.partStatMom1.x     = 0.0
part_beam.partStatMom1.y     = 0.0
part_beam.partStatMom1.z     = -2.6250000000000004
part_beam.partStatMom1.xp    = 0.0
part_beam.partStatMom1.yp    = 0.0
part_beam.partStatMom1.gamma = 11741.70710144324
part_beam.arStatMom2[0]      = 2.7888999999999995e-10
part_beam.arStatMom2[1]      = 0.0
part_beam.arStatMom2[2]      = 2.7225e-12
part_beam.arStatMom2[3]      = 2.6010000000000002e-11
part_beam.arStatMom2[4]      = 0.0
part_beam.arStatMom2[5]      = 2.809e-13
part_beam.arStatMom2[10]     = 1.2321000000000002e-06

magnetic_fields = []
magnetic_fields.append(SRWLMagFldH(1, 'v', 
                                   _B=0.6012270899961967, 
                                   _ph=0.0, 
                                   _s=-1, 
                                   _a=1.0))
magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.035, _nPer=142.0)
magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure], 
                                       _arXc=array('d', [0.0]), 
                                       _arYc=array('d', [0.0]), 
                                       _arZc=array('d', [0.0]))

mesh = SRWLRadMesh(_eStart=9999.998465105169,
                   _eFin  =9999.998465105169,
                   _ne    =1,
                   _xStart=-0.00025,
                   _xFin  =0.00025,
                   _nx    =256,
                   _yStart=-0.00025,
                   _yFin  =0.00025,
                   _ny    =256,
                   _zStart=44.09)

stk = SRWLStokes()
stk.allocate(1,256,256)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam

initial_mesh = deepcopy(wfr.mesh)
srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1,0.01,0.0,0.0,50000,1,0.0])

mesh0 = deepcopy(wfr.mesh)
arI = array('f', [0]*mesh0.nx*mesh0.ny)
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
arIx = array('f', [0]*mesh0.nx)
srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
arIy = array('f', [0]*mesh0.ny)
srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot2d1d (arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []

oe_0=SRWLOptA(_shape='c',
               _ap_or_ob='a',
               _Dx=0.0005,
               _Dy=0.0005,
               _x=0.0,
               _y=0.0)

pp_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_0)
srw_pp_array.append(pp_oe_0)

oe_1=srwl_opt_setup_CRL(_foc_plane=3,
                _delta=3.408273456573241e-06,
                _atten_len=0.00836915564073516,
                _shape=1,
                _apert_h=0.0005,
                _apert_v=0.0005,
                _r_min=5e-05,
                _n=15,
                _wall_thick=3e-05,
                _xc=0.0,
                _yc=0.0,
                _void_cen_rad=None,
                _e_start=0.0,
                _e_fin=0.0,
                _nx=1001,
                _ny=1001)

pp_oe_1 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_1)
srw_pp_array.append(pp_oe_1)

drift_before_oe_2 = SRWLOptD(0.4945)
pp_drift_before_oe_2 = [0,0,1.0,2,0,4.0,2,4.0,2,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_2)
srw_pp_array.append(pp_drift_before_oe_2)



####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)
srwl.PropagElecField(wfr, optBL)

mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny)
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0)
arI1x = array('f', [0]*mesh1.nx)
srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, 0)
arI1y = array('f', [0]*mesh1.ny)
srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, 0, 0)

#%%
#画二维图
plt.figure(2)
I2D = numpy.reshape(arI1,(2048,2048))
x0 = numpy.linspace(1e6*mesh1.xStart, 1e6*mesh1.xFin, mesh1.nx)
y0 = numpy.linspace(1e6*mesh1.yStart, 1e6*mesh1.yFin, mesh1.ny)
im =plt.imshow(I2D, interpolation='bicubic', cmap='jet',
               extent=[min(y0),max(y0),min(x0),max(x0)], aspect='auto')#,vmin=vmin, vmax=vmax)
plt.axvline(x=0.0,color = 'r', linestyle = '--')
plt.xlim([-2, 2])
plt.ylim([-2, 2])
plt.xlabel('x(um)')
plt.ylabel('y(um)')
plt.colorbar(im, shrink=0.5)
plt.show()  
plt.savefig('A1.png',dpi = 300)
#%%
#选择中间一条线输出
plt.figure(3)
plt.plot(x0,I2D[len(y0)//2,:])
plt.xlabel('position(um)')
plt.ylabel('Intensity(a.u.)')  
plt.legend()
plt.show()
plt.savefig('ideal_yz_Section_1D.png',dpi = 300)


