import os, sys
sys.path.append(r".//SRW_base//")  #srw base 文件位置

try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

if not srwl_uti_proc_is_master(): exit()

import pickle
import matplotlib.pyplot as plt
from scipy import interpolate


def caclcuation_SR(dx0):
    ####################################################
    # LIGHT SOURCE

    part_beam = SRWLPartBeam()
    part_beam.Iavg               = 0.2
    part_beam.partStatMom1.x     = 0.0
    part_beam.partStatMom1.y     = 0.0
    part_beam.partStatMom1.z     = -2.555
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
                    _xStart=-0.00045,
                    _xFin  =0.00045,
                    _nx    =256,
                    _yStart=-0.00024,
                    _yFin  =0.00024,
                    _ny    =256,
                    _zStart=37.0)

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

    oe_0 = SRWLOptMirEl(_size_tang=0.54,
                        _size_sag=0.01,
                        _p=37.0,
                        _q=100000000.0,
                        _ang_graz=0.0016999506914423595,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1,
                        _nvx=0,
                        _nvy=0.9999985550841713,
                        _nvz=-0.0016999498726803933,
                        _tvx=0,
                        _tvy=-0.0016999498726803933,
                        _x=0.0,
                        _y=0.0)
    oe_0.set_dim_sim_meth(_size_tang=0.54,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1)
    oe_0.set_orient(_nvx=0,
                    _nvy=0.9999985550841713,
                    _nvz=-0.0016999498726803933,
                    _tvx=0,
                    _tvy=-0.0016999498726803933,
                    _x=0.0,
                    _y=0.0)


    pp_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_0)
    srw_pp_array.append(pp_oe_0)

    drift_after_oe_0 = SRWLOptD(2.0)
    pp_drift_after_oe_0 = [0,0,1.0,1,0,2.0,1.0,2.0,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_0)
    srw_pp_array.append(pp_drift_after_oe_0)

    oe_1=SRWLOptCryst(_d_sp=3.1355,
                            _psi0r=-9.7631e-06,
                            _psi0i=1.4871e-07,
                            _psi_hr=5.1546e-06,
                            _psi_hi=1.0368e-07,
                            _psi_hbr=4.7516e-06,
                            _psi_hbi=9.5183e-08,
                            _tc=0.01,
                            _ang_as=0.0,
                            _nvx=0,
                            _nvy=0.980260823963851,
                            _nvz=-0.19770866698683703,
                            _tvx=0,
                            _tvy=-0.19770866698683703,
                            _uc=1)

    pp_oe_1 = [0,0,1.0,2,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_1)
    srw_pp_array.append(pp_oe_1)

    drift_after_oe_1 = SRWLOptD(0.02)
    pp_drift_after_oe_1 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_1)
    srw_pp_array.append(pp_drift_after_oe_1)

    oe_2=SRWLOptCryst(_d_sp=3.1355,
                            _psi0r=-9.7631e-06,
                            _psi0i=1.4871e-07,
                            _psi_hr=5.1546e-06,
                            _psi_hi=1.0368e-07,
                            _psi_hbr=4.7516e-06,
                            _psi_hbi=9.5183e-08,
                            _tc=0.01,
                            _ang_as=0.0,
                            _nvx=0,
                            _nvy=-0.980260823963851,
                            _nvz=-0.19770866698683703,
                            _tvx=0,
                            _tvy=-0.19770866698683703,
                            _uc=1)

    pp_oe_2 = [0,0,1.0,2,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_2)
    srw_pp_array.append(pp_oe_2)

    drift_after_oe_2 = SRWLOptD(2.48)
    pp_drift_after_oe_2 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_2)
    srw_pp_array.append(pp_drift_after_oe_2)

    oe_3=SRWLOptA(_shape='r',
                _ap_or_ob='a',
                _Dx=0.0008,
                _Dy=0.000747,
                _x=0.0,
                _y=0.0)

    pp_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_3)
    srw_pp_array.append(pp_oe_3)

    drift_after_oe_3 = SRWLOptD(1.0)
    pp_drift_after_oe_3 = [0,0,1.0,2,0,1.0,1.0,1.0,2.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_3)
    srw_pp_array.append(pp_drift_after_oe_3)

    oe_4 = SRWLOptMirPl(_size_tang=0.45,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1,
                        _nvx=0,
                        _nvy=0.9999985550841713,
                        _nvz=-0.0016999498726803933,
                        _tvx=0,
                        _tvy=-0.0016999498726803933,
                        _x=0.0,
                        _y=0.0)
    oe_4.set_dim_sim_meth(_size_tang=0.45,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1)
    oe_4.set_orient(_nvx=0,
                    _nvy=0.9999985550841713,
                    _nvz=-0.0016999498726803933,
                    _tvx=0,
                    _tvy=-0.0016999498726803933,
                    _x=0.0,
                    _y=0.0)


    pp_oe_4 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_4)
    srw_pp_array.append(pp_oe_4)

    drift_after_oe_4 = SRWLOptD(1.176)
    pp_drift_after_oe_4 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_4)
    srw_pp_array.append(pp_drift_after_oe_4)

    oe_5 = SRWLOptMirPl(_size_tang=0.45,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1,
                        _nvx=0,
                        _nvy=-0.9999985550841713,
                        _nvz=-0.0016999498726803933,
                        _tvx=0,
                        _tvy=-0.0016999498726803933,
                        _x=0.0,
                        _y=0.0)
    oe_5.set_dim_sim_meth(_size_tang=0.45,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1)
    oe_5.set_orient(_nvx=0,
                    _nvy=-0.9999985550841713,
                    _nvz=-0.0016999498726803933,
                    _tvx=0,
                    _tvy=-0.0016999498726803933,
                    _x=0.0,
                    _y=0.0)


    pp_oe_5 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_5)
    srw_pp_array.append(pp_oe_5)

    drift_after_oe_5 = SRWLOptD(1.4)
    pp_drift_after_oe_5 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_5)
    srw_pp_array.append(pp_drift_after_oe_5)

    oe_6 = SRWLOptMirEl(_size_tang=0.6,
                        _size_sag=0.01,
                        _p=100000000.0,
                        _q=4.5,
                        _ang_graz=0.0016999506914423595,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1,
                        _nvx=0,
                        _nvy=-0.9999985550841713,
                        _nvz=-0.0016999498726803933,
                        _tvx=0,
                        _tvy=-0.0016999498726803933,
                        _x=0.0,
                        _y=0.0)
    oe_6.set_dim_sim_meth(_size_tang=0.6,
                        _size_sag=0.01,
                        _ap_shape='r',
                        _sim_meth=2,
                        _treat_in_out=1)
    oe_6.set_orient(_nvx=0,
                    _nvy=-0.9999985550841713,
                    _nvz=-0.0016999498726803933,
                    _tvx=0,
                    _tvy=-0.0016999498726803933,
                    _x=0.0,
                    _y=0.0)


    pp_oe_6 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_6)
    srw_pp_array.append(pp_oe_6)

    drift_after_oe_6 = SRWLOptD(4.5+dx0)
    pp_drift_after_oe_6 = [0,0,1.0,1,0,1.0,1.0,2.0,2.0,0,0.0,0.0]

    srw_oe_array.append(drift_after_oe_6)
    srw_pp_array.append(pp_drift_after_oe_6)



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

    x0 = numpy.linspace(1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx)
    y0 = numpy.linspace(1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny)
    return arI1x, x0, arI1y, y0

def main():
    Position = numpy.linspace(-1000, 1000, 31)*1e-3
    image2Dx = []
    image2Dy = []
    x1D0 = numpy.linspace(-0.6,0.6,200)
    y1D0 = numpy.linspace(-0.1,0.1,1000)
    for i, dx0 in enumerate(Position):
        I1Dx, x1D, I1Dy, y1D = caclcuation_SR(dx0)
        fx = interpolate.interp1d(x1D, I1Dx, kind='cubic')
        fy = interpolate.interp1d(y1D, I1Dy, kind='cubic')
        I1Dx0 = fx(x1D0)
        I1Dy0 = fy(y1D0)
        image2Dx.append(I1Dx0)
        image2Dy.append(I1Dy0)
        print(i)
    dump = [x1D0, y1D0,Position, image2Dx, image2Dy]   
    pickleName = '2D_section.pickle'
    with open(pickleName, 'wb') as f:
        pickle.dump(dump, f, protocol=2)

if __name__ == '__main__':
    main()
