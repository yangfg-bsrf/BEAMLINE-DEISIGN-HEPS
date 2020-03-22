
#将误差数据文件放在上一层目录下的surface error_DATA文件夹下
#访问途径
Error_path = '../surface error_DATA/'
""" 第二光学元件：准直镜@37m """  
if False:
    p_WCM = 37000
    theta_WCM = 1.7e-3
    oe_error_WCM = 0
    Surface_name=  Error_path + 'VCM45_'
    kwargs_OE = dict(
                name='WCM', 
                center=[0, p_WCM, 0],
                pitch= theta_WCM,  
                material = Si,
                limPhysX=[-1, 1], limPhysY=[-275,275],
                p=p_WCM, q=5000000000,
                isCylindrical=True)    
    if oe_error_WCM == 0:        
        beamLine.WMirr = raycing.oes.EllipticalMirrorParam(bl=beamLine,**kwargs_OE)
    else:            
        kwargs_OE['get_distorted_surface'] ='error'
        kwargs_OE['fname1'] =Surface_name+'X.txt'
        kwargs_OE['fname2'] =Surface_name+'Y.txt'
        kwargs_OE['fname3'] =Surface_name+'Z.txt'       
        beamLine.WMirr = Distorted.EllipMirrorDistorted(beamLine, **kwargs_OE)  