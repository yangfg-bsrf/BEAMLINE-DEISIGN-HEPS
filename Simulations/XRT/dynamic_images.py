#从保存的图片中读取图片，并保存动态图片
#这里的文件名是与主追迹软件中的plot相关联；因此调用plot.savename
import imageio
def dynamic_images(Plots):
    gif_images = []    
    for plot in Plots:
        filename = plot.saveName
        gif_images.append(imageio.imread(filename)) 

    imageio.mimsave("test.gif",gif_images,fps=1.5) 
    print("dynamic images : Done!")   