#from data_utils import affinity_utils
import math
import numpy as np
from multiprocessing import Pool
from skimage.morphology import skeletonize
#import data_utils.graph_utils as graph_utils
import data_utils.sknw as sknw
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
import numpy as np
from numpy import ma
import skimage
import os
import cv2
import torch
from torch.utils.data import Dataset
from osgeo import gdal
from collections import OrderedDict
import os.path as osp
from pathlib import Path
import time
import argparse
from multiprocessing import Process
from osgeo import gdal
from torch.nn.modules.utils import _pair
import multiprocessing
class RSInferenceDataset(Dataset):
    def __init__(
            self,
            rs_file,
            patch_size=(1024, 1024),
            slide_step=(1024, 1024),
            pad=256,
    ):
        super().__init__()
        self.rs_file = rs_file
        self.patch_size = _pair(patch_size)
        self.slide_step = _pair(slide_step)

        # get data info

        ds = gdal.Open(rs_file)
        self.data_info = self._get_data_info(ds)
        self.ids = self._get_patch_ids()
        ds = None
        self.pad=pad

    def __getitem__(self, idx):
        img,center,center_ = self._read_patch(idx)
        ids=np.array((self.ids[idx],center,center_))
       # print(ids)
        return img, ids

    def __len__(self):
        return len(self.ids)

    def _get_data_info(self, src):
        return {
            'width': src.RasterXSize,
            'height': src.RasterYSize,
            'driver': src.GetDriver().ShortName,
            'dtype': gdal.GetDataTypeName(src.GetRasterBand(1).DataType),
            'bands': src.RasterCount,
            'proj': src.GetProjection(),
            'geotransform': src.GetGeoTransform(),
        }

    def _get_patch_ids(self):
        left, top = 0, 0
        width, height = self.data_info['width'], self.data_info['height']
        if width%self.patch_size[0]!=0:
            width=(width//self.patch_size[0]+1)*self.patch_size[0]
        if height%self.patch_size[1]!=0:
            height=(height//self.patch_size[1]+1)*self.patch_size[1]
        left_top_xy = []  # left-top corner coordinates (xmin, ymin)
        while left < width:
            top = 0
            while top < height:
                left_top_xy.append((left, top))
                if top + self.patch_size[1] >= height:
                    break
                else:
                    top += self.slide_step[1]

            if left + self.patch_size[0] >= width:
                break
            else:
                left += self.slide_step[0]

        return left_top_xy# 

    def _read_patch(self, idx):
        xmin, ymin = self.ids[idx]
        width, height = self.data_info['width'], self.data_info['height']
        center_xsize,center_ysize=self.patch_size[0],self.patch_size[1]
        center_x,center_y=self.pad,self.pad
        xsize,ysize=self.patch_size[0]+2*self.pad,self.patch_size[1]+2*self.pad
        x_min,y_min=xmin-self.pad,ymin-self.pad
        left=0#x
        top=0#left-top corner coordinates (xmin, ymin)
        if x_min<0:
            xsize=xsize+x_min
            center_x=center_x+x_min
            x_min=0
        elif x_min+self.patch_size[0]+2*self.pad>width:
            xsize=width-x_min
            center_xsize=min(xsize-self.pad,self.patch_size[0])
        if y_min<0:
            ysize=ysize+y_min
            center_y=center_y+y_min
            y_min=0
        elif y_min+self.patch_size[1]+2*self.pad>height:
            ysize=height-y_min
            center_ysize=min(ysize-self.pad,self.patch_size[1])
            
        # to use multi-processing
        ds = gdal.Open(self.rs_file)
       
        band = ds.GetRasterBand(1)
        
        img = band.ReadAsArray(
            xoff=x_min,
            yoff=y_min,
            win_xsize=xsize,
            win_ysize=ysize,
        )
        img[img!=3]=0
        img[img==3]=1
        center=(center_xsize,center_ysize)
        center_=(center_x,center_y)
        ds = None
        return img,center,center_

    @property
    def width(self):
        return self.data_info['width']

    @property
    def height(self):
        return self.data_info['height']

def distance(x,y):
    return np.sqrt(np.sum(np.square([x[0]-y[0],x[1]-y[1]])))
def patch_regular(gt,tau,thickness):
    #gt=tifffile.imread(img_path)
    ske = skeletonize(gt).astype(np.uint16)
    graph = sknw.build_sknw(ske)
    points=[]
    nodes=set()
    # draw edges by pts
    for (s,e) in graph.edges():
        ps = graph[s][e]['pts']
        p1=[float(ps[0,0]),float(ps[0,1])]
        p2=[float(ps[-1,0]),float(ps[-1,1])]
        nodes.add(str(p1))
        nodes.add(str(p2))
        points.append({str(p1),str(p2)})
        for i in range(0,len(ps)-1):
            cv2.line(gt,(int(ps[i,1]),int(ps[i,0])), (int(ps[i+1,1]),int(ps[i+1,0])), 1,thickness=thickness)
    ps=[eval(i) for i in list(nodes)]
    for num in range(len(ps)):
        mindis=float("inf")
        for other in range(len(ps)):
            if other!=num and {str(ps[num]),str(ps[other])} not in points:
                dis= distance(ps[num],ps[other])
                if dis<mindis:
                    mindis=dis
                    mindis_point=other
        if mindis<tau:
            cv2.line(gt,(int(ps[num][1]),int(ps[num][0])), (int(ps[mindis_point][1]),int(ps[mindis_point][0])), 1,thickness=thickness)
    return gt
#      if (abs(ps[num]-ps[num+1])).sum()<50:
#          n+=1
#          plt.plot([ps[num,1],ps[num+1,1]], [ps[num,0],ps[num+1,0]], 'black',linewidth=2)
    
        
    # draw node by o
#    nodes = graph.nodes()
#    ps = np.array([nodes[i]['o'] for i in nodes])
    
    
    # title and show
#    name=img_path.split('/')[-1][:-4]
#    out_path=os.path.join("/data1/qiuxinling/figure/",name+'.tif')
#    tifffile.imsave(out_path,gt)
    #    all_data = []
    #    for k, v in data:
    #        for val in v:
    #            all_data.append((k, val))
    #        all_data.append((np.NaN,np.NaN))
    #    df = pd.DataFrame(all_data, columns=['ImageId', 'WKT_Pix'])
    #    df.to_csv('xk1_gray.csv', index=False)
def regular(img, offset,remove_obj_holes,tau,thickness,q):
    offset=np.squeeze(offset)
    center_xsize,center_ysize=offset[1]
    center_x, center_y = offset[2]
    offset=offset[0]
    x_offset, y_offset = offset
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (remove_obj_holes, remove_obj_holes))
    img = cv2.morphologyEx(img, cv2.MORPH_CLOSE, kernel)
    img_regular=patch_regular(img,tau,thickness)
    #img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

        # pred = result['out']  # nB,nC,H,W
    img_regular = img_regular[center_y:(center_y+center_ysize), center_x:(center_x+center_xsize)]
    q.put(img_regular)
                                        
                                        
                                        
def large_image_regular(image_file,out_dir,patch_size=1024,remove_obj_holes=11,tau=50,thickness=10):
    suffixs = ('.tif', '.tiff', '.img')
    if osp.isfile(image_file) and image_file.endswith(suffixs):
        image_list = [image_file]
    else:
        image_list = [
            osp.join(args.image_file, f)
            for f in scandir(args.image_file, suffix=suffixs)
        ]
    for img_file in image_list:
        dataset = RSInferenceDataset(img_file,
                                     patch_size=patch_size,
                                     slide_step=patch_size,
                                     pad=patch_size//4)
        basename = Path(img_file).stem
        out_file = osp.join(out_dir, f'{basename}_regular.tif')
        driver = gdal.GetDriverByName('GTiff')
        src_ds = gdal.Open(img_file)
        out_raster = driver.Create(out_file, dataset.width, dataset.height,
                                   1, gdal.GDT_Byte)
        gt = src_ds.GetGeoTransform()
        if gt is not None:
            out_raster.SetGeoTransform(gt)
        out_raster.SetProjection(src_ds.GetProjection())
        src_ds = None
        out_band = out_raster.GetRasterBand(1)
        pbar = tqdm(dataset)
        q = multiprocessing.Queue()
        for img, offset in pbar:
            p=multiprocessing.Process(target=regular, args=(img, offset,remove_obj_holes,tau,thickness,q))
            p.start()
            x_offset, y_offset=np.squeeze(offset)[0]
            out_band.WriteArray(q.get(),xoff=x_offset.item(),
                                        yoff=y_offset.item())
        out_raster = None
def parse_args():
    parser = argparse.ArgumentParser(description='road regular')

    parser.add_argument('image_file', help='input file path or directory')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('--patch_size',
                        type=int,
                        default=1024,
                        help='patch size')
    parser.add_argument('--remove_obj_holes',
                        type=int,
                        default=11,
                        help='remove small objects and holes')
    parser.add_argument('--tau',
                        type=int,
                        default=100,
                        help='connect value')
    parser.add_argument('--thickness',
                        type=int,
                        default=10,
                        help='thickness value')         
                        
    args = parser.parse_args()
    return args
if __name__ == "__main__":
    args = parse_args()
    large_image_regular(args.image_file,args.out_dir,remove_obj_holes=args.remove_obj_holes,tau=args.tau,thickness=args.thickness)


    
