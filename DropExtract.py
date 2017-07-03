"""

@author: zhaojie


"""


from __future__ import division
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
import numpy as np
import tifffile as T
from skimage import color
from skimage import filters
import skimage
import scipy.ndimage



def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue
        if 'df1_x' in str(globals()[var]): continue
        if 'df1' in str(globals()[var]): continue
        if 'size1' in str(globals()[var]): continue
        if 'df2_x' in str(globals()[var]): continue
        if 'df2' in str(globals()[var]): continue
        if 'size2' in str(globals()[var]): continue
        if 'df3_x' in str(globals()[var]): continue
        if 'df3' in str(globals()[var]): continue
        if 'size3' in str(globals()[var]): continue
        if 'df4_x' in str(globals()[var]): continue
        if 'df4' in str(globals()[var]): continue
        if 'size4' in str(globals()[var]): continue
        if 'df5_x' in str(globals()[var]): continue
        if 'df5' in str(globals()[var]): continue
        if 'size5' in str(globals()[var]): continue
        
        
        del globals()[var]
        
        

def cut_ROI(img,a,b,c,d):
    im1 = color.rgb2gray(img)
    im2 = im1[a:b,c:d]
    return im2

    
#first set of image files    
#read image stacks and cut

a=166 ## upper edge position
b=400 ## lower edge position
c=30  ## right starting point postion
d=680 ## ending edge postion
a1=189
b1=423
c1=21
d1=671


img=T.imread( u'/Users/username/Desktop/Droplets separation/3.44-0.16-45/f-35.tif')
img_m=T.imread(u'/Users/username/Desktop/Droplets separation/3.44-0.16-45/fm-35.tif')

img_cut=np.ones((img.shape[0],b-a,650))  # 650 is the length of the ROI
img_m_cut=np.ones((img.shape[0],b1-a1,650))
for i1 in range(img.shape[0]):
    img_cut[i1]=cut_ROI(img[i1],a,b,c,d)
for i2 in range(img_m.shape[0]):
    img_m_cut[i2]=cut_ROI(img_m[i2],a1,b1,c1,d1)    
    
#get the size and center of the droplets
droplet_size=[]
droplet_centroid=[]
droplet_centroid_m=[]

for j in range(img_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_cut[j])
    mask=threshold<img_cut[j]
    label_img=skimage.measure.label(mask)
    props1=skimage.measure.regionprops(label_img)
    for k in range(len(props1)):
        if props1[k].equivalent_diameter*1023/(b-a)>30 and  props1[k].equivalent_diameter*1023/(b-a)<70: #convert scale from pixel unit to actual unit in micrometer exclude outliners
            droplet_size.append(props1[k].equivalent_diameter*1023/234)
            droplet_centroid.append(props1[k].centroid)
        else:
            continue

for h in range(img_m_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_m_cut[h])
    mask=threshold<img_m_cut[h]
    label_img_m=skimage.measure.label(mask)
    props2=skimage.measure.regionprops(label_img_m)
    for l in range(len(props2)):
        if props2[l].equivalent_diameter*1023/(b1-a1)>30 and  props2[l].equivalent_diameter*1023/(b1-a1)<70:
            droplet_centroid_m.append(props2[l].centroid)
        else:
            continue
droplet_size_arr=np.array(droplet_size)
plt.figure(figsize=(12,10))
plt.hist(droplet_size_arr,bins=10)
plt.xlabel('Droplet_Equivalent_Diameter (um)')
plt.ylabel('Count')
plt.show()


droplet_centroid_arr=np.array(droplet_centroid)
droplet_centroid_m_arr=np.array(droplet_centroid_m)
stack=np.max(img_cut,axis=0)
stack_m=np.max(img_m_cut,axis=0)
stack_temp=np.zeros((2,img_cut.shape[1],img_cut.shape[2]))
stack_temp[0]=stack
stack_temp[1]=stack_m
stack_all=np.max(stack_temp,axis=0)
stack_all_zoom=scipy.ndimage.zoom(stack_all,1023/(b-a))


plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
plt.gca().invert_yaxis()
plt.imshow(stack_all_zoom,'gray')
plt.plot(droplet_centroid_arr[:,1]*1023/(b-a),droplet_centroid_arr[:,0]*1023/(b-a),'ob',markersize=2)
plt.plot(droplet_centroid_m_arr[:,1]*1023/(b1-a1),droplet_centroid_m_arr[:,0]*1023/(b1-a1),'or',markersize=2)
plt.axis([0,650*1023/(b-a),1023,0])  #234 and 650 is the width and length of the cutting region (ROI) of the image

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)

plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Y Position (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
plt.savefig(u'/Users/username/Desktop/1.tiff',bbox_inches='tight')

#first set of droplets
size1=np.mean(droplet_size_arr)
dis=np.zeros((1,len(droplet_centroid)))
df1=np.zeros((1,len(droplet_centroid_m)))
df1_x=np.zeros((1,len(droplet_centroid_m)))
for m in range(len(droplet_centroid_m)):
    for n in range(len(droplet_centroid)):
        dis[0][n]=abs(droplet_centroid[n][1]-droplet_centroid_m[m][1])
    N=dis.argmin()
    df1[0][m]=abs(droplet_centroid[N][0]-droplet_centroid_m[m][0])
    df1_x[0][m]=droplet_centroid_m[m][1]    

plt.figure(figsize=(12,10))
plt.plot(df1_x*1023/(b-a),df1*1023/(b-a),'ob')
plt.xlabel('X Position (\u03bcm)')
plt.ylabel('Deflection (\u03bcm)')
plt.show()


#second set of image files
clear_all

a=162 ## upper edge position
b=397 ## lower edge position
c=18 ## right starting point postion
d=668 ## ending edge postion
a1=149
b1=384
c1=20
d1=670
    
#read image stacks and cut
img=T.imread( u'/Users/username/Desktop/3.12-0.48-45/f-50.tif')
img_m=T.imread(u'/Users/username/Desktop/3.12-0.48-45/fm-50.tif')

img_cut=np.ones((img.shape[0],b-a,650)) 
img_m_cut=np.ones((img.shape[0],b1-a1,650))
for i1 in range(img.shape[0]):
    img_cut[i1]=cut_ROI(img[i1],a,b,c,d)
for i2 in range(img_m.shape[0]):
    img_m_cut[i2]=cut_ROI(img_m[i2],a1,b1,c1,d1)    
    
#get the size and center of the droplets
droplet_size=[]
droplet_centroid=[]
droplet_centroid_m=[]

for j in range(img_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_cut[j])
    mask=threshold<img_cut[j]
    label_img=skimage.measure.label(mask)
    props1=skimage.measure.regionprops(label_img)
    for k in range(len(props1)):
        if props1[k].equivalent_diameter*1023/(b-a)>40 and  props1[k].equivalent_diameter*1023/(b-a)<110:
            droplet_size.append(props1[k].equivalent_diameter*1023/234)
            droplet_centroid.append(props1[k].centroid)
        else:
            continue

for h in range(img_m_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_m_cut[h])
    mask=threshold<img_m_cut[h]
    label_img_m=skimage.measure.label(mask)
    props2=skimage.measure.regionprops(label_img_m)
    for l in range(len(props2)):
        if props2[l].equivalent_diameter*1023/(b1-a1)>40 and  props2[l].equivalent_diameter*1023/(b1-a1)<110:
            droplet_centroid_m.append(props2[l].centroid)
        else:
            continue
droplet_size_arr=np.array(droplet_size)
plt.figure(figsize=(12,10))
plt.hist(droplet_size_arr,bins=10)
plt.xlabel('Droplet_Equivalent_Diameter (um)')
plt.ylabel('Count')
plt.show()


droplet_centroid_arr=np.array(droplet_centroid)
droplet_centroid_m_arr=np.array(droplet_centroid_m)
stack=np.max(img_cut,axis=0)
stack_m=np.max(img_m_cut,axis=0)
stack_temp=np.zeros((2,img_cut.shape[1],img_cut.shape[2]))
stack_temp[0]=stack
stack_temp[1]=stack_m
stack_all=np.max(stack_temp,axis=0)
stack_all_zoom=scipy.ndimage.zoom(stack_all,1023/(b-a))



plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
plt.gca().invert_yaxis()
plt.imshow(stack_all_zoom,'gray')
plt.plot(droplet_centroid_arr[:,1]*1023/(b-a),droplet_centroid_arr[:,0]*1023/(b-a),'ob',markersize=2)
plt.plot(droplet_centroid_m_arr[:,1]*1023/(b1-a1),droplet_centroid_m_arr[:,0]*1023/(b1-a1),'or',markersize=2)
plt.axis([0,650*1023/(b-a),1023,0])  #234 and 650 is the width and length of the cutting region (ROI) of the image, this is only in pixel unit but it will slight change according to the actual cut of the image

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)
plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Y Position (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
plt.savefig(u'/Users/username/Desktop/2.tiff',bbox_inches='tight')

#second set of droplets
size2=np.mean(droplet_size_arr)
dis=np.zeros((1,len(droplet_centroid)))
df2=np.zeros((1,len(droplet_centroid_m)))
df2_x=np.zeros((1,len(droplet_centroid_m)))
for m in range(len(droplet_centroid_m)):
    for n in range(len(droplet_centroid)):
        dis[0][n]=abs(droplet_centroid[n][1]-droplet_centroid_m[m][1])
    N=dis.argmin()
    df2[0][m]=abs(droplet_centroid[N][0]-droplet_centroid_m[m][0])
    df2_x[0][m]=droplet_centroid_m[m][1]    

plt.figure(figsize=(12,10))
plt.plot(df2_x*1023/(b-a),df2*1023/(b-a),'ob')
plt.xlabel('X Position (u03bcm)')
plt.ylabel('Deflection (u03bcm)')
plt.show()


#third set of image files
clear_all
 
a=163 ##upper edge position
b=399 ##lower edge position
c=18  ##starting point position
d=668 ## ending edge postion
a1=177
b1=413
c1=20
d1=670
    
#read image stacks and cut
img=T.imread( u'/Users/username/Desktop/2.8-0.8-45/f-50.tif')
img_m=T.imread(u'/Users/username/Desktop/2.8-0.8-45/fm-50.tif')

img_cut=np.ones((img.shape[0],b-a,650)) 
img_m_cut=np.ones((img.shape[0],b1-a1,650))
for i1 in range(img.shape[0]):
    img_cut[i1]=cut_ROI(img[i1],a,b,c,d)
for i2 in range(img_m.shape[0]):
    img_m_cut[i2]=cut_ROI(img_m[i2],a1,b1,c1,d1)    
    
#get the size and center of the droplets
droplet_size=[]
droplet_centroid=[]
droplet_centroid_m=[]

for j in range(img_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_cut[j])
    mask=threshold<img_cut[j]
    label_img=skimage.measure.label(mask)
    props1=skimage.measure.regionprops(label_img)
    for k in range(len(props1)):
        if props1[k].equivalent_diameter*1023/(b-a)>50 and  props1[k].equivalent_diameter*1023/(b-a)<140:
            droplet_size.append(props1[k].equivalent_diameter*1023/234)
            droplet_centroid.append(props1[k].centroid)
        else:
            continue

for h in range(img_m_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_m_cut[h])
    mask=threshold<img_m_cut[h]
    label_img_m=skimage.measure.label(mask)
    props2=skimage.measure.regionprops(label_img_m)
    for l in range(len(props2)):
        if props2[l].equivalent_diameter*1023/(b1-a1)>50 and  props2[l].equivalent_diameter*1023/(b1-a1)<140:
            droplet_centroid_m.append(props2[l].centroid)
        else:
            continue
droplet_size_arr=np.array(droplet_size)
plt.figure(figsize=(12,10))
plt.hist(droplet_size_arr,bins=10)
plt.xlabel('Droplet_Equivalent_Diameter (\u03bcm)')
plt.ylabel('Count')
plt.show()


droplet_centroid_arr=np.array(droplet_centroid)
droplet_centroid_m_arr=np.array(droplet_centroid_m)
stack=np.max(img_cut,axis=0)
stack_m=np.max(img_m_cut,axis=0)
stack_temp=np.zeros((2,img_cut.shape[1],img_cut.shape[2]))
stack_temp[0]=stack
stack_temp[1]=stack_m
stack_all=np.max(stack_temp,axis=0)
stack_all_zoom=scipy.ndimage.zoom(stack_all,1023/(b-a))


plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
plt.gca().invert_yaxis()
plt.imshow(stack_all_zoom,'gray')
plt.plot(droplet_centroid_arr[:,1]*1023/(b-a),droplet_centroid_arr[:,0]*1023/(b-a),'ob',markersize=2)
plt.plot(droplet_centroid_m_arr[:,1]*1023/(b1-a1),droplet_centroid_m_arr[:,0]*1023/(b1-a1),'or',markersize=2)
plt.axis([0,650*1023/(b-a),1023,0])  

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)
plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Y Position (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
plt.savefig(u'/Users/username/Desktop/3.tiff',bbox_inches='tight')

#third set of droplets
size3=np.mean(droplet_size_arr)
dis=np.zeros((1,len(droplet_centroid)))
df3=np.zeros((1,len(droplet_centroid_m)))
df3_x=np.zeros((1,len(droplet_centroid_m)))
for m in range(len(droplet_centroid_m)):
    for n in range(len(droplet_centroid)):
        dis[0][n]=abs(droplet_centroid[n][1]-droplet_centroid_m[m][1])
    N=dis.argmin()
    df3[0][m]=abs(droplet_centroid[N][0]-droplet_centroid_m[m][0])
    df3_x[0][m]=droplet_centroid_m[m][1]    

plt.figure(figsize=(12,10))
plt.plot(df3_x*1023/(b-a),df3*1023/(b-a),'ob')
plt.xlabel('X Position (\u03bcm)')
plt.ylabel('Deflection (\u03bcm)')
plt.show()


#fourth set of image files
clear_all
 
a=161 ##upper edge position
b=396 ##lower edge position
c=23  ##starting point position
d=673 ## ending edge postion
a1=144
b1=379
c1=21
d1=671
    
#read image stacks and cut
img=T.imread( u'/Users/username/Desktop/2.48-1.12-45/f-50.tif')
img_m=T.imread(u'/Users/username/Desktop/2.48-1.12-45/fm-50.tif')

img_cut=np.ones((img.shape[0],b-a,650)) 
img_m_cut=np.ones((img.shape[0],b1-a1,650))
for i1 in range(img.shape[0]):
    img_cut[i1]=cut_ROI(img[i1],a,b,c,d)
for i2 in range(img_m.shape[0]):
    img_m_cut[i2]=cut_ROI(img_m[i2],a1,b1,c1,d1)    
    
#get the size and center of the droplets
droplet_size=[]
droplet_centroid=[]
droplet_centroid_m=[]

for j in range(img_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_cut[j])
    mask=threshold<img_cut[j]
    label_img=skimage.measure.label(mask)
    props1=skimage.measure.regionprops(label_img)
    for k in range(len(props1)):
        if props1[k].equivalent_diameter*1023/(b-a)>60 and  props1[k].equivalent_diameter*1023/(b-a)<150:
            droplet_size.append(props1[k].equivalent_diameter*1023/234)
            droplet_centroid.append(props1[k].centroid)
        else:
            continue

for h in range(img_m_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_m_cut[h])
    mask=threshold<img_m_cut[h]
    label_img_m=skimage.measure.label(mask)
    props2=skimage.measure.regionprops(label_img_m)
    for l in range(len(props2)):
        if props2[l].equivalent_diameter*1023/(b1-a1)>60 and  props2[l].equivalent_diameter*1023/(b1-a1)<150:
            droplet_centroid_m.append(props2[l].centroid)
        else:
            continue
droplet_size_arr=np.array(droplet_size)
plt.figure(figsize=(12,10))
plt.hist(droplet_size_arr,bins=10)
plt.xlabel('Droplet_Equivalent_Diameter (\u03bcm)')
plt.ylabel('Count')
plt.show()


droplet_centroid_arr=np.array(droplet_centroid)
droplet_centroid_m_arr=np.array(droplet_centroid_m)
stack=np.max(img_cut,axis=0)
stack_m=np.max(img_m_cut,axis=0)
stack_temp=np.zeros((2,img_cut.shape[1],img_cut.shape[2]))
stack_temp[0]=stack
stack_temp[1]=stack_m
stack_all=np.max(stack_temp,axis=0)
stack_all_zoom=scipy.ndimage.zoom(stack_all,1023/(b-a))


plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
plt.gca().invert_yaxis()
plt.imshow(stack_all_zoom,'gray')
plt.plot(droplet_centroid_arr[:,1]*1023/(b-a),droplet_centroid_arr[:,0]*1023/(b-a),'ob',markersize=2)
plt.plot(droplet_centroid_m_arr[:,1]*1023/(b1-a1),droplet_centroid_m_arr[:,0]*1023/(b1-a1),'or',markersize=2)
plt.axis([0,650*1023/(b-a),1023,0]) 

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)
plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Y Position (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
plt.savefig(u'/Users/username/Desktop/4.tiff',bbox_inches='tight')

#fouth set of droplets
size4=np.mean(droplet_size_arr)
dis=np.zeros((1,len(droplet_centroid)))
df4=np.zeros((1,len(droplet_centroid_m)))
df4_x=np.zeros((1,len(droplet_centroid_m)))
for m in range(len(droplet_centroid_m)):
    for n in range(len(droplet_centroid)):
        dis[0][n]=abs(droplet_centroid[n][1]-droplet_centroid_m[m][1])
    N=dis.argmin()
    df4[0][m]=abs(droplet_centroid[N][0]-droplet_centroid_m[m][0])
    df4_x[0][m]=droplet_centroid_m[m][1]    

plt.figure(figsize=(12,10))
plt.plot(df4_x*1023/(b-a),df4*1023/(b-a),'ob')
plt.xlabel('X Position (um)')
plt.ylabel('Deflection (um)')
plt.show()


#process fifth set of image files

clear_all
 
a=162 ##upper edge position
b=395 ##lower edge position
c=19  ##starting point position
d=669 ## ending edge postion
a1=163
b1=396
c1=19
d1=669
    
#read image stacks and cut
img=T.imread( u'/Users/username/Desktop/2.16-1.44-45/f-50.tif')
img_m=T.imread(u'/Users/username/Desktop/2.16-1.44-45/fm-50.tif')

img_cut=np.ones((img.shape[0],b-a,650)) 
img_m_cut=np.ones((img.shape[0],b1-a1,650))
for i1 in range(img.shape[0]):
    img_cut[i1]=cut_ROI(img[i1],a,b,c,d)
for i2 in range(img_m.shape[0]):
    img_m_cut[i2]=cut_ROI(img_m[i2],a1,b1,c1,d1)    
    
#get the size and center of the droplets
droplet_size=[]
droplet_centroid=[]
droplet_centroid_m=[]

for j in range(img_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_cut[j])
    mask=threshold<img_cut[j]
    label_img=skimage.measure.label(mask)
    props1=skimage.measure.regionprops(label_img)
    for k in range(len(props1)):
        if props1[k].equivalent_diameter*1023/(b-a)>70 and  props1[k].equivalent_diameter*1023/(b-a)<180:
            droplet_size.append(props1[k].equivalent_diameter*1023/234)
            droplet_centroid.append(props1[k].centroid)
        else:
            continue

for h in range(img_m_cut.shape[0]):
    threshold=1*filters.threshold_yen(img_m_cut[h])
    mask=threshold<img_m_cut[h]
    label_img_m=skimage.measure.label(mask)
    props2=skimage.measure.regionprops(label_img_m)
    for l in range(len(props2)):
        if props2[l].equivalent_diameter*1023/(b1-a1)>70 and  props2[l].equivalent_diameter*1023/(b1-a1)<180:
            droplet_centroid_m.append(props2[l].centroid)
        else:
            continue
droplet_size_arr=np.array(droplet_size)
plt.figure(figsize=(12,10))
plt.hist(droplet_size_arr,bins=10)
plt.xlabel('Droplet_Equivalent_Diameter (\u03bcm)')
plt.ylabel('Count')
plt.show()


droplet_centroid_arr=np.array(droplet_centroid)
droplet_centroid_m_arr=np.array(droplet_centroid_m)
stack=np.max(img_cut,axis=0)
stack_m=np.max(img_m_cut,axis=0)
stack_temp=np.zeros((2,img_cut.shape[1],img_cut.shape[2]))
stack_temp[0]=stack
stack_temp[1]=stack_m
stack_all=np.max(stack_temp,axis=0)
stack_all_zoom=scipy.ndimage.zoom(stack_all,1023/(b-a))


plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
plt.gca().invert_yaxis()
plt.imshow(stack_all_zoom,'gray')
plt.plot(droplet_centroid_arr[:,1]*1023/(b-a),droplet_centroid_arr[:,0]*1023/(b-a),'ob',markersize=2)
plt.plot(droplet_centroid_m_arr[:,1]*1023/(b1-a1),droplet_centroid_m_arr[:,0]*1023/(b1-a1),'or',markersize=2)
plt.axis([0,650*1023/(b-a),1023,0]) 

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)
plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Y Position (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
plt.savefig(u'/Users/username/Desktop/5.tiff',bbox_inches='tight')

#fifth set of droplets
size5=np.mean(droplet_size_arr)
dis=np.zeros((1,len(droplet_centroid)))
df5=np.zeros((1,len(droplet_centroid_m)))
df5_x=np.zeros((1,len(droplet_centroid_m)))
for m in range(len(droplet_centroid_m)):
    for n in range(len(droplet_centroid)):
        dis[0][n]=abs(droplet_centroid[n][1]-droplet_centroid_m[m][1])
    N=dis.argmin()
    df5[0][m]=abs(droplet_centroid[N][0]-droplet_centroid_m[m][0])
    df5_x[0][m]=droplet_centroid_m[m][1]    

plt.figure(figsize=(12,10))
plt.plot(df5_x*1023/(b-a),df5*1023/(b-a),'ob')
plt.xlabel('X Position (um)')
plt.ylabel('Deflection (um)')
plt.show()


plt.figure(figsize=(12,10))
rc('axes', linewidth=2)
df1_plot=plt.scatter(df1_x*1023/234,df1*1023/234,s=100,marker='o',c='b')
df2_plot=plt.scatter(df2_x*1023/234,df2*1023/234,s=100,marker='^',c='g')
df3_plot=plt.scatter(df3_x*1023/234,df3*1023/234,s=100,marker='v',c='r')
df4_plot=plt.scatter(df4_x*1023/234,df4*1023/234,s=100,marker='<',c='y')
df5_plot=plt.scatter(df5_x*1023/234,df5*1023/234,s=100,marker='>',c='k')
plt.legend((df1_plot,df2_plot,df3_plot,df4_plot,df5_plot),(u'49 \u03bcm',u'78 \u03bcm',u'95 \u03bcm',u'106 \u03bcm',u'123 \u03bcm'),loc='upper left',scatterpoints=2,markerscale=2,handletextpad=1,fontsize=30)
plt.axis([0,650*1023/234,0,1023])

fontProperties={'family':'sans-serif','sans-serif':['Arial'],'weight' : 'black', 'size' :50}
rc('font',**fontProperties) 

ax = gca()
ax.set_xticklabels(ax.get_xticks(), fontProperties)
ax.set_yticklabels(ax.get_yticks(), fontProperties)
plt.xlabel(u'X Position (\u03bcm)',fontsize=50,fontweight='bold')
plt.ylabel(u'Deflection (\u03bcm)',fontsize=50,fontweight='bold')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0 + 0.03,  pos1.width, pos1.height] 
ax.set_position(pos2) # set a new position
plt.savefig(u'/Users/username/Desktop/deflection.tiff',bbox_inches='tight')