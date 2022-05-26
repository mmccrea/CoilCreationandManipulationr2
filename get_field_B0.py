#!/usr/bin/python3

# Sorting points by Jeff and Rosie
# Jeff updating to parse u3 and u2-u3 for double-cos box coil
# Jeff updated to read data from comsol mesh
# Jeff updated to triangulate and plot contours in matplotlib
# Jeff updated to extract contours
# June 16, 2019 Jeff updated to use patch.py classes
# June 17, 2019 now working properly
# June 25, 2019 updates for different graphing
# February 27, 2021 Added rerouting for side pipes
# May 21, 2021 Divide up into multiple coils for deformation studies
# June 1, 2021 Start to add pipes class, add mayavi use for traces
# March 2022: changed pipe class to pipefitting class, and updated some of the graphicals display options.

import numpy as np
import math
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from scipy import interpolate 
import io

from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
import matplotlib.cm as cm

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from patch import *
from PipesFitting import *

from scipy.optimize import curve_fit
from scipy.optimize import minimize

#for getting current date
import datetime

# parse command line arguments
parser = OptionParser()

parser.add_option("-f", "--file", dest="infile",
                  default="data.txt", help="read data from file",
                  metavar="FILE")

parser.add_option("-m", "--mesh", dest="plotmesh",
                  default=False, action="store_true",
                  help="plot where the mesh points are")

parser.add_option("-c", "--contours", dest="contours",
                  default=False, action="store_true",
                  help="show extracted contours")

parser.add_option("-g", "--graph", dest="graph",
                  default=False, action="store_true",
                  help="show graph of quarter magnet and contours")

parser.add_option("-i", "--current", dest="current", default=0.1,
                  help="design current (A) (default = 0.1)")

parser.add_option("-p", "--planes", dest="planes", default=False,
                  action="store_true",
                  help="show field maps in three cut planes")

parser.add_option("-t", "--traces", dest="traces", default=False,
                  action="store_true",
                  help="show 3D view of coil traces")

parser.add_option("-a", "--nou1", dest="nou1", default=False,
                  action="store_true",
                  help="no solution for u1; just use u3 alone")

parser.add_option("-r", "--roi", dest="roi", default=False,
                  action="store_true",
                  help="calculate statistics on ROI")

parser.add_option("-v", "--verbose", dest="verbose", default=False,
                  action="store_true",
                  help="verbose option")

parser.add_option("-s", "--simplify", dest="simplify",
                  default=-1, help="factor for VW simplification",
                  metavar="factor")
# simplify is used to remove points, making it easier to draw and
# faster to calculate the field

parser.add_option("-x", dest="x", default=False, action="store_true",
                  help="use file containing parameters for scalar potential contour locations")

parser.add_option("-1", "--and1", dest="and1", default=False,
                  action="store_true", help="add one parameter to optimization fitting function")

parser.add_option("-w", "--wiggle", dest="wiggle",
                  default=-1, help="sigma to wiggle each point (m)",
                  metavar="sigma")

parser.add_option("-o", "--optimize", dest="optimize",
                  default=False,
                  action="store_true",
                  help="if invoked the the optimization process will be run from the given -i or -x and x.txt will be overwritten.")

(options, args) = parser.parse_args()

# line simplification
if(float(options.simplify)>0): #only import if simplification is to be done
    from simplification.cutil import simplify_coords, simplify_coords_vw

#-------------------------------------------------------------------------------
#open data file and set initial parameters for coil size
with open(options.infile) as stream:
    print("opening file: " , options.infile)
    d=np.loadtxt(stream,comments="%",unpack=True)
    
du1=d[:,~np.isnan(d[2])] # remove NaN's in u1
x_inner,y_inner,u1_inner,u2_inner,u3_inner=du1

du2=d[:,~np.isnan(d[3])] # remove NaN's in u2
x_outer,y_outer,u1_outer,u2_outer,u3_outer=du2

# geometry factors from COMSOL model
a_out = 2.2 # m, outer of double layer
a_in = 2.0 # m, inner of double layer
coil_length = 2.2 #m, length to extend the coil from 2D

# specify levels
current=float(options.current) # amperes; design current = step in scalar potential

#sets to use optimized results
if options.and1:
    if options.x:
        n=len(np.loadtxt("x.txt",ndmin=1))+1
    else:
        n=2
else:
    n=10 # set what you want here, for the order of the polynomial 2n+1
xarray=np.zeros(n)
xarray[0]=current

if options.x:
    xarray_loaded=np.loadtxt("x.txt",ndmin=1)
    print("Loaded file x.txt: ",xarray_loaded)
    current=xarray_loaded[0]
    for i in range(min(n,len(xarray_loaded))):
        xarray[i]=xarray_loaded[i]
    
#def get_levels(current,alpha=0,beta=0):
#    maxu1=np.max(u1_inner) # u1 always goes positive
#    maxu23=np.max(u2_outer-u3_outer)
#    minu23=np.min(u2_outer-u3_outer)
#    maxu3=np.max(u3_outer)
#    minu3=np.min(u3_outer)
#    maxphi=maxu1
#    minphi=minu23

#    nmax=int(maxphi/current+.5)
#    nmin=int(minphi/current-.5)
#    maxlevel=(nmax+.5)*current
#    minlevel=(nmin+.5)*current
#    levels=[]
#    for n in range(nmin,nmax):
#        levels.append((n+.5)*current+(n+.5)**3*alpha/1e6+(n+.5)**5*beta/1e10)
#    return levels

def get_levels(parameters):
    current=parameters[0]
    
    maxu1=np.max(u1_inner) # u1 always goes positive
    maxu23=np.max(u2_outer-u3_outer)
    minu23=np.min(u2_outer-u3_outer)
    maxu3=np.max(u3_outer)
    minu3=np.min(u3_outer)
    maxphi=maxu1
    minphi=minu23

    nmax=int(maxphi/current+.5)
    nmin=int(minphi/current-.5)
    maxlevel=(nmax+.5)*current
    minlevel=(nmin+.5)*current
    levels=[]
    for n in range(nmin,nmax):
        thislevel=(n+.5)*current
        for order in range(1,len(parameters)):
            thislevel=thislevel+(n+.5)**(2*order+1)*parameters[order]/10**(2*(2*order+1))
        levels.append(thislevel)
    return levels

mylevels=get_levels(xarray)

# mask out bad triangles that will be created when automatically
# triangulating outer (concave) region.

polygon_inner=Polygon([(0,0),(a_in/2,0),(a_in/2,a_in/2),(0,a_in/2)])
polygon_outer=Polygon([(a_in/2,0),(a_out/2,0),(a_out/2,a_out/2),(0,a_out/2),(0,a_in/2),(a_in/2,a_in/2)])

tri=Triangulation(x_outer,y_outer)
ntri=tri.triangles.shape[0]

# example of masking a region from https://matplotlib.org/examples/pylab_examples/tripcolor_demo.html

xmid=x_outer[tri.triangles].mean(axis=1)
ymid=y_outer[tri.triangles].mean(axis=1) # finds the center points of each triangle
mask=np.zeros(ntri,dtype=bool)
i=0
for x,y in zip(xmid,ymid):
    if not polygon_outer.contains(Point(x,y)):
        mask[i]=True
    i=i+1
print(mask)
tri.set_mask(mask)

## refiner from https://matplotlib.org/3.1.0/gallery/images_contours_and_fields/tricontour_smooth_delaunay.html
#subdiv=3
#refiner = UniformTriRefiner(tri)
#tri_refi, u3_refi = refiner.refine_field(u3_outer, subdiv=subdiv)
#tri_refi, u23_refi = refiner.refine_field(u2_outer-u3_outer, subdiv=subdiv)
## refine inner coils too... seems to affect levels(?)
tri_inner=Triangulation(x_inner,y_inner)
##refiner_inner = UniformTriRefiner(tri_inner)
##tri_inner_refi, u1_refi = refiner.refine_field(u1_inner, subdiv=subdiv)


# scipy.interpolate attempt
#from scipy.interpolate import griddata
#grid_z2 = griddata(tri, u2_outer-u3_outer, (grid_x, grid_y), method='cubic')

#-------------------------------------------------------------------------------
# function for extracting wire lines from 2D scalara potential difference and winding 3D coils
#-------------------------------------------------------------------------------
def wind_coils(levels):
    '''
    levels is the contour plot levels to plot with and extract wires at.
    This function relies on other defined variables in this script that are not passed explicitly so it should not be moved to an imported module.
    pipe re-routes added in this function before mirroring to make symmetric coil re-routes.
    
    '''
    # define penetrating pipes outside function for plotting
    mypipes=pipelist()
    QuickPipes(mypipes)
    pipe_density = 4#14
    #----------------------------------------------------------------------------
    #making contour plot to extract wires from
    fig,ax1=plt.subplots()

    #ax1.triplot(tri,color='0.7') # if you want to see the triangulation

    #u23_contours=ax1.tricontour(tri_refi,u23_refi,levels=levels)
    u23_contours=ax1.tricontour(tri,u2_outer-u3_outer,levels=levels)
    if (options.plotmesh):
        ax1.plot(x_outer,y_outer,'k.')
    ax1.axis((0,a_out/2,0,a_out/2))
    fig.colorbar(u23_contours,ax=ax1)

    #u3_contours=ax1.tricontour(tri_refi, u3_refi, levels=levels)
    u3_contours=ax1.tricontour(tri,u3_outer,levels=levels)
    if (options.plotmesh):
        ax1.plot(x_outer,y_outer,'k.')

    #u1_contours=ax1.tricontour(x_inner, y_inner, u1_inner, levels=levels)
    u1_contours=ax1.tricontour(tri_inner,u1_inner,levels=levels)
    if (options.plotmesh):
        ax1.plot(x_inner,y_inner,'k.')


    # In python3, there might be empty contours, which I am cutting out
    # using the following commands
    u23_contours.allsegs=[x for x in u23_contours.allsegs if x]
    u3_contours.allsegs=[x for x in u3_contours.allsegs if x]
    u1_contours.allsegs=[x for x in u1_contours.allsegs if x]
    
    if not options.graph:
        plt.close()
    if options.graph:
        fig.suptitle("Wire Winding and Feed Thrus")
        ax1.set_xlabel("x(meter)")
        ax1.set_ylabel("y(meter)")
        ax1.set_aspect('equal', adjustable='datalim')
        plt.savefig('get_field-OutputFiles/plot-WireWinding-2D-%s.png'%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"), dpi=300, bbox_inches='tight')
        #plt.show()

    # extracting all the contours and graphing them
    if (options.contours):
        print("Plotting Contours Starting")
        fig2, (ax3, ax4) = plt.subplots(nrows=2)

        # nseg=len(u23_contours.allsegs)

        for i,cnt in enumerate(u23_contours.allsegs):
            seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
            x=seg[:,0]
            y=seg[:,1]
            ax3.plot(x,y,'.-',color='black',ms=1)
            if(float(options.simplify)>0):
                segsimp=simplify_coords_vw(seg,float(options.simplify))
                xsimp=segsimp[:,0]
                ysimp=segsimp[:,1]
                ax3.plot(xsimp,ysimp,'.-',color='red')
        ax3.axis((0,a_out/2,0,a_out/2))

        for i,cnt in enumerate(u3_contours.allsegs):
            seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
            x=seg[:,0]
            y=seg[:,1]
            ax4.plot(x,y,'.-',color='black',ms=1)
            if(float(options.simplify)>0):
                segsimp=simplify_coords_vw(seg,float(options.simplify))
                xsimp=segsimp[:,0]
                ysimp=segsimp[:,1]
                ax4.plot(xsimp,ysimp,'.-',color='red')
        for i,cnt in enumerate(u1_contours.allsegs):
            seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
            x=seg[:,0]
            y=seg[:,1]
            ax4.plot(x,y,'.-',color='black',ms=1)
            if(float(options.simplify)>0):
                segsimp=simplify_coords_vw(seg,float(options.simplify))
                xsimp=segsimp[:,0]
                ysimp=segsimp[:,1]
                ax4.plot(xsimp,ysimp,'.-',color='red')
        ax4.axis((0,a_out/2,0,a_out/2))
        
        plt.savefig('get_field-OutputFiles/plot-Contours-%s.png'%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"), dpi=300, bbox_inches='tight')
        #plt.show()
        
    # conclusion:  the contours are all there and are ordered in the same way relative to each other
    #----------------------------------------------------------------------------
    # assemble 3D coils
    body_t=coilset()
    body_b=coilset()
    body_tr=coilset()
    body_tl=coilset()
    body_br=coilset()
    body_bl=coilset()

    body_list=[body_t,body_b,body_tr,body_tl,body_br,body_bl]

    #print("There are %d outer contours."%len(u23_contours.allsegs))
    #re-work starting with pipes, then checking intersections, then making loops that step outward.

    for i,cnt in enumerate(u23_contours.allsegs):
        seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
        # these go from inner to outer!
        if(float(options.simplify)<0):
            x=seg[:,0]
            y=seg[:,1]
        else:
            segsimp=simplify_coords_vw(seg,float(options.simplify))
            xsimp=segsimp[:,0]
            ysimp=segsimp[:,1]
            x=xsimp
            y=ysimp
        #setting coil length
        z=[-coil_length/2]*len(y)
        xs=np.flip(x,0)
        ys=np.flip(y,0)
        zs=[coil_length/2]*len(y)
        xnew=x
        ynew=y
        znew=z
        
        xnew=np.concatenate((xnew,xs))
        ynew=np.concatenate((ynew,ys))
        znew=np.concatenate((znew,zs))
    
        mirrored=False
        if (x[0]<0.000001):
            # mirror through yz-plane
            xnew=np.concatenate((xnew,np.delete(-x,0))) # avoid repeating point
            ynew=np.concatenate((ynew,np.delete(y,0)))
            znew=np.concatenate((znew,np.delete(zs,0)))

            xnew=np.concatenate((xnew,-xs))
            ynew=np.concatenate((ynew,ys))
            znew=np.concatenate((znew,z))
            mirrored=True
            points=np.column_stack((xnew,ynew,znew))
            body_t.add_coil(points)
        else:
            # complete loop
            xnew=np.append(xnew,xnew[0])
            ynew=np.append(ynew,ynew[0])
            znew=np.append(znew,znew[0])
            points=np.column_stack((xnew,ynew,znew))
            body_tr.add_coil(points)

        # reflect through xz-plane (vertical reflection)
        ynew=-ynew
        points=np.column_stack((xnew,ynew,znew))
        if mirrored:
            body_b.add_coil(points)
        else:
            body_br.add_coil(points)
            ynew=-ynew # put it back for a sec
            # reflect separate trace through yz-plane
            xnew=-xnew
            # reverse the windings
            xnew=np.flip(xnew,0)
            ynew=np.flip(ynew,0)
            znew=np.flip(znew,0)
            points=np.column_stack((xnew,ynew,znew))
            body_tl.add_coil(points)
            # and through the xz-plane
            ynew=-ynew
            points=np.column_stack((xnew,ynew,znew))
            body_bl.add_coil(points)
        #end of pipe checking loop
        # now for the face plates
        front_face_coil=coilset()
        back_face_coil=coilset()

        #print("There are %d face coils in %d levels."%(len(u1_contours.allsegs),len(u1_contours.levels)))

    if (not options.nou1):
        for i,cnt in enumerate(u1_contours.allsegs):
            seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
            # these go from outer to inner
            x=seg[:,0]
            y=seg[:,1]
            # get the corresponding contour from u3
            seg=u3_contours.allsegs[i][0]
            xs=seg[:,0]
            ys=seg[:,1]
            xnew=np.concatenate((x,np.delete(xs,0))) # xs on the boundary is usually a repeated point
            ynew=np.concatenate((y,np.delete(ys,0)))
            # mirror through yz-plane
            xnew=np.concatenate((xnew,np.delete(np.flip(-xs,0),[0,len(xs)-1]))) # first and last point are repeated
            xnew=np.concatenate((xnew,np.flip(-x,0)))
            ynew=np.concatenate((ynew,np.delete(np.flip(ys,0),[0,len(ys)-1])))
            ynew=np.concatenate((ynew,np.flip(y,0)))
            znew=[-coil_length/2]*len(ynew)#[-a_out/2]*len(ynew)
            points=np.column_stack((xnew,ynew,znew))
            back_face_coil.add_coil(points)
            # mirror through xy-plane
            znew=[coil_length/2]*len(ynew)#[a_out/2]*len(ynew)
            xnew=np.flip(xnew,0)
            ynew=np.flip(ynew,0)
            points=np.column_stack((xnew,ynew,znew))
            front_face_coil.add_coil(points)
            # mirror through xz-plane
            ynew=-ynew
            points=np.column_stack((xnew,ynew,znew))
            front_face_coil.add_coil(points)
            znew=[-coil_length/2]*len(ynew)#[-a_out/2]*len(ynew)
            xnew=np.flip(xnew,0)
            ynew=np.flip(ynew,0)    
            points=np.column_stack((xnew,ynew,znew))
            back_face_coil.add_coil(points)
    else:
        for i,cnt in enumerate(u3_contours.allsegs):
            seg=cnt[0] # if there are multiple contours at same level there will be more than one seg
            if(float(options.simplify)<0):
                xs=seg[:,0]
                ys=seg[:,1]
            else:
                segsimp=simplify_coords_vw(seg,float(options.simplify))
                xsimp=segsimp[:,0]
                ysimp=segsimp[:,1]
                xs=xsimp
                ys=ysimp
            xnew=np.concatenate((xs,np.delete(np.flip(-xs,0),0))) # delete repeated point on axis
            ynew=np.concatenate((ys,np.delete(np.flip(ys,0),0)))
            xnew=np.concatenate((xnew,[xnew[0]])) # really force closing the loop
            ynew=np.concatenate((ynew,[ynew[0]]))
            znew=[-coil_length/2]*len(ynew)#[-a_out/2]*len(ynew)
            points=np.column_stack((xnew,ynew,znew))
            back_face_coil.add_coil(points)
            # mirror through xy-plane
            znew=[coil_length/2]*len(ynew)#[a_out/2]*len(ynew)
            xnew=np.flip(xnew,0)
            ynew=np.flip(ynew,0)
            points=np.column_stack((xnew,ynew,znew))
            front_face_coil.add_coil(points)
            # mirror through xz-plane
            ynew=-ynew
            points=np.column_stack((xnew,ynew,znew))
            front_face_coil.add_coil(points)
            znew=[-coil_length/2]*len(ynew)#[-a_out/2]*len(ynew)
            xnew=np.flip(xnew,0)
            ynew=np.flip(ynew,0)    
            points=np.column_stack((xnew,ynew,znew))
            back_face_coil.add_coil(points)
    
    
    all_coil_list=body_list
    all_coil_list.append(front_face_coil)
    all_coil_list.append(back_face_coil)
	
	#----------------------------------------------------------------------------
    #coil extraction post processing
    
    #rotating z to be vertical in all coil sets
    for section in all_coil_list:
        section.rotate(np.pi/2,0,0)
    
    # re-route all pipe groups
    for section in all_coil_list:
        mypipes.reroute_wires(section,pipe_density=pipe_density)

    # Apply coil deformations, if to be included in design optimization
    # body_tr.move(.1,0,0)
    #body_tr.move(.001,0,0)
	
    #move the front and end faces away from end cap to allow construction.
    front_face_coil.move(0,-.002,0)
    back_face_coil.move(0,.002,0)



    return all_coil_list

all_coil_list=wind_coils(mylevels)




for coil in all_coil_list:
    print("There are %d coils in this coil with %i segments"%(coil.ncoils,coil.GetNumSegments()))

if(options.wiggle>0):
    for coil in all_coil_list:
        coil.wiggle(float(options.wiggle))

#Apply coil deformations post winding
#body_tr.move(.1,0,0)
#body_tr.move(.001,0,0)

if(options.traces):
    plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14) 
    plt.rc('ytick', labelsize=14)
    plt.rc('axes', titlesize=14)     # fontsize of the axes title
    
    fig3 = plt.figure()
    ax5 = fig3.add_subplot(111, projection='3d')
    
    fig3.set_size_inches(6.5,6.5)

    #fig3.suptitle("Wire Paths")
    fig3.tight_layout()
    
    ax5.set_xlabel("z(meter)")
    ax5.set_ylabel("x(meter)")
    ax5.set_zlabel("y(meter)")
    ax5.xaxis.labelpad=29
    ax5.yaxis.labelpad=29
    ax5.zaxis.labelpad=20

    
    ax5.view_init(elev=15., azim=30)
    #colors=['black','grey','darkgrey','silver','lightgrey','whitesmoke','blue','green']
    colors=['black','black','black','black','black','black','blue','green']
    for i,coil in enumerate(all_coil_list):
        coil.draw_coils(ax5,linestyle='-',linewidth=1.2,color=colors[i])
        # coil.draw_coils_mayavi()
    # mlab.show()
    ax5.axis('auto')

    #testing length method
    
    #front_face_coil.draw_coil(3,ax5,'-','red')
    #print(front_face_coil.coils[3].length())

    #body_tr.draw_coil(100,ax5,'-','red')
    #print(body_tr.coils[100].length())

    l=0
    for coil in all_coil_list:
        thisl=coil.length()
        print('Length %f'%thisl)
        l=l+thisl
    print('Total length %f'%l)

    # table from https://bulkwire.com/magnet-wire
    
    ohmkm=21.37
    print('18 AWG is %f Ohm/km'%ohmkm)

    resistance=l*ohmkm/1000
    print('Total resistance %f Ohm'%resistance)

    voltage=current*resistance
    print('Voltage %f V'%voltage)

    power=current*voltage
    print('Power %f W'%power)

    kgkm=7.47
    print('18 AWG is %f kg/km'%kgkm)

    weight=l*kgkm/1000
    print('Weight %f kg'%weight)
    plt.savefig('get_field-OutputFiles/plot-Traces-%s.png'%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"), dpi=300, bbox_inches='tight')
    #plt.show()

def vecb(coil_list,x,y,z):
	# calculate the total magnetic field created by all coils in the list
    bx=0.*x
    by=0.*y
    bz=0.*z
    for coil in coil_list:
        bx_tmp,by_tmp,bz_tmp=coil.b_prime(x,y,z)
        bx=bx+bx_tmp
        by=by+by_tmp
        bz=bz+bz_tmp
    return bx,by,bz

#this sets the current in all wires of the coil set to a common value.
for coil in all_coil_list:
    coil.set_common_current(current)

#sets field plotting range in Tesla
design_field=-4*pi/10*1.e-6
bx,by,bz=vecb(all_coil_list,0.,0.,0.)
central_field=bz
delta_field=1.e-9
min_field=central_field-delta_field
max_field=central_field+delta_field
print("Field Plot Range, central: ", central_field, "+/-" ,delta_field)

if (options.planes):
    # show the inner field

    figtest,(axtest1,axtest2,axtest3)=plt.subplots(nrows=3)
    figtest.set_size_inches(5.5,11)
    figtest.suptitle("Field Magnitude In Coil Interior")
    
    x2d,y2d=np.mgrid[-1.0:1.0:101j,-1.0:1.0:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,x2d,y2d,0.)
    im=axtest1.pcolormesh(x2d,y2d,np.sqrt(bx2d**2+by2d**2+bz2d**2),vmin=abs(min_field),vmax=abs(max_field),shading="nearest")
    #im=axtest1.pcolormesh(x2d,y2d,bx2d,vmin=-3e-6,vmax=3e-6)
    
    figtest.colorbar(im,ax=axtest1,format='%.3e',label="Tesla")
    axtest1.set_xlabel("x (m)")
    axtest1.set_ylabel("y (m)")
    axtest1.set_aspect('equal', adjustable='datalim')

    x2d,z2d=np.mgrid[-1.0:1.0:101j,-1.0:1.0:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,x2d,0.,z2d)
    im=axtest2.pcolormesh(z2d,x2d,np.sqrt(bx2d**2+by2d**2+bz2d**2),vmin=abs(min_field),vmax=abs(max_field),shading="nearest")
    #im=axtest2.pcolormesh(z2d,x2d,by2d,vmin=-3e-6,vmax=3e-6)
    
    figtest.colorbar(im,ax=axtest2,format='%.3e',label="Tesla")
    axtest2.set_xlabel("z (m)")
    axtest2.set_ylabel("x (m)")
    axtest2.set_aspect('equal', adjustable='datalim')

    y2d,z2d=np.mgrid[-1.0:1.0:101j,-1.0:1.0:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,0.,y2d,z2d)
    im=axtest3.pcolormesh(z2d,y2d,np.sqrt(bx2d**2+by2d**2+bz2d**2),vmin=abs(min_field),vmax=abs(max_field),shading="nearest")
    #im=axtest3.pcolormesh(z2d,y2d,by2d,vmin=-3e-6,vmax=3e-6)
    
    figtest.colorbar(im,ax=axtest3,format='%.3e',label="Tesla")
    axtest3.set_xlabel("z (m)")
    axtest3.set_ylabel("y (m)")
    axtest3.set_aspect('equal',adjustable='datalim')

    # show the outer field

    figouter,(axouter1,axouter2,axouter3)=plt.subplots(nrows=3)
    figouter.set_size_inches(5.5,11)
    figouter.suptitle("Field Magnitude outside 2.4 m")

    outer_roi=1.5
    inner_roi=1.2
    figouter.suptitle("Field Magnitude outside %4.2f m"%(inner_roi*2))
    
    x2d,y2d=np.mgrid[-outer_roi:outer_roi:101j,-outer_roi:outer_roi:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,x2d,y2d,0.)
    bmod=np.sqrt(bx2d**2+by2d**2+bz2d**2)
    mask=((abs(x2d)<inner_roi)&(abs(y2d)<inner_roi))
    x2d_masked=np.ma.masked_where(mask,x2d)
    y2d_masked=np.ma.masked_where(mask,y2d)
    im=axouter1.pcolor(x2d_masked,y2d_masked,bmod)

    figtest.colorbar(im,ax=axouter1,format='%.3e',label="Tesla")
    axouter1.set_xlabel("x (m)")
    axouter1.set_ylabel("y (m)")
    axouter1.set_aspect('equal', adjustable='datalim')

    x2d,z2d=np.mgrid[-outer_roi:outer_roi:101j,-outer_roi:outer_roi:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,x2d,0.,z2d)
    bmod=np.sqrt(bx2d**2+by2d**2+bz2d**2)
    mask=((abs(x2d)<inner_roi)&(abs(z2d)<inner_roi))
    x2d_masked=np.ma.masked_where(mask,x2d)
    z2d_masked=np.ma.masked_where(mask,z2d)
    im=axouter2.pcolor(z2d_masked,x2d_masked,bmod)
    
    figtest.colorbar(im,ax=axouter2,format='%.3e',label="Tesla")
    axouter2.set_xlabel("z (m)")
    axouter2.set_ylabel("x (m)")
    axouter2.set_aspect('equal', adjustable='datalim')

    y2d,z2d=np.mgrid[-outer_roi:outer_roi:101j,-outer_roi:outer_roi:101j]
    bx2d,by2d,bz2d=vecb(all_coil_list,0.,y2d,z2d)
    bmod=np.sqrt(bx2d**2+by2d**2+bz2d**2)
    mask=((abs(y2d)<inner_roi)&(abs(z2d)<inner_roi))
    y2d_masked=np.ma.masked_where(mask,y2d)
    z2d_masked=np.ma.masked_where(mask,z2d)
    print("y2d = \n" , y2d)
    print("y2d_masked = \n" , y2d_masked)
    bmod_masked=np.ma.masked_where(mask,bmod)
    im=axouter3.pcolor(z2d_masked,y2d_masked,bmod)
    
    figtest.colorbar(im,ax=axouter3,format='%.3e',label="Tesla")
    axouter3.set_xlabel("z (m)")
    axouter3.set_ylabel("y (m)")
    axouter3.set_aspect('equal', adjustable='datalim')

    #tight_layout() reduces white space and text overlap
    figtest.tight_layout()
    figouter.tight_layout()
    plt.savefig('get_field-OutputFiles/plot-InnerOuterField-%s.png'%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"), dpi=300, bbox_inches='tight')
    #plt.show()

def fitfunc(x,p0,p2,p4,p6):
    return p0+p2*x**2+p4*x**4+p6*x**6

def fitgraph(xdata,ydata,ax):
    popt,pcov=curve_fit(fitfunc,xdata[abs(xdata)<.5],ydata[abs(xdata)<.5])
    ax.plot(xdata,fitfunc(xdata,*popt),'r--',label='$p_0$=%2.1e,$p_2$=%2.1e,$p_4$=%2.1e,$p_6$=%2.1e'%tuple(popt))
    #print(popt)
    return(popt)

def fitnograph(xdata,ydata):
    popt,pcov=curve_fit(fitfunc,xdata[abs(xdata)<.5],ydata[abs(xdata)<.5])
    return(popt)

def dofit(coil_list):

    fig7, (ax71) = plt.subplots(nrows=1)

    points1d=np.mgrid[-1:1:101j]
    
    bx1d,by1d,bz1d=vecb(coil_list,points1d,0.,0.)
    # print("Bz on x = " , bz1d)
    fitgraph(points1d,bz1d,ax71)
    ax71.plot(points1d,bz1d,label='$B_z(x,0,0)$')
    
    bx1d,by1d,bz1d=vecb(coil_list,0.,points1d,0.)
    # print("Bz on x = " , bz1d)
    popt=fitgraph(points1d,bz1d,ax71)
    ax71.plot(points1d,bz1d,label='$B_z(0,y,0)$')
    
    bx1d,by1d,bz1d=vecb(coil_list,0.,0.,points1d)
    # print("Bz on x = " , bz1d)
    fitgraph(points1d,bz1d,ax71)
    ax71.plot(points1d,bz1d,label='$B_z(0,0,z)$')

    ax71.axis((-.5,.5,min_field,max_field))
    ax71.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax71.legend()

    if not options.graph:
        plt.close()
        plt.savefig('plot-dofit-%s.png'%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"), dpi=300, bbox_inches='tight')
    #plt.show()

    return popt[1]

def dofity(coil_list):

    points1d=np.mgrid[-1:1:101j]
    bx1d,by1d,bz1d=vecb(coil_list,0.,points1d,0.)
    popt=fitnograph(points1d,by1d)

    return popt[1]

def dofitxyz(coil_list):
    points1d=np.mgrid[-1:1:101j]
    bx1d,by1d,bz1d=vecb(coil_list,0.,points1d,0.)
    
    popty=fitnograph(points1d,by1d)
    bx1d,by1d,bz1d=vecb(coil_list,points1d,0.,0.)
    
    poptx=fitnograph(points1d,by1d)
    bx1d,by1d,bz1d=vecb(coil_list,0.,0.,points1d)
    
    poptz=fitnograph(points1d,by1d)
    wt2=1/3*.5**3
    wt4=1/5*.5**5
    wt6=1/7*.5**7
    return sqrt((popty[1]**2+poptx[1]**2+poptz[1]**2)*wt2**2
                +(popty[2]**2+poptx[2]**2+poptz[2]**2)*wt4**2
                +(popty[3]**2+poptx[3]**2+poptz[3]**2)*wt6**2)
    #return sqrt((popty[1]**2+poptx[1]**2+poptz[1]**2)*wt2**2)


dofit(all_coil_list)
# Statistics in ROI

if options.roi:
    print
    print('Statistics on the ROI')
    print

    x,y,z=np.mgrid[-.3:.3:61j,-.3:.3:61j,-.3:.3:61j]

    rcell=0.3 # m, cell radius
    hcell=0.1601 # m, cell height
    dcell=0.08 # m, bottom to top distance of cells
    mask=(abs(z)>=dcell/2)&(abs(z)<=dcell/2+hcell)&(x**2+y**2<rcell**2)
    mask_upper=mask&(z>0)
    mask_lower=mask&(z<0)

    # This is used to test the cell dimensions.

    #fig=plt.figure()
    #ax=fig.add_subplot(111,projection='3d')
    #scat=ax.scatter(x[mask_upper],y[mask_upper],z[mask_upper])
    #plt.show()

    bx_roi,by_roi,bz_roi=vecb(all_coil_list,x,y,z)

    print('Both cells')
    by_mask_max=np.amax(by_roi[mask])
    by_mask_min=np.amin(by_roi[mask])
    by_mask_delta=by_mask_max-by_mask_min
    print('The max/min/diff By masks are %e %e %e'%(by_mask_max,by_mask_min,by_mask_delta))
    by_std=np.std(by_roi[mask])
    print('The masked standard deviation of By is %e'%by_std)
    print

    print('Upper cell')
    by_mask_max=np.amax(by_roi[mask_upper])
    by_mask_min=np.amin(by_roi[mask_upper])
    by_mask_delta=by_mask_max-by_mask_min
    print('The max/min/diff By masks are %e %e %e'%(by_mask_max,by_mask_min,by_mask_delta))
    by_std=np.std(by_roi[mask_upper])
    print('The masked standard deviation of By is %e'%by_std)
    print

    print('Lower cell')
    by_mask_max=np.amax(by_roi[mask_lower])
    by_mask_min=np.amin(by_roi[mask_lower])
    by_mask_delta=by_mask_max-by_mask_min
    print('The max/min/diff By masks are %e %e %e'%(by_mask_max,by_mask_min,by_mask_delta))
    by_std=np.std(by_roi[mask_lower])
    print('The masked standard deviation of By is %e'%by_std)
    print

    print('Both cells BT2')
    bt2_roi=bx_roi**2+bz_roi**2
    bt2_ave=np.average(bt2_roi[mask])
    print('The BT2 is %e'%bt2_ave)
    print

def fun(x):
    current=x[0]
    levels=get_levels(x)
    all_coil_list=wind_coils(levels)
    for coil in all_coil_list:
        coil.set_common_current(current)
    fitresult=dofitxyz(all_coil_list)
    err=fitresult**2*1e20
    print('FUN',x,err)
    return err

from scipy.optimize import Bounds

if options.optimize:
    print("Beginning Wire Optmization Process.")
    #sweeping current
    #current_step=0.0001
    #theis=[]
    #errs=[]
    #for i in arange(current-20*current_step,current+20*current_step,current_step):
    #    err=fun([i])
    #    theis.append(i)
    #    errs.append(err)
    #fig8, (ax81) = plt.subplots(nrows=1)
    #ax81.plot(theis,errs)
    #ax81.set_yscale('log')
    #plt.show()

    # fitting current
    boundsarray = []
    for i in xarray:
        boundsarray.append([i-0.01*abs(i),i+0.01*abs(i)])
    print("Initial x array = \n" , xarray)
    print("Fitting bounds array = \n" , boundsarray)
    res=minimize(fun,xarray,bounds=boundsarray,method='Nelder-Mead')
    # res=minimize(fun,[current],method='Nelder-Mead')
    print('res.x',res.x)
    np.savetxt('x.txt',res.x)

coil_colours = ["b", "g", "r", "c", "m", "y","tab:orange","tab:grey"]
if options.graph:
    print("saving all unclosed plots")
    # plot on one figure
    fig1 = plt.figure(figsize=(9.5,9.5))
    for color,coil in zip(coil_colours,all_coil_list):
        ax4 = coil.draw_layout(fig1,arrow=False,poslabel=False,color=color,linewidth=1)
	    
    ax4[3].set_xlim3d(-1.4, 1.4)
    ax4[3].set_ylim3d(-1.4, 1.4)
    ax4[3].set_zlim3d(-1.4, 1.4)
    
    MultiPagePdf("get_field-OutputFiles/AllPlots-%s.pdf"%datetime.datetime.now().strftime("%Y_%m_%d-%H_%M"))
    plt.show()

# output points in csv file for comsol:
# addon = ["a","b","c","d","e","f","g","h"]
# for add,coil in zip(addon,all_coil_list):
  # print("coil"+add)
  # coil.output_csv("get_field-OutputFiles/CoilPathOutput-MountsInitial","coil"+add)
  
#output into openscad format
# for add,coil in zip(addon,all_coil_list):
    # coil.output_scad("output-get_field_B0-%s-snone.scad"%add,thickness=0.0015)
