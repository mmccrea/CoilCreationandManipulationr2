'''the goal of the pipes classes is to define a series of cylindrical volumes, then take as an input a patch.py coilset, process over it, identify interference between the coil and pipes and then reroute the coils appropriately

Based on but not backward compatible with pipes.py created by Jeff Martin.

create by M. McCrea

'''
import numpy as np
from scipy.constants import mu_0, pi
import numpy as np
from numpy import sqrt,cos,sin
import math

class pipe:
    '''
      The pipe class allows for checking if a straight line parallel to either the x,y, or z axis, as defined by haxis, offset by a height vpipe in the direction given by vaxis, will intersect with a circle in the plane defined by the two axis.
      
       _     ^ v-axis
    --/-\--  |-> h-axis
      \_/
       Note: This code is limited in that the lines must approach along the given horizontal axis as diagrammed above with the ----- line crossing the circle.
    '''
    def __init__(self,rpipe,hpipe,vpipe,haxis,vaxis,p_min=0,p_max=0):
        self.rpipe=rpipe #radius
        self.vpipe=vpipe #vertical position
        self.hpipe=hpipe #horizontal position
        self.haxis=haxis # 'x', 'y', or 'z' , the real axis for the horizontal direction
        self.vaxis=vaxis # 'x', 'y', or 'z', the real axis for the vertical direction
        #min and max distances along the axis perpendicular to the v and h axes to check for intersections and draw distance.
        self.pmin = p_min
        self.pmax = p_max
        self.inters = []  #list of pipe_intersects() class objects that collects the intersecting lines for the pipe
        
        #determine axis directions
        if self.vaxis == "x":
            self.index_v  = 0
        elif self.vaxis == "y":
            self.index_v  = 1
        elif self.vaxis == "z":
            self.index_v  = 2
        if self.haxis == "x":
            self.index_h  = 0
        elif self.haxis == "y":
            self.index_h  = 1
        elif self.haxis == "z":
            self.index_h  = 2
        #set the perpendicular axis to the axis not used for the other two.
        self.index_p = [value for value in [0,1,2] if (value != self.index_h and value !=self.index_v)][0]
        
    def __str__(self):
        #sets results when an object of the class is used in print()
        return "class pipe, cntr = (%f,%f), rad=%f , haxis= %s, vaxis = %s, height=%f to %f intersects = %i"%(self.hpipe,self.vpipe,self.rpipe,self.haxis,self.vaxis,self.pmin,self.pmax,len(self.inters))

    def check_intersect(self,line):
        '''
            checks the intersection of a line with the pipe.
            
            assumes the central axis of the pipe and the wire plane are perpendicular and the wire is along the horizontal axis of the pipe.
            
            returns true if intersects, false if it does not.
        '''
        #trouble shooting diagnostics
        # print(self)
        # print("index_h = " , self.index_h)
        # print("index_v = " , self.index_v)
        # print("index_p = " , self.index_p)
        # print("vertical:" , self.vpipe-self.rpipe, line[0,self.index_v] , self.vpipe+self.rpipe)
        # print("horizontal:" , self.hpipe-self.rpipe, line[0,self.index_v] , self.hpipe+self.rpipe)
        # print("line = " , line)

        
        if np.array_equal(line[0],line[1]):
            #if point is itself, don't check for intersections.
            return False

        #if the point is the correct plane
        if(abs(line[0,self.index_p] - line[-1,self.index_p])<0.000001):
            #if the line doesn't overlap the height of the cylinder
            if line[0,self.index_p] < self.pmin or line[0,self.index_p] > self.pmax:
                # print("outside axial ends of cylinder")
                return False
            # if the point has the correct height
            # print("height check = " , line[0,self.index_v]<self.vpipe+self.rpipe, " and " , line[0,self.index_v]>self.vpipe-self.rpipe)
            if(line[0,self.index_v]<self.vpipe+self.rpipe and line[0,self.index_v]>self.vpipe-self.rpipe):
                #if it is prefered that the line must either start or end outside the circle:( untested code)
                # if (np.square(line[0,self.index_h]-self.hpipe)+np.square(line[0,self.index_v]-self.vpipe))<np.square(self.rpipe) or (np.square(line[1,self.index_h]-self.hpipe)+np.square(line[1,self.index_v]-self.vpipe))<np.square(self.rpipe):
                    # print("line inside circle")
                    # return False
                #and opposite sides of it along the horizontal:
                # if line[0,self.index_h]>(self.hpipe+self.rpipe) and line[1,self.index_h]>(self.hpipe-self.rpipe):
                    # print("line behind circle")
                    # return False
                # if (line[0,self.index_h]<(self.hpipe+self.rpipe) and line[1,self.index_h]<(self.hpipe-self.rpipe)):
                    # print("line in front of circle")
                    # return False
                # print("   check_intersect:intersection found")
                return True
        # print("general failure to intersect")
        return False


    def reroute_wire(self,pipe_density=14):
        '''
        after pipelist:check_intersects() has been used this function will re-route around all intersections identified for this pipe.
        
        plane_min != plane_max, only re-route wires in that given plane.
        '''
        
        # print("\nstart reroute_wire")
        #no intersects, nothing to do, exit function 
        if len(self.inters)==0:
            # print("   0 wires intersected pipe")
            return 0,0


        pipenew = pipelist()#list of new pipes to check intersections based on re-routed wire paths
        #one or more coil loops or lines in a loop intersect a given pipe:
        loops = []
        lines = []
        delta_v = []
        planes = []
        for inter in self.inters:
            inter.coil.rerouted = inter.coil.rerouted+1 #indicate wire has been re-routed.
            for i1p1,i1p2 in zip(inter.points1,inter.points2):
                loops.append(inter.coil)
                lines.append(np.array([i1p1,i1p2]))
                delta_v.append(i1p1[self.index_v] - self.vpipe)
                planes.append(i1p1[self.index_p])
        
        sort_index = np.argsort(delta_v)

        # print("loops =\n " , loops)
        # print("lines =\n " , lines)
        # print("delta_v =\n " , delta_v)
        # print("planes =\n " , planes)
        # print("np.argsort(delta_v) = " , np.argsort(delta_v))
        sorting = np.argsort(delta_v)
        loops = np.array(loops)[sorting]
        lines = np.array(lines)[sorting]
        delta_v = np.array(delta_v)[sorting]
        planes = np.array(planes)[sorting]
        # print("sorted")
        # print("loops =\n " , loops)
        # print("lines =\n " , lines)
        # print("delta_v =\n " , delta_v)
        # print("planes =\n " , planes)
        
        
        # find all intersection lines in the same plane
        for pl in np.unique(planes):
            loopst = loops[pl==planes]
            linest = lines[pl==planes]
            delta_vt = delta_v[pl==planes]
            # print("pl = " , pl)
            # print("pl==planes = " , pl==planes)
            # print("loopst.size = " , loopst.size)
            # print("linest.size = " , linest.size)
            # print("loopst = " , loopst)
            # print("linest = " , linest)
            
            #for rerouting split into reroutes going over and under
            loopUnder = loopst[delta_vt<=0]
            lineUnder = linest[delta_vt<=0]
            loopOver = loopst[delta_vt>0]
            lineOver = linest[delta_vt>0]
            # print("pl = " , pl)
            # print("loopOver.size = " , loopOver.size)
            # print("lineOver.size = " , lineOver.size)
            # print("loopOver = " , loopOver[::-1])
            # print("lineOver = " , lineOver[::-1])
            wire_space = 0.008
            rad_add_over = 0
            for loop,line in zip(loopOver,lineOver):
                # print("over: loop = " , loop)
                # print("over: line = \n" , line)
                # print("over: rad_add_over = " , rad_add_over)
                self.insert_arc_points(loop,line,rad_add_over,pipe_density=pipe_density)
                rad_add_over=rad_add_over+wire_space
            
            # print("loopUnder.size = " , loopUnder.size)
            # print("lineUnder.size = " , lineUnder.size)
            # print("loopUnder = " , loopUnder[::-1])
            # print("lineUnder = " , lineUnder[::-1])
            rad_add_under = 0
            for loop,line in zip(loopUnder[::-1],lineUnder[::-1]):
                # print("under: loop = " , loop)
                # print("under: line = \n" , line)
                # print("under: rad_add_under = " , rad_add_under)
                self.insert_arc_points(loop,line,rad_add_under,pipe_density=pipe_density)
                rad_add_under=rad_add_under+wire_space
            
            #make new pipe to encompass newly re-routed wires
            
            if rad_add_over != 0 or rad_add_under !=0:
                rad_add = max(rad_add_over,rad_add_under)
                pipenew.add_pipe(self.rpipe+rad_add,
                                 self.vpipe,
                                 self.hpipe,
                                 self.haxis,
                                 self.vaxis,
                                 line[0,self.index_p]-0.0000001,
                                 line[0,self.index_p]+0.0000001)
        #sorting reroutes between end points of lines
        for loop,line in zip(loops,lines):
            self.sort_between_ends(loop,line)
            
        return len(self.inters),pipenew
        
    def gen_arc_points(self,line,radius, pipe_density = 14):
      '''
      taking a pipe object from pipes.py and line defined by an two points in an array of 3d points as used in patch.py, create a numpy array of a circular arc 
      pipe_density = number of points around primary arc
      
      returns points to be inserted into loop.
      '''
      h_around=[] #horizontal motion
      v_around=[] #vertical motion
      p_around=[] #in plane coordinates
      
      #load pipe information
      rpipe = radius
      hpipe = self.hpipe
      vpipe = self.vpipe
      index_h  = self.index_h
      index_v  = self.index_v
      
      #check if re-route needs to be created
      #if(line[0,index_v]<vpipe+rpipe and line[0,index_v]>vpipe-rpipe):
      p_around=[line[0,self.index_p]]*pipe_density #coordinate in plane is constant
      #start calculating 2D point locations for shortest path around pipe:
      if(line[0,index_v]>vpipe): # if above center, go around the top side
        # print("gen_arc_points: going over")
        theta_start=math.atan2(line[0,index_v]-vpipe,sqrt(rpipe**2-(line[0,index_v]-vpipe)**2))
        theta_end=math.pi-theta_start
        theta_around=[theta_start+(theta_end-theta_start)*i/(pipe_density-1) for i in range(0,pipe_density)]
        v_around=vpipe+rpipe*sin(theta_around)
      else: # go around the bottom side if below or equal
        # print("gen_arc_points: going under")
        theta_start=math.atan2(vpipe-line[0,index_v],sqrt(rpipe**2-(line[0,index_v]-vpipe)**2))
        theta_end=math.pi-theta_start
        theta_around=[theta_start+(theta_end-theta_start)*i/(pipe_density-1) for i in range(0,pipe_density)]
        v_around=vpipe-rpipe*sin(theta_around)
      
      h_around=hpipe+rpipe*cos(theta_around)
      
      #recombine into [x,y,z] points appropriately for the horizontal and vertical axes:
      if self.haxis == "x":
        if self.vaxis == "y":
          new_arc = np.array([h_around, v_around, p_around])
        if self.vaxis == "z":
          new_arc = np.array([h_around, p_around, v_around])
      if self.haxis == "y":
        if self.vaxis == "x":
          new_arc = np.array([v_around, h_around, p_around])
        if self.vaxis == "z":
          new_arc = np.array([p_around, h_around, v_around])
      if self.haxis == "z":
        if self.vaxis == "x":
          new_arc = np.array([v_around, p_around, h_around])
        if self.vaxis == "y":
          new_arc = np.array([p_around, v_around, h_around])
      new_arc = np.transpose(new_arc.reshape(3,-1))
      
      #return new arc ensuring start point is closer first point than second point in line to keep direction correct.
      
      #diagnostic checks, uncomment for testing.
      '''
      print("line = " , line)
      print("rpipe = " , rpipe)
      print("hpipe = " , hpipe)
      print("vpipe = " , vpipe)
      print("haxis = " , self.haxis)
      print("vaxis = " , self.vaxis)
      
      print("index_h = " , index_h)
      print("index_v = " , index_v)
      print("line[0,index_v] = " , line[0,index_v])
      print("line[0,index_h] = " , line[0,index_h])
        
      print("h_around = " , h_around)
      print("v_around = " , v_around)
      print("p_around = " , p_around)
      
      print("distance to start " , np.linalg.norm(line[0,:]-new_arc[0,:]))
      print("distance to end " , np.linalg.norm(line[0,:]-new_arc[-1,:]))
      '''
      
      if(np.linalg.norm(line[0,:]-new_arc[0,:]) > np.linalg.norm(line[0,:]-new_arc[-1,:])):
        #print("true, flipped")
        return np.flip(new_arc,axis=0)
      else:
        #print("unflipped")
        return new_arc
    
    def insert_arc_points(self,loop,line,rad_add,pipe_density=14):
        new_arc = self.gen_arc_points(line=line,radius=self.rpipe+rad_add)
        
        #remove points from outside the original line from the arc
        line_start = min(line[0, self.index_h],line[-1, self.index_h])
        line_end = max(line[0, self.index_h],line[-1, self.index_h])
        bool_diag = False
        if np.any(np.logical_and(new_arc[:,self.index_h]>line_start, new_arc[:,self.index_h]<line_end)):
            # bool_diag = True
            # print("line = \n" , line)
            # print("line_start = " , line_start)
            # print("line_end = " , line_end)
            # print("new_arc = " , new_arc)
            # print("new_arc[:,self.index_h]>line_start = " , new_arc[:,self.index_h]>line_start)
            # print("new_arc[:,self.index_h]<line_end = " , new_arc[:,self.index_h]<line_end)
            # print("np.and = " , np.logical_and(new_arc[:,self.index_h]>line_start, new_arc[:,self.index_h]<line_end))
            new_arc = new_arc[np.logical_and(new_arc[:,self.index_h]>line_start, new_arc[:,self.index_h]<line_end)]
            # print("new_arc shortened = " , new_arc , "\n")
        
        ind = np.where((loop.points == line[0]).all(axis=1))[0][0]
        # print("loop = " , loop)
        # print("loop.points =\n" , loop.points)
        # print("line = " , line)
        # print("ind = " , ind)
        # print("full = \n" , full)
        # print("self.gen_arc_points(line=line,radius=self.rpipe+rad_add) =\n" , new_arc)
        #add arc to loop
        loop.points = np.insert(
            loop.points,
            ind+1,
            new_arc,
            axis=0)
        # if bool_diag:print("loop.points after = \n" , loop.points)
    
    def sort_between_ends(self,loop,line):
        '''
        takes as an input a line segment, and sorts all points between the endpoints along one axis to make that line either uniformly increasing or decreasing
        '''
        if np.all(line[0]==loop.points[0]): #check if line starts at begining of loop, and use that as the index to prevent np.where confusion with closed loop.
            start = 0
        else:
            start = np.where((loop.points == line[0]).all(axis=1))[0][0]
        # print("start = " , start)
        if np.all(line[1]==loop.points[-1]):#deal with end of line being end of closed loop
            end = len(loop.points)
        else:
            end = np.where((loop.points == line[1]).all(axis=1))[0][0]
        # print("end = " , end)
        
        resort = loop.points[start+1:end]
        #if there are inserted points between the line segment endpoints, sort them along the horizontal axis:
        if not len(resort) == 0:
            resort = resort[resort[:, self.index_h].argsort()]
            #if ascending sort is incorrect, flip to descending
            if(np.linalg.norm(line[0,:]-resort[0,:]) > np.linalg.norm(line[0,:]-resort[-1,:])):
                # print("need to invert order")
                # print("resorting resort = \n" , np.flip(resort,axis=0))
                resort = np.flip(resort,axis=0)
            loop.points[start+1:end] = resort
        
    def draw_pipe(self,ax,trans=0.2, div_length = 2, div_rad = 14,**kwargs):
      '''
      on the given matplotlib axis ax, draws the cyinder for the pipe in the given transparanecy, color, and divisions of line segments.
      The div arguments give the number of flat panels to divide the circle into.
      '''
      #create cylinder along z-axis
      theta = np.linspace(0, 2*np.pi, div_rad)
      length = np.linspace(self.pmin, self.pmax, div_length)
      theta_grid, h_grid=np.meshgrid(theta, length)
      xs = self.rpipe*np.cos(theta_grid)
      ys = self.rpipe*np.sin(theta_grid)
      zs = h_grid
      # print(" xs = \n" , xs)
      # print(" ys = \n" , ys)
      # print(" zs = \n" , zs)
      
      if self.index_p == 0:
          if self.haxis == "y":
              ax.plot_surface(zs, xs+self.hpipe, ys+self.vpipe,alpha=trans,**kwargs)
          if self.haxis == "z":
              ax.plot_surface(zs, xs+self.vpipe, ys+self.hpipe,alpha=trans,**kwargs)
      if self.index_p == 1:
          cyl= np.array([xs,zs,ys])
          if self.haxis == "x":
              ax.plot_surface(xs+self.hpipe, zs, ys+self.vpipe,alpha=trans,**kwargs)
          if self.haxis == "z":
              ax.plot_surface(xs+self.vpipe, zs, ys+self.hpipe,alpha=trans,**kwargs)
      if self.index_p == 2:
          cyl= np.array([xs,ys,zs])
          if self.haxis == "x":
              ax.plot_surface(xs+self.hpipe, ys+self.vpipe, zs,alpha=trans,**kwargs)
          if self.haxis == "y":
              ax.plot_surface(xs+self.vpipe, ys+self.hpipe, zs,alpha=trans,**kwargs)
      

      return 1

class pipe_intersects:
    '''
    A list of lines segments from one coil that intersect the pipe to which this is attached.
    Used to store the output of check_intersects() and then read by reroute_wire() to make the required wire re-routes.
    '''
    def __init__(self,coil):
        self.coil = coil  #mutable coil loop link
        self.points1 = [] #first point in intersecting line
        self.points2 = [] #second point in intersecting line
    #def __str__(self):
    
    def add_inter(self, point1, point2):
        '''
        add a line segment from the given coil for the associated pipe.s
        '''
        self.points1.append(point1)
        self.points2.append(point2)

class pipelist:
    '''
    class for having a list of pipes and managing their properties with respect to a coil
    '''
    def __init__(self):
        self.pipes=[] #list of pipes to be checked
        self.npipes=len(self.pipes)

    def __str__(self):
        return "class pipelist\n    num_pipes = %i\n    id = %i"%(self.npipes,id(self))

    def add_pipe(self,rpipe,vpipe,hpipe,axish,axisv,pmin,pmax):
        newpipe=pipe(rpipe,vpipe,hpipe,axish,axisv,pmin,pmax)
        self.pipes.append(newpipe)
        self.npipes=len(self.pipes)
        
    def join_lists(self,new_pipes):
        for pipe in new_pipes.pipes:
            self.add_pipe(pipe.rpipe,pipe.vpipe,pipe.hpipe,pipe.haxis,pipe.vaxis,pipe.pmin,pipe.pmax)
        
    def draw_pipes(self,pltax,**kwargs):
        for pipe in self.pipes:
            pipe.draw_pipe(pltax,**kwargs)
            
    def draw_xy(self,ax,div_rad = 14,**plt_kwargs):
        for pipe in self.pipes:
            if pipe.index_p == 2: #only pipes  with axis perpendicular to draw plane
                theta = np.linspace(0, 2*np.pi, div_rad)
                theta_grid = np.meshgrid(theta)[0]
                xs = pipe.rpipe*np.cos(theta_grid)
                ys = pipe.rpipe*np.sin(theta_grid)
                if pipe.haxis == "x":
                    ax.plot(xs+pipe.hpipe,ys+pipe.vpipe,**plt_kwargs)
                if pipe.haxis == "y":
                    ax.plot(xs+pipe.vpipe,ys+pipe.hpipe,**plt_kwargs)
            ax.set_xlabel('x (m)')
            ax.set_ylabel('y (m)')
    def draw_zy(self,ax, div_rad = 14,**plt_kwargs):
        for pipe in self.pipes:
            if pipe.index_p == 0: #only pipes  with axis perpendicular to draw plane
                theta = np.linspace(0, 2*np.pi, div_rad)
                theta_grid = np.meshgrid(theta)[0]
                ys = pipe.rpipe*np.cos(theta_grid)
                xs = pipe.rpipe*np.sin(theta_grid)
                if pipe.haxis == "z":
                    ax.plot(xs+pipe.hpipe,ys+pipe.vpipe,**plt_kwargs)
                if pipe.haxis == "y":
                    ax.plot(xs+pipe.vpipe,ys+pipe.hpipe,**plt_kwargs)
            ax.set_xlabel('z (m)')
            ax.set_ylabel('y (m)')
    
    def draw_xz(self,ax, div_rad = 14,**plt_kwargs):
        for pipe in self.pipes:
            if pipe.index_p == 1: #only pipes  with axis perpendicular to draw plane
                theta = np.linspace(0, 2*np.pi, div_rad)
                theta_grid = np.meshgrid(theta)[0]
                xs = pipe.rpipe*np.cos(theta_grid)
                ys = pipe.rpipe*np.sin(theta_grid)
                if pipe.haxis == "x":
                    ax.plot(xs+pipe.hpipe,ys+pipe.vpipe,**plt_kwargs)
                if pipe.haxis == "z":
                    ax.plot(xs+pipe.vpipe,ys+pipe.hpipe,**plt_kwargs)
            ax.set_xlabel('x (m)')
            ax.set_ylabel('z (m)')
    def draw_layout(self,fig,div_rad=14,**plt_kwargs):
        '''
        drawings front, side, top, and isometric view on the passed figure.
        the built in figure axes list must be either blank in which case 4 subplots are added, or if there are 4 axis, draws on them
        no error handling for bad cases is included.
        '''
        if not fig.axes:
            ax3=[]
            ax3.append(fig.add_subplot(2, 2, 3)) #lower left, xy-plane
            ax3.append(fig.add_subplot(2, 2, 4)) #lower right,yz-plane
            ax3.append(fig.add_subplot(2, 2, 1)) #upper left, xz-plane
            ax3.append(fig.add_subplot(2, 2, 2, projection='3d')) #upper right isometric view
        else:
            ax3 = fig.axes
        self.draw_xy(ax3[0],div_rad=div_rad,**plt_kwargs)
        self.draw_zy(ax3[1],div_rad=div_rad,**plt_kwargs)
        self.draw_xz(ax3[2],div_rad=div_rad,**plt_kwargs)
        self.draw_pipes(ax3[3],**plt_kwargs)
        return ax3

    def check_intersects(self,coils,level=0,plane=0):
        '''
        for each coil go through all pipe and for each coil loop in a coil set (coils) add intersect to each pipe
        '''
        print("Starting pipelist:check_intersects of ", len(self.pipes), " pipes.")
        for index, pipe in enumerate(self.pipes):
            # print("pipelist:check_intersects: pipe " , index," of ", len(self.pipes))
            # print("pipelist:check_intersects: Total Coils = " ,len(coils.coils))
            for index, coil in enumerate(coils.coils):
                # if index%5 ==0:#uncomment to print current coils coil being processed.
                    # print(index,end=" ")
                # print("coil = " , coil)
                # print("coil.rerouted = " , coil.rerouted)
                if coil.rerouted == level:
                    add_flag = False
                    coil_ints = pipe_intersects(coil)
                    for j in range(len(coil.points)-1): #for closed loops
                        line = [coil.points[j],coil.points[j+1]]
                        # print("      line[0] = " , line[0])
                        # print("      line[1] = " , line[1])
                        if abs(line[0][pipe.index_p] - plane)<0.0000001:
                            if np.all(line[0] == line[1]):
                                continue #skip repeated points
                            if pipe.check_intersect(np.array(line)):
                                # print("      intersection found")
                                add_flag = True
                                coil_ints.add_inter(line[0],line[1])
                        # else:
                            # print("Wrong Plane")
                    if add_flag:
                        pipe.inters.append(coil_ints)
                        # print("        check_intersects: pipe.inters added")
    
    def reroute_wires(self,coils,pipe_density=14):
        '''
        performs the re-routing for all pipes that have flagged intersects from check_intersects that must be run first.
        '''
        #clean up coils to prepare for checking for intersections
        coils.make_closed()
        coils.round_all()
        
        #pl sets the plane positions to check for re-routes along the axial plane of each cylindrical pipe.
        for pl in [-1.1, -1, 1, 1.1]:
            #reset variables for new plane.
            newpipes = []
            for coil in coils.coils:
                coil.rerouted=0
            newpipes.append(self)
            go_on = True
            counter = 0
            
            #reroute wires and continue to do so until a rerouting make not changes.
            while go_on:
                # print("\n\n - - - -pl = " , pl , ":  start loop = " , counter , "\n\n")
                newpipes[-1].check_intersects(coils, level = 0,plane = pl)
                
                temp = pipelist()
                changes = 0 #track how many lines are changed
                # print("newpipes[-1] = " , newpipes[-1])
                for index, pipe in enumerate(newpipes[-1].pipes):
                    # print("pipelist:reroute_wires: pipe " , index, " of " , len(self.pipes)-1)
                    # print("   pipe = " , pipe)
                    total_reroutes, temppipes = pipe.reroute_wire(pipe_density=pipe_density)
                    changes = changes + total_reroutes
                    # print("      Total Reroutes:" , total_reroutes)
                    if temppipes != 0:
                        temp.join_lists(temppipes)
                newpipes.append(temp)
                counter = counter+1
                go_on = False #uncomment this line to only do the first resets
                if changes == 0 :go_on = False
        
        #clean up after re-routing
        # coils.make_open()
    
    def insert_arcs(self,coils):
        '''
        generate arcs without regard to actual placement requirements without overlap and  proper ordering.
        superseded by reroute_wires().
        '''
        coil_flag = []
        point_flag = []
        pipe_flag = []
        for pipe in self.pipes:
          print(pipe)
          for coil in coils:
            for line in zip(coil.points, np.roll(coil.points,-1,axis=0)):
              if pipe.check_intersect(np.array(line)):
                  coil_flag.append(coil)
                  pipe_flag.append(pipe)
                  print(line[0])
                  print("where = " , np.where((coil.points == line[0]).all(axis=1)))
                  #point_flag.append(np.where((coil.points == line[0]).all(axis=1)))
                  point_flag.append(line[0])
        for index,coil,pipe in zip(point_flag,coil_flag,pipe_flag):
          ind = np.where((coil.points == index).all(axis=1))
          #print("type arc:" , type(pipe.gen_arc_points(coil.points[ind])))
          print("type points:" , type(coil.points))
          #print("arc" , pipe.gen_arc_points(coil.points[ind]))
          print("points:" , coil.points)
          print("index = " , index)
          print("ind = " , ind)
          print(pipe)
          
          coil.points = np.insert(coil.points,ind[0]+1,pipe.gen_arc_points(coil.points[ind[0]]),axis=0)
          print("new points:\n" , coil.points)

def QuickPipes(mypipes):
    '''
    the standard MSR pipe layouts for easy reference
    use this to add all current feedthrus to a pipelist.
    '''
    #back wall pipes
    # gcc_x=0.62 #m, guide center-to-center in x direction
    # gcc_y=.64 #m, guide center-to-center in y direction
    # gdia=.15 #m, guide diameter
    # mypipes.add_pipe(0.1,0,0,'z','x',-1.21,1.21)
    # mypipes.add_pipe(gdia/2,gcc_x/2,gcc_y/2,'z','x',-1.21,1.21)
    # mypipes.add_pipe(gdia/2,gcc_x/2,-gcc_y/2,'z','x',-1.21,1.21)
    # mypipes.add_pipe(gdia/2,-gcc_x/2,gcc_y/2,'z','x',-1.21,1.21)
    # mypipes.add_pipe(gdia/2,-gcc_x/2,-gcc_y/2,'z','x',-1.21,1.21)

    #all pipe radii below are for the feed throughs themselves. r_add will be added to the radii for clearance on alignment and feed throughs:
    r_add = 0.000
    
    #side wall pipes
    rpipes=[0.0508,           0.0508, #m
            0.03,0.03,0.03,0.03,0.03,
                  0.015,0.015,
            0.03,0.03,0.03,0.03,0.03,
                  0.015,0.015,
            0.03,0.03,0.03,0.03,0.03,
            0.0508,            0.0508
            ] #m
    ypipes=[ 0.840,               0.840,
            -0.4,-0.4,-0.4,-0.4,-0.4,
                  0.14,0.185/2,
             0,   0,   0,   0,   0,
                  -0.14,-0.185/2,
             0.4, 0.4, 0.4, 0.4, 0.4,
            -0.840,             -0.840
            ] #m
    zpipes=[-1.050,            1.050,
            -0.44,-0.22,0,0.22,0.44,
                       0,0,
            -0.44,-0.22,0,0.22,0.44,
                       0,0,
            -0.44,-0.22,0,0.22,0.44,
            -1.050,            1.050
            ] #m
    
    rpipes = np.array(rpipes)+r_add
    for j in range(len(rpipes)):
        mypipes.add_pipe(rpipes[j],zpipes[j],ypipes[j],'y',"z",-1.21,1.21)
    
    #floor supports
    rpipes_floor = [ 
                           0.03/2 ,             0.03/2,
                                 0.03/2 , 0.03/2,
                    0.03/2 ,                           0.03/2,
                           0.03/2 ,             0.03/2,
                    0.03/2 , 0.03/2, 0.03/2 , 0.03/2, 0.03/2 , 0.03/2,
                    0.03/2 , 0.03/2, 0.03/2 , 0.03/2, 0.03/2 , 0.03/2,
                    0.03/2 , 0.03/2, 0.03/2 , 0.03/2, 0.03/2 , 0.03/2,
                    0.03/2 , 0.03/2, 0.03/2 , 0.03/2, 0.03/2 , 0.03/2,
                           0.03/2 ,             0.03/2,
                                 0.03/2 , 0.03/2,
                    0.03/2 ,                           0.03/2,
                           0.03/2 ,             0.03/2]
    
    xpipes_floor = [
                            -1.025,         -1.025,
                                 -0.9025,-0.9025,
                    -0.9,                            -0.9,
                            -0.775,         -0.775,
             -0.5025 ,-0.5025 ,-0.5025 ,-0.5025 ,-0.5025 ,-0.5025 ,
             -0.0625 ,-0.0625 ,-0.0625 ,-0.0625 ,-0.0625 ,-0.0625 ,
              0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 ,
              0.5025 , 0.5025 , 0.5025 , 0.5025 , 0.5025 , 0.5025 ,
                             0.775,          0.775,
                                  0.9025, 0.9025,
                     0.9,                             0.9,
                             1.025,          1.025]
    
    ypipes_floor = [
                                    -0.380 ,           0.380,
                                        -0.0675, 0.0675,
                    -0.9225 ,                                       0.9225,
                                    -0.380 ,         0.380,
                    -0.9225 ,-0.4025 , -0.0675 , 0.0675 ,  0.4025 , 0.9225,
                    -0.9225 ,-0.4025 , -0.0675 , 0.0675 ,  0.4025 , 0.9225,
                    -0.9225 ,-0.4025 , -0.0675 , 0.0675 ,  0.4025 , 0.9225,
                    -0.9225 ,-0.4025 , -0.0675 , 0.0675 ,  0.4025 , 0.9225,
                                    -0.380 ,         0.380,
                    -0.9225 ,                                       0.9225,
                                        -0.0675, 0.0675,
                                    -0.380 ,           0.380]
    
    rpipes_floor = np.array(rpipes_floor)+r_add
    for rad,y,x in zip(rpipes_floor,ypipes_floor,xpipes_floor):
        mypipes.add_pipe(rad,x,y,'y','x',-1.21,1.21)
    
    #central 5 feedthrus
    rpipes_feedthru = [0.0697/2,0.0697/2,0.0697/2,0.0697/2,0.0697/2]
    ypipes_feedthru=[-0.2,0,0,0,0.2]
    xpipes_feedthru=[0,-0.2,0,0.2,0]
    rpipes_feedthru = np.array(rpipes_feedthru)+r_add
    for rad,y,x in zip(rpipes_feedthru,ypipes_feedthru,xpipes_feedthru):
        mypipes.add_pipe(rad,x,y,'y','x',-1.21,1.21)

    #side adjustment feed throughs
    oval_length = 00.06 #m
    oval_radius = 0.02#m
    ypipes_sidesprt = [-0.55,-0.55, 0.55, 0.55,-0.55-oval_length,-0.55-oval_length, 0.55+oval_length, 0.55+oval_length]
    zpipes_sidesprt = [-0.11, 0.11,-0.11, 0.11,-0.11, 0.11,-0.11, 0.11]
    for y,z in zip(ypipes_sidesprt,zpipes_sidesprt):
        mypipes.add_pipe(oval_radius,y,z,'y','z',-1.21,1.21)

    #old floor feed thrus for earlier MSL design.
    # rpipes_floor=[0.0697/2,0.0697/2,0.0697/2,0.0697/2,0.0697/2,
                  # .015,.015,.015,.015,.015,.015,.015,
                  # .015,.015,.015,.015,.015,.015,.015,
                  # .015,.015,.015,.015,.015,.015,.015,
                  # .015,.015,.015,.015,.015,.015,.015] #m
    # ypipes_floor=[-0.2,0,0,0,0.2,
                  # .38,.38,.38,.38,.38,.38,.38,
                  # .795,.795,.795,.795,.795,.795,.795,
                  # -.38,-.38,-.38,-.38,-.38,-.38,-.38,
                  # -.795,-.795,-.795,-.795,-.795,-.795,-.795] #m
    # xpipes_floor=[0,-0.2,0,0.2,0,
                  # -.41-.35-.25,-.41-.35,-.41,0,.41,.41+.35,.41+.35+.25,
                  # -.41-.35-.25,-.41-.35,-.41,0,.41,.41+.35,.41+.35+.25,
                  # -.41-.35-.25,-.41-.35,-.41,0,.41,.41+.35,.41+.35+.25,
                  # -.41-.35-.25,-.41-.35,-.41,0,.41,.41+.35,.41+.35+.25] #m

    # for j in range(len(rpipes_floor)):
        # mypipes.add_pipe(rpipes_floor[j],xpipes_floor[j],ypipes_floor[j],'y','x',-1.21,1.21)