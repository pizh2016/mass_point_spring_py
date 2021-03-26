import numpy as np

#~ *********************************************************************
#~ ****************************  CLASS MassPoint ***********************
#~ *********************************************************************
class MassPoint(object):

    def __init__(self):        
	    self.mass = 0.01               #质量
	    self.fixed = False             #是否固定
	    self.position = np.zeros(3)    #位置坐标
	    self.velocity = np.zeros(3)    #速度
	    self.normal = np.array([0,0,1])       #法向
	    self.forces = np.zeros(3)            #受力
	    self.ID = None            
                 


#~ *********************************************************************
#~ ****************************  CLASS Spring **************************
#~ *********************************************************************
class Spring(object):
	def __init__(self):        
	    self.stiffness = 10.0          #弹簧刚度
	    self.initLength = 1.0          #弹簧原长度
	    self.elasticForce = 0.0
	    self.startPoint = None    #起点
	    self.endPoint = None      #终点
	




