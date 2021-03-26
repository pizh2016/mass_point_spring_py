from massspring import *
import numpy as np
import numba as nb
import time

def normalize(arr,axis=None):
    """
    Normalize an array by L2 norm
    """
    return arr/np.linalg.norm(arr,axis=axis)

#创建质点和弹簧数组
point_array = np.zeros(0,dtype=MassPoint)
spring_array = np.zeros(0,dtype=Spring)
pointCount = 0
springCount = 0

initLength = 1.0
dampFactor = 2.0
Kstructure, Kshear, Kbending = 10,10,50
gravity = np.array([0., 0., -0.98])

def initSquare(gridx,gridy):
    global pointCount, springCount, point_array, spring_array
    #计算质点的数量
    pointCount = gridy * gridx
    #计算弹簧的数量,S1——结构，S2——剪切，S3——弯曲，
    S1 = (gridx - 1) * gridy + (gridy - 1) * gridx 
    S2 = (gridx - 1) * (gridy - 1) * 2
    S3 = (gridx - 2) * gridy + (gridy - 2) * gridx
    springCount = S1 + S2 + S3
    point_array = np.zeros(pointCount,dtype=MassPoint)
    spring_array = np.zeros(springCount,dtype=Spring)
    
    #初始化质点位置
    for i in range(gridy):
        for j in range(gridx):
            point_array[i*gridx + j] = MassPoint()
            point_array[i*gridx + j].position = np.array([i, j, 2.0])
            point_array[i*gridx + j].ID = i*gridx + j
    
    point_array[0].fixed = True
    point_array[gridx-1].fixed = True
    point_array[gridx*(gridy-1)].fixed = True
    point_array[pointCount-1].fixed = True
    sp = 0
    #创建结构弹簧  左右 上下
    for i in range(gridy):
        for j in range(gridx-1):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = i*gridx + j + 1
            spring_array[sp].stiffness = Kstructure          
            spring_array[sp].initLength = initLength
            sp += 1
    for i in range(gridy-1):
        for j in range(gridx):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = (i+1)*gridx + j
            spring_array[sp].stiffness = Kstructure          
            spring_array[sp].initLength = initLength
            sp += 1 

    #创建剪切弹簧 斜下 斜上
    for i in range(gridy-1):
        for j in range(gridx-1):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = (i+1)*gridx + j + 1
            spring_array[sp].stiffness = Kshear          
            spring_array[sp].initLength = initLength*(2)**0.5
            sp += 1  
    for i in range(gridy-1):
        for j in range(1,gridx):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = (i+1)*gridx + j - 1
            spring_array[sp].stiffness = Kshear          
            spring_array[sp].initLength = initLength*(2)**0.5
            sp += 1  

    #创建弯曲弹簧 左右 上下
    for i in range(gridy):
        for j in range(gridx-2):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = i*gridx + j + 2
            spring_array[sp].stiffness = Kbending          
            spring_array[sp].initLength = initLength*2
            sp += 1  
    for i in range(gridy-2):
        for j in range(gridx):
            spring_array[sp] = Spring()
            spring_array[sp].startPoint = i*gridx + j
            spring_array[sp].endPoint = (i+2)*gridx + j
            spring_array[sp].stiffness = Kbending          
            spring_array[sp].initLength = initLength*2
            sp += 1  
    return

def calculatesystem(T,Nt):
    """
    Solve MassPoint-Spring System

    Parameters
    ----------
    T  : Total calculating time. Unit: s 
    Nt : Segment number of Total time

    """
    global point_array, spring_array
    dt = T / Nt
    for k in range(Nt+1):
        tt = k*dt
        #计算弹力
        for i in range(springCount):
            springLength = np.linalg.norm(point_array[spring_array[i].startPoint].position - point_array[spring_array[i].endPoint].position)
            elongation = springLength - spring_array[i].initLength
            spring_array[i].elasticForce = spring_array[i].stiffness * elongation
        #计算下一时刻运动状态
        for i in range(pointCount):
            if point_array[i].fixed:
                point_array[i].velocity = np.zeros(3)
                point_array[i].position = point_array[i].position
            else:
                joinforce = gravity*point_array[i].mass
                for j in range(springCount):
                    if spring_array[j].startPoint == i :
                        #弹力方向归一化向量
                        vectorDirection = normalize(point_array[spring_array[j].endPoint].position - point_array[i].position)
                        joinforce += spring_array[j].elasticForce * vectorDirection
                        FVD = point_array[spring_array[j].endPoint].velocity - point_array[i].velocity
                        joinforce += dampFactor*FVD
                    if spring_array[j].endPoint == i:
                        vectorDirection = normalize(point_array[spring_array[j].startPoint].position - point_array[i].position)
                        joinforce += spring_array[j].elasticForce * vectorDirection
                        FVD = point_array[spring_array[j].endPoint].velocity - point_array[spring_array[j].startPoint].velocity
                        joinforce += dampFactor*FVD
                #acceleraton
                point_array[i].forces = joinforce;
                acceleration = joinforce / point_array[i].mass;
                #update the velocity
                point_array[i].velocity = point_array[i].velocity + acceleration*dt;
                #update the position
                point_array[i].position = point_array[i].position + point_array[i].velocity*dt;
                #Bottom(Workbench) collision detection
                if point_array[i].position[2] < 0:
                    point_array[i].position[2] = 0
                    point_array[i].velocity[2] = 0
    return

if __name__=='__main__':
    t1 = time.time()
    initSquare(20,20)
    t2 = time.time()
    print(t2-t1)

    print (point_array)
    t1 = time.time()
    calculatesystem(0.1,100)
    t2 = time.time()
    print(t2-t1)

    position = np.array([p.position for p in point_array])
    #print(position)