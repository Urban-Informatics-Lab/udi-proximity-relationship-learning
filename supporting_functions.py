#Urban Data Integration(UDI) program is used to learn proximity 
#relationships between elements of different classes in a city.    

#Copyright (C) 2017-2018  Karan Gupta,Zheng Yang, Rishee Jain      

#This program is free software: you can redistribute it and/or modify 
#it under the terms of the GNU Affero General Public License as published 
#by the Free Software Foundation, either version 3 of the License, or 
#(at your option) any later version.      

#This program is distributed in the hope that it will be useful, but 
#WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
#Affero General Public License for more details. You should have 
#received a copy of the GNU Affero General Public License along with 
#this program.  If not, see <http://www.gnu.org/licenses/>



from __future__ import absolute_import
import sys
sys.setrecursionlimit(10000)


import math
from math import cos, sin, acos, atan2,asin

from haversine import haversine

from vectors import Point, Vector

def checkIfPointWithinBox(lat_list,long_list,max_lat,min_lat,max_long,min_long): #Checks if the given set of points are within the box defined by the given latitudes and longitudes
	for i in range(0,len(lat_list)):
		if lat_list[i]<=max_lat and lat_list[i]>=min_lat:
			if long_list[i]<=max_long and long_list[i]>=min_long:
				return True 
	return False
def perDistancePointLine_2(a,s1,s2): #return [min_distance,final_point] - Gives minimum distance between point and a line and the point of intersection of perpendicular from the point onto the line
	
	Va=Vector(a[0],a[1],0)
	
	Vs1=Vector(s1[0],s1[1],0)
	Vs2=Vector(s2[0],s2[1],0)
	Vas1=Va.sum(Vs1.multiply(-1))
	Vs1s2=Vs2.sum(Vs1.multiply(-1))
	

	min_distance=(Vas1.cross(Vs1s2)).magnitude()/Vs1s2.magnitude()
	
	perpendicular_point_vector=Vs1s2.multiply(Vas1.dot(Vs1s2)/Vs1s2.magnitude()/Vs1s2.magnitude())
	final_point=perpendicular_point_vector.sum(Vs1)
	return [min_distance,final_point]
	
	
def perDistancePointLine(Va,Vs1,Vs2,Vs1s2):#return [min_distance,final_point] - Gives minimum distance between a point and a line and the point of intersection of a perpendicular from the point to the line
	
	
	
	Vas1=Va.sum(Vs1.multiply(-1))
	
	

	min_distance=(Vas1.cross(Vs1s2)).magnitude()/Vs1s2.magnitude()

	k1=abs(Vas1.dot(Vs1s2))
	k2=Vs1s2.magnitude()
	k=k1/k2/k2
	perpendicular_point_vector=Vs1s2.multiply(k)

	final_point=perpendicular_point_vector.sum(Vs1)

	return [min_distance,final_point]
	
def distanceBetweenSegmentandBuilding(ab,a,af,s1,s2):#return final_result[0][1],final_result[1] - Gives min distance between two adjacent sides of a polygon and a line segment and the point of of shortest distance
	
	Vab=Vector(ab[0],ab[1],0)
	Va=Vector(a[0],a[1],0)
	Vaf=Vector(af[0],af[1],0)
	Vs1=Vector(s1[0],s1[1],0)
	Vs2=Vector(s2[0],s2[1],0)
	Vs1s2=Vs2.sum(Vs1.multiply(-1))
	result=[]
	result.append([perDistancePointLine(Vab,Vs1,Vs2,Vs1s2),Vab])
	result.append([perDistancePointLine(Va,Vs1,Vs2,Vs1s2),Va])
	result.append([perDistancePointLine(Vaf,Vs1,Vs2,Vs1s2),Vaf])
	distances=[result[0][0][0],result[1][0][0],result[2][0][0]]
	final_result=result[distances.index(min(distances))]
	

	return final_result[0][1],final_result[1]



def findNearestRoofPoint(polygon,a,b):#return	min_distance.index(min(min_distance)) - Finds the point on a polygon which is closest to a line segment given by two points
	Va=Vector(a[0],a[1],0)
	Vb=Vector(b[0],b[1],0)
	Vab=Vb.sum(Va.multiply(-1))
	min_distance=[]
	for point in polygon:
		p=Vector(point[0],point[1],0)
		ap=p.sum(Va.multiply(-1))
		min_distance.append((Vab.cross(ap).magnitude())/Vab.magnitude())
		
	#print min(min_distance), "in roof point"
	return	min_distance.index(min(min_distance))
def defineSideOfSegment(building,a, b):# return 1 or 2 - Gives the side of line segment on which an element lies. The element is given by the coordinates of the center of the element. 
	
	m=[]
	try:
		m=(b[1]-a[1])/(b[0]-a[0])
	except ZeroDivisionError:
		
		if (a[0]>0 and a[0]>building[0])or (a[0]<0 and a[0]<building[0]):
			return 1
		else:
			return 2
		
		
	c=a[1]-m*a[0]
	if c>=0:
		if (building[1]-m*building[0]-c)>0:
			return 1
		else:
			return 2
	 #if returning 1 and positive then side 1 ,if returning 1 and negative side 2
	else:
		if (building[1]-m*building[0]-c)<=0:
			return 1
		else:
			return 2
		

def findAdjacentSegmentIndex(segment,index):#return index-1 or index+1 or 1 - Gives the index in the array representing a segment which is next to the given index. 
	if index==0:
		return 1
	elif index==(len(segment)-1):
		return index-1
	else:
		return index+1

def findAdjacentIndex(temp_roof1, vertex_a_index):# return vertex_a_index+1,len(temp_roof1)-1 or (vertex_a_index+1)%len(temp_roof1),vertex_a_index-1 - Finds the indices of the vertices adjacent to the given vertex for a polygon
	if vertex_a_index==0:
		return vertex_a_index+1,len(temp_roof1)-1
	else: 
		return (vertex_a_index+1)%len(temp_roof1),vertex_a_index-1
	
	
def projectOnLine(arg1,arg2,arg3): #return t which is the projection amount between 0 and 1 - find the projection of a point on a line segment
	a=Vector(arg1[0],arg1[1],0)
	b=Vector(arg2[0],arg2[1],0)
	k=Vector(arg3[0],arg3[1],0)
	#a=Vector(1,2,0)
	#b=Vector(3,4,0)
	#k=Vector(2,-8,0)
	axis_vector=b.sum(a.multiply(-1))
	per_vector=Vector(axis_vector.y,-1*axis_vector.x,0)

	lhs=a.sum(k.multiply(-1))

	lhsx=float(lhs.x)
	lhsy=float(lhs.y)
	axisx=float(axis_vector.x)
	axisy=float(axis_vector.y)
	perx=float(per_vector.x)
	pery=float(per_vector.y)
	t=(lhsx-perx/pery*lhsy)/(perx*axisy/pery-axisx)
	return t

def findOverlap(polygon, point1,point2):# returns the amount of overlap a polygon has with a line segment with two points
	overlap=[]
	for i in range(0,len(polygon)):
		overlap.append(projectOnLine(point1,point2,polygon[i]))
	min1=min(overlap)
	max1=max(overlap)
	if min1>=1 or max1<=0:
		return 0
	elif min1<=0 :
		if max1>=1:
			return 1
		elif max1<1:
			return max1
	elif min1>0:
		if max1>=1:
			return 1-min1
		elif max1<1:
			return max1-min1
	else :
		return 0
		

def isInternal(a,b,c):# returns True if a is internal to a and b
	Va=Vector(a[0],a[1],0)
	Vb=Vector(b[0],b[1],0)
	Vc=Vector(c[0],c[1],0)
	Vab=Va.sum(Vb.multiply(-1))
	Vac=Va.sum(Vc.multiply(-1))
	Vbc=Vc.sum(Vb.multiply(-1))
	zero_vector=Point(0,0,0)
	min_distance=100
	final_point=[]
	type=[]
	if Vbc.magnitude()<.0001:
		return False
	if abs(Vab.dot(Vbc)/Vbc.magnitude())<= Vbc.magnitude() and abs(Vac.dot(Vbc)/Vbc.magnitude())<= Vbc.magnitude():
		return True
	else :
		return False
	
	return False
		
def findDistanceBetweenSegments(ab,a,af,bb,b,bf):# returns geo distance between two adjacent line segments for two polygons in m
	
	z_point=Point(0,0,0)
	Paf=Point(af[0],af[1],0)
	Pa=Point(a[0],a[1],0)
	Pab=Point(ab[0],ab[1],0)
	Pbf=Point(bf[0],bf[1],0)
	Pb=Point(b[0],b[1],0)
	Pbb=Point(bb[0],bb[1],0)
	compare=[]
	dist=[]
	
	temp_result=findShortestDistanceSegments(Pab, Pa,Pbb, Pb) #[[min_distance,[final_point.x,final_point.y]],a2]
	compare.append([temp_result[0][1],temp_result[1]])
	dist.append(temp_result[0][0])
	
	temp_result=findShortestDistanceSegments(Pab, Pa,Pb, Pbf)
	compare.append([temp_result[0][1],temp_result[1]])
	dist.append(temp_result[0][0])
	
	temp_result=findShortestDistanceSegments(Pa, Paf,Pb, Pbb)
	compare.append([temp_result[0][1],temp_result[1]])
	dist.append(temp_result[0][0])
	
	temp_result=findShortestDistanceSegments(Pa, Paf,Pb, Pbf)
	compare.append([temp_result[0][1],temp_result[1]])
	dist.append(temp_result[0][0])
	
	min_index=dist.index(min(dist))
	
	final_result=compare[min_index]
	
	return haversine(cartesianToGeo(final_result[0][0],final_result[0][1]),cartesianToGeo(final_result[1].x,final_result[1].y))*1000
	
def ccw(A,B,C):
    return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x)


def intersect(A,B,C,D):# Returns true if line segments AB and CD intersect
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


	
def findShortestDistanceSegments(a1,a2,b1,b2): #returns the min_distance between two segments, perpendicular intersection point and the point from which the distance is calculated
	
	Va=Vector.from_points(a1, a2)
	Vb=Vector.from_points(b1, b2)
	
	if intersect(a1,a2,b1,b2):
		return[[0,[0,0]],a1]
	result=[]
	result.append([distancePointSegment(a1,Vb,b1,b2),a1])
	result.append([distancePointSegment(a2,Vb,b1,b2),a2])
	result.append([distancePointSegment(b1,Va,a1,a2),b1])
	result.append([distancePointSegment(b2,Va,a1,a2),b2])
	distances_list=[result[0][0],result[1][0],result[2][0],result[3][0]]
	final_result=result[distances_list.index(min(distances_list))]
	return final_result


def distancePointSegment(a,Vbc,b,c): #returns min_distance between bc and a and the point of intersection of perpendicular from a on bc
	
	Vab=Vector.from_points(b, a)
	Vac=Vector.from_points(c, a)
	zero_vector=Point(0,0,0)
	min_distance=100
	final_point=[]
	type=[]
	if Vbc.magnitude()<.000001:
		return [9999999,[a.x,a.y]]
	if abs(Vab.dot(Vbc)/Vbc.magnitude())<= Vbc.magnitude() and abs(Vac.dot(Vbc)/Vbc.magnitude())<= Vbc.magnitude():
		min_distance=(Vab.cross(Vbc)).magnitude()/Vbc.magnitude()
		
		perpendicular_point_vector=Vbc.multiply(abs(Vab.dot(Vbc))/Vbc.magnitude()/Vbc.magnitude())
		final_point=perpendicular_point_vector.sum(Vector.from_points(zero_vector,b))
		type="middle"
	else :
		min_distance=min([Vab.magnitude(),Vac.magnitude()])
		if min_distance==Vab.magnitude():
			final_point=Vector.from_points(zero_vector, b)
			type="vertex"
		else:
			final_point=Vector.from_points(zero_vector, c)
			type="vertex"

	return [min_distance,[final_point.x,final_point.y]]

	
def geoToCartesian(temp_lat,temp_long):# returns x and y coordinates in m
	temp_lat=float(temp_lat)*float(math.pi/180)
	temp_long=float(temp_long)*float(math.pi/180)
	return 6371* cos(temp_lat)*cos(temp_long)*1000, 6371*cos(temp_lat)* sin(temp_long)*1000

def cartesianToGeo(tempx,tempy): #returns lat and long from cartesian coordinates
	lon = float(atan2(tempy, tempx))
	temp = float(acos(tempx/6371000/cos(lon)))*float(180/math.pi)
	lon=lon*float(180/math.pi)
	if 90< temp and temp <=180:
		temp =90-temp
	return temp,lon	
def findlatlong(sheet):# returns the column in sheet which has lat and long
	lattemp=[]
	longtemp=[]
	
	for j in range(1, sheet.max_column):
		temp=sheet.cell(row=1,column=j).value.lower()
		if lattemp!=[] and longtemp !=[]:
			break
		if u'latitude' in temp:
			lattemp=j
			continue
		if u'longitude' in temp: 
			longtemp=j
			continue
	return lattemp, longtemp
		
def findSides(polygon):# identifies sides of a polygon and the order of vertices
	temp_polygon=polygon
	temp_x=list(zip(*temp_polygon)[0])
	temp_y=list(zip(*temp_polygon)[1])
	starting_point=temp_x.index(min(temp_x))
	length_x=len(temp_x)
	count=0
	sides=[-1]*length_x
	next_index=(starting_point+1)%length_x
	type="clock"
	if temp_y[next_index] >= temp_y[starting_point]:
		count=0
		while count<length_x:
			sides[(count+starting_point)%length_x]=count+1
			count+=1
	else:
		count=0
		while count<length_x:
			sides[(count+starting_point)%length_x]=length_x-count
			count+=1
		type="anti"
	return [sides,type]

		
		
def extractSide(b_sides,index1,index2):# returns the calculated side of the polygon
		
	k=[]
	if b_sides[1]=="clock" :
		if min(index1,index2)==0 and max(index1,index2)==len(b_sides[0])-1:
			k= b_sides[0][max(index1,index2)]
		else:
			k =b_sides[0][min(index1,index2)]
	else:
		if min(index1,index2)==0 and max(index1,index2)==len(b_sides[0])-1:
			k= b_sides[0][min(index1,index2)]
		else:
			k= b_sides[0][max(index1,index2)]

	return k
							
		
def findattribute(sheet,string):# returns the column which is the given attribute in the table

	for j in range(1, sheet.max_column):
		if string.strip().lower() in sheet.cell(row=1,column=j).value.lower():
			return j
	