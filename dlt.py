from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cv2
import numpy as np

edge = 6.3333	#length in mm

x_2d = [[1316, 1236], [928, 1303], [1455, 884], [668, 979], [1067, 699], [989, 503]] 

x_3d = [[edge, edge, 0], [0, edge, -1*edge], [edge*2, edge*2, 0], [0, edge*2, -1*edge*2], [edge, edge*3, -1*edge], [edge*2, edge*3, -1*edge*2]]




def compute_mx(i, x_2d, x_3d):
	axi = [ -1*x_3d[i][0], -1*x_3d[i][1], -1*x_3d[i][2], -1, 0, 0, 0, 0, x_2d[i][0]*x_3d[i][0],  x_2d[i][0]*x_3d[i][1],  x_2d[i][0]*x_3d[i][2], x_2d[i][0]]
	return axi;

def compute_my(i, x_2d, x_3d):
	ayi = [ 0, 0, 0, 0, -1*x_3d[i][0], -1*x_3d[i][1], -1*x_3d[i][2], -1, x_2d[i][1]*x_3d[i][0],  x_2d[i][1]*x_3d[i][1],  x_2d[i][1]*x_3d[i][2], x_2d[i][1]]
	return ayi;
 
def callibrateDlt(x_2d, x_3d):
	ax0 = compute_mx(0, x_2d, x_3d) 
	ax1 = compute_mx(1, x_2d, x_3d) 
	ax2 = compute_mx(2, x_2d, x_3d) 
	ax3 = compute_mx(3, x_2d, x_3d) 
	ax4 = compute_mx(4, x_2d, x_3d) 
	ax5 = compute_mx(5, x_2d, x_3d) 

	ay0 = compute_my(0, x_2d, x_3d) 
	ay1 = compute_my(1, x_2d, x_3d) 
	ay2 = compute_my(2, x_2d, x_3d) 
	ay3 = compute_my(3, x_2d, x_3d) 
	ay4 = compute_my(4, x_2d, x_3d) 
	ay5 = compute_my(5, x_2d, x_3d) 

	M = [ax0, ay0, ax1, ay1, ax2, ay2, ax3, ay3, ax4, ay4, ax5, ay5]	
	u, s, vh = np.linalg.svd(M, full_matrices=False)
	print(len(M))
	return vh;

j = callibrateDlt(x_2d, x_3d)
p = j[11]


P = [[p[0], p[1], p[2], p[3]], [p[4], p[5], p[6], p[7]], [p[8], p[9], p[10], p[11]] ]
Y = [[0], [0], [0], [1]]
result = [[0], [0], [0]]

def multiply(X, Y):
	for i in range(len(X)):
			for j in range(len(Y[0])):
				for k in range(len(Y)):
					result[i][j] += X[i][k] * Y[k][j]
	return result;

v = multiply(P, Y)

# for r in v:
# 		print(r)

def decomposeRQ(P):
	H = [ [p[0], p[1], p[2]], [p[4], p[5], p[6]], [p[8], p[9], p[10]] ]
	h = [ [p[3]], [p[7]], [p[11]] ]
	H_inv = np.linalg.inv(H) 
	Xo = multiply(-1*H_inv, h)
	q, r = np.linalg.qr(H_inv)
	q1 = np.linalg.inv(q)
	r1 = np.transpose(r)
	return Xo, q1, r1;

Xo, R, K = decomposeRQ(P)
print(Xo)
print(R)
print(K)


###### I AM UPDATING THE CODE. WILL SEND IT TO YOU.

