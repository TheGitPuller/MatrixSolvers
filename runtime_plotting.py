# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 18:52:01 2020

@author: Alexander Campbell
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 10:12:52 2020

@author: Alexander Campbell
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
  
df_dense_dgemm = pd.read_csv('Matrix/matMatMult.txt', sep =",", header=None)
df_dense_dgemv = pd.read_csv('Matrix/matVecMult.txt', sep =",", header=None)

#df_dense_chol = pd.read_csv('cholesky_dense.txt', sep =",", header=None)
#df_dense_gauss_seid = pd.read_csv('gauss_seidel_dense.txt', sep =",", header=None)
#df_dense_jacobi = pd.read_csv('jacobi_dense.txt', sep =",", header=None)
#df_dense_cg = pd.read_csv('conjugate_gradient_dense.txt', sep =",", header=None)



x_dgemm = df_dense_dgemm[0]
y_mat_dgemm = df_dense_dgemm[1]
y_blas_dgemm = df_dense_dgemm[2]

x_dgemv = df_dense_dgemv[0]
y_mat_dgemv = df_dense_dgemv[1]
y_blas_dgemv = df_dense_dgemv[2]


fig, ax = plt.subplots(1, 2, sharey=False, figsize=(10, 5))
ax[0].plot(x_dgemm[:], y_mat_dgemm[:], 'bo', label='Matrix')
ax[0].plot(x_dgemm[:], y_blas_dgemm[:], 'ro', label='BLAS')
ax[0].set_title('DGEMM', fontsize=12)
ax[0].set_xlabel('array size (rows)', fontsize=12)
ax[0].set_ylabel('computation time (ms)', fontsize=12)
ax[0].legend(loc='best')


ax[1].plot(x_dgemv[:], y_mat_dgemv[:], 'bo', label='Matrix')
ax[1].plot(x_dgemv[:], y_blas_dgemv[:], 'ro', label='BLAS')
ax[1].set_title('DGEMV', fontsize=12)
ax[1].set_xlabel('array size (rows)', fontsize=12)
ax[1].legend(loc='best')


#degree = 3
#poly_coeffs = np.polyfit(xchold[:], ychold[:], degree)
#print('poly_coeffs: ',poly_coeffs)
#p1 = np.poly1d(poly_coeffs)
#x = np.linspace(0., np.max(xchold), 500)
##fig = plt.figure(figsize=(6, 6))
#ax1 = fig.add_subplot(244)
#ax1.set_xlabel('array size (rows)', fontsize=12)
#ax1.set_ylabel('computation time (ms)', fontsize=12)
#ax1.plot(xchold[:-1], ychold[:-1], 'bo', label='Data')
#ax1.plot(x, p1(x), 'r', label=r'$y = {0:.4f}x+{1:.4f}$'.format(poly_coeffs[0], poly_coeffs[1]))
#ax1.set_title('Cholesky dense, poly. deg. %d' %degree, fontsize=12)
##plt.show()
#
#
#
#
#degree = 3
#poly_coeffs = np.polyfit(xgsd[:], ygsd[:], degree)
#print('poly_coeffs: ',poly_coeffs)
#p1 = np.poly1d(poly_coeffs)
#x = np.linspace(0., np.max(xgsd), 500)
##fig = plt.figure(figsize=(6, 6))
#ax1 = fig.add_subplot(245)
#ax1.set_xlabel('array size (rows)', fontsize=12)
#ax1.set_ylabel('computation time (ms)', fontsize=12)
#ax1.plot(xgsd[:-1], ygsd[:-1], 'bo', label='Data')
#ax1.plot(x, p1(x), 'r', label=r'$y = {0:.4f}x+{1:.4f}$'.format(poly_coeffs[0], poly_coeffs[1]))
#ax1.set_title('Dense Gauss-Seidel, poly. deg. %d' %degree, fontsize=12)
##plt.show()
#
#
#
#
#degree = 4
#poly_coeffs = np.polyfit(xjd[:], yjd[:], degree)
#print('poly_coeffs: ',poly_coeffs)
#p1 = np.poly1d(poly_coeffs)
#x = np.linspace(0., np.max(xjd), 500)
##fig = plt.figure(figsize=(6, 6))
#ax1 = fig.add_subplot(246)
#ax1.set_xlabel('array size (rows)', fontsize=12)
#ax1.set_ylabel('computation time (ms)', fontsize=12)
#ax1.plot(xjd[:-1], yjd[:-1], 'bo', label='Data')
#ax1.plot(x, p1(x), 'r', label=r'$y = {0:.4f}x+{1:.4f}$'.format(poly_coeffs[0], poly_coeffs[1]))
#ax1.set_title('Dense Jacobi, poly. deg. %d' %degree, fontsize=12)
##plt.show()
#
#degree = 2
#poly_coeffs = np.polyfit(xcgd[:], ycgd[:], degree)
#print('poly_coeffs: ',poly_coeffs)
#p1 = np.poly1d(poly_coeffs)
#x = np.linspace(0., np.max(xcgd), 500)
##fig = plt.figure(figsize=(6, 6))
#ax1 = fig.add_subplot(247)
#ax1.set_xlabel('array size (rows)', fontsize=12)
#ax1.set_ylabel('computation time (ms)', fontsize=12)
#ax1.plot(xcgd[:-1], ycgd[:-1], 'bo', label='Data')
#ax1.plot(x, p1(x), 'r', label=r'$y = {0:.4f}x+{1:.4f}$'.format(poly_coeffs[0], poly_coeffs[1]))
#ax1.set_title('Dense Conjugate Gradient, poly. deg. %d' %degree, fontsize=12)
#
#
#
#plt.tight_layout()
plt.show()



