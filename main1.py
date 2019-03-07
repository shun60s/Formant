#coding:utf-8

#
#   LPC対数スペクトを３Ｄプロットして
#   ホルマントの候補位置を赤●でマークする
#

import argparse
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from get_fp import *

#Check version
# Python 3.6.4, 64bit on Win32 (Windows 10)
# numpy (1.14.0)
# matplotlib (2.1.1)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='3D draw of spectrum and formants')
	parser.add_argument("wav_file", help="python main1.py wav-file-name(16khz,mono,16bit)")
	args = parser.parse_args()
	
	# instance class
	fp0=Class_get_fp()
	
	#
	INPUT_WAVE_FILE0=args.wav_file  # ="wav/aiueo1.wav"  # 「あいうえお」声の例
	
	# calculate lpc log-spectrum and  formant, pitch
	spec_out, fout_index,pout= fp0.get_fp(INPUT_WAVE_FILE0)
	
	# 3D plot
	x = np.arange(0,  spec_out.shape[1] , 1) * fp0.df0
	y = np.arange(0,  spec_out.shape[0] , 1) * fp0.dt0
	X, Y = np.meshgrid(x, y)
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.set_xlabel('frequency[Hz]')
	ax.set_ylabel('time[sec]')
	ax.set_zlabel('level')
	ax.plot_wireframe(X,Y,spec_out)
	
	# formant marked as a color
	x0=np.zeros(spec_out.shape[0] * fp0.max_num_formants)
	y0=np.zeros(spec_out.shape[0] * fp0.max_num_formants)
	z0=np.zeros(spec_out.shape[0] * fp0.max_num_formants)
	colorlist = ["r", "g", "b", "c", "m", "y", "k", "w"]
	
	for j in range( fp0.max_num_formants):
		for i in range( spec_out.shape[0] ):
			if fout_index[i][j] >= 0: # if detected
				x0[i + j * spec_out.shape[0]]= fout_index[i][j] * fp0.df0
				y0[i + j * spec_out.shape[0]]= fp0.dt0 * i
				z0[i + j * spec_out.shape[0]]= spec_out[i][ int(fout_index[i][j]) ]
	
	ax.plot( x0, y0 , z0, "o", color=colorlist[0], ms=2, mew=0.5)
	
	# plot show
	plt.show()
	

#This file uses TAB
