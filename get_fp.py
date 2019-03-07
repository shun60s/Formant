#coding: utf-8


import numpy as np
import scipy.signal
import wave
from LPC import *

#Check version
# Python 3.6.4, 64bit on Win32 (Windows 10)
# numpy (1.14.0)
# scipy (1.0.0)

class Class_get_fp(object):
	def __init__(self,NFRAME=640, NSHIFT=320, lpcOrder=32, FreqPoints=1024, max_num_formants=5):
		
		self.NFRAME=NFRAME     # 640 sr=16Khz 40mS  # 400 sr=16Khz 25mS 
		self.NSHIFT=NSHIFT     # 320 sr=16Khz 20mS  # 160 sr=16khz 10mS
		self.lpcOrder=lpcOrder
		self.FreqPoints = FreqPoints  # need many number for precise analysis especially pitch detect
		self.window = np.hamming(self.NFRAME)  # Windows is Hamming
		self.preemph=0.97  # pre-emphasis 
		self.max_num_formants =max_num_formants   # maximum number of formant candidate to detect


	def get_fp(self,file_name ):
		#   入力：wave ファイル　mono 16bit
		#
		#   出力：LPC対数スペクト周波数の行列
		#         ホルマント周波数の候補のインデックス
		#         ピッチ周波数の候補
		#         
		
    	# read wave file
		waveFile= wave.open( file_name, 'r')
		
		nchannles= waveFile.getnchannels()
		samplewidth = waveFile.getsampwidth()
		sampling_rate = waveFile.getframerate()
		nframes = waveFile.getnframes()
		self.df0 = (sampling_rate /2.) / self.FreqPoints
		self.dt0 = 1.0 / sampling_rate
		
		# check input wave file condition
		assert nchannles == 1, ' channel is not MONO '
		assert samplewidth==2, ' sample width is not 16bit '
		
		buf = waveFile.readframes(-1) # read all at oance
		
		waveFile.close()
		
		# 16bit integer to float32
		data = np.frombuffer(buf, dtype='int16')
		fdata = data.astype(np.float32)
		
		count= int(((nframes - ( self.NFRAME - self.NSHIFT)) / self.NSHIFT))
		
		
		# prepare output
		spec_out= np.zeros([count,self.FreqPoints])
		fout = np.zeros([count,self.max_num_formants])
		fout_index = np.ones([count,self.max_num_formants]) * -1
		pout = np.zeros(count)
		pout_index = np.ones(count) * -1
		
		pos = 0  # position
		countr=0
		
		for loop in range(count):
			
			## copy to avoid original over-change
			frame = fdata[pos:pos + self.NFRAME].copy()
			
			## pre-emphasis
			frame -= np.hstack((frame[0], frame[:-1])) * self.preemph
			## do window
			windowed = self.window * frame
			## get lpc coefficients
			a,e=lpc(windowed, self.lpcOrder)
			## get lpc spectrum
			w, h = scipy.signal.freqz(np.sqrt(e), a, self.FreqPoints)  # from 0 to the Nyquist frequency
			lpcspec = np.abs(h)
			lpcspec[lpcspec < 1.0] = 1.0  # to avoid log(0) error
			loglpcspec = 20 * np.log10(lpcspec)
			spec_out[loop]=loglpcspec # store to output
			## get formant candidate
			f_result, i_result=self.formant_detect(loglpcspec, self.df0)
			
			if len(f_result) > self.max_num_formants:
				fout[loop]=f_result[0:self.max_num_formants]
				fout_index[loop]=i_result[0:self.max_num_formants]
			else:
				fout[loop]=f_result[0:len(f_result)]
				fout_index[loop]=i_result[0:len(f_result)]
				
			
			## calcuate lpc residual error (= input source)
			r_err=residual_error(a, windowed)
			## autocorrelation of lpc residual error (= input source)
			a_r_err=autocorr(r_err)
			a_f_result, a_i_result = self.pitch_detect(a_r_err, self.dt0)
			if len(a_f_result) > 0: # if candidate exist,
				pout[loop]=a_f_result[0]
				pout_index[loop]=a_i_result[0]
			
			## print output of candidates of [formants] pitch,  frequency[Hz]
			if countr == 0:
				print ('candidates of [formants] pitch,  frequency[Hz] ')
			print (fout[loop], pout[loop])
			
			# index count up
			countr +=1
			# next
			pos += self.NSHIFT
		
		return spec_out, fout_index, pout


	def formant_detect(self,input0, df0, f_min=250):
		#   対数スペクトルから
		#   山型（凸）のピークポイントを見つける
		#
		#   入力：対数スペクトル
		#         周波数単位
		#         （オプション）最低の周波数
		#
		#   出力：ピークのインデックス
		#         ピークの周波数
		is_find_first= False
		f_result=[]
		i_result=[]
		for i in range (1,len(input0)-1):
			if f_min is not None and  df0 * i <= f_min :
				continue
			if input0[i] > input0[i-1] and input0[i] > input0[i+1] :
				if not is_find_first :
					f_result.append( df0 * i)
					i_result.append(i)
					is_find_first =True
				else:
					f_result.append( df0 * i)
					i_result.append(i)

		return f_result, i_result

	def pitch_detect(self, input0, dt0, ratio0=0.2, f_min=100, f_max=500):
		# 　自己相関の
		# 　山と谷の両方のピークを求める
		#
		#   入力：lpc予測残差の自己相関
		#         時間単位
		#         （オプション）自己エネルギー0次成分に対する比率（これ以上を対象とする）
		#         （オプション）最低の周波数
		#         （オプション）最大の周波数
		#
		#   出力：最大ピークのインデックス
		#         最大ピークの周波数の値
		#
		#
		is_find_first= False
		f_result=[]
		i_result=[]
		v_result=[]
		for i in range (1,len(input0)-1):
			if np.abs(input0[i]) < np.abs(input0[0] * ratio0):
				continue
			fp= 1.0 / (dt0 * i)
			if f_max is not None  and fp >= f_max :
				continue
			if f_min is not None and  fp <= f_min :
				continue
			if input0[i] > input0[i-1] and input0[i] > input0[i+1] :
				if not is_find_first :
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
					is_find_first =True
				else:
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
			elif input0[i] < input0[i-1] and input0[i] < input0[i+1] :
				if not is_find_first :
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i] )
					is_find_first =True
				else:
					f_result.append( fp)
					i_result.append(i)
					v_result.append( input0[i])
		
		if is_find_first:  # 最大のピークを探す
			a=np.argmax( np.array(v_result))
			f_result2= [ f_result[np.argmax( np.array(v_result))] ]
			i_result2= [ i_result[np.argmax( np.array(v_result))] ]
		else: #　候補なし
			f_result2=[]
			i_result2=[]
			
		return f_result2, i_result2


#This file uses TAB
