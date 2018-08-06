python

import math
from collections import defaultdict
import time
import pandas as pd

a_name = '5t04'
i_name = '3zev'
tm = 5
ALIGNMENT = True

data = pd.read_csv('GPCR-TM-table-identity-resis-pair.csv', index_col=None)

def get_row_index(data, a_name, i_name):
	for idx in data.index:
		cur_a_name = data['pdb_active'][idx]
		cur_i_name = data['pdb_inactive'][idx]
		if cur_i_name[:4] == i_name[:4] and cur_a_name[:4] == a_name[:4]:
			return idx	

def clear_rms_array(Y):
	T = Y
	for idx, y in enumerate(T):

		if y == -1 or y == 0:
			if (idx - 1) >= 0 and (idx + 1) < len (T):
				if T[idx + 1] != -1 and T[idx + 1] != 0:
					T[idx] = (T[idx + 1] + T[idx - 1]) / 2.
				else:
					T[idx] = T[idx - 1]
			elif (idx - 1) > 0:
				T[idx] = T[idx -1]
			elif (idx + 1) < len(T):
				if T[idx + 1] != -1 and T[idx + 1] != 0:
					T[idx] = T[idx + 1]
				else:
					T[idx] = [i for i in Y if i > 0][0]
	return T
                
def smooth_data(Y):
    sum_rms = []
    for idx, y in enumerate(Y):
        if (idx + 4) < len(Y):
            s = sum(Y[idx:idx + 5]) / 5.
            sum_rms.append([idx, s])
    return sum_rms

def get_min_rms(smoothed_data):
    y = [x[1] for x in smoothed_data]
    min_rms = min(y)
    for element in smoothed_data:
        if element[1] == min_rms:
            return element
        
def get_min_resi(min_rms, X, Y):
    idx_base = min_rms[0]
    y_min = None
    if (idx_base + 4) < len(Y):
        y_min = min(Y[idx_base:idx_base + 5])
    else:
        y_min = min(Y[idx_base:])

    for x, y in zip(X, Y):
        if y == y_min:
            return x, y
        
def get_resn_by_resi(tg_resi, target_resn):
    
    cmd.iterate('model gpcr_i and resi ' + str(tg_resi), 'target_resn.add(resn)')
    return list(target_resn)[0]

def get_angle(X, tg_resi):
    angle_b = cmd.get_angle('model tm6_i and resi ' + str(X[0]) + ' and name CA', 'model tm6_i and resi ' + str(tg_resi) + ' and name CA', 'model tm6_a and resi ' + str(X[0]) + ' and name CA')
    angle_a = cmd.get_angle('model tm6_i and resi ' + str(X[-1]) + ' and name CA', 'model tm6_i and resi ' + str(tg_resi) + ' and name CA', 'model tm6_a and resi ' + str(X[-1]) + ' and name CA')
    
    return angle_b, angle_a

cmd.reinitialize()

cmd.fetch(i_name, async=1)
cmd.fetch(a_name, async=1)

time.sleep(3)

idx = get_row_index(data, a_name, i_name)
a_chain = data['chain_a'][idx]
i_chain = data['chain_i'][idx]
tm_i_begin = int(data['TM' + str(tm) + '_inactive'][idx].split('-')[0])
tm_i_end = int(data['TM' + str(tm) + '_inactive'][idx].split('-')[1])
tm_a_begin = int(data['TM' + str(tm) + '_active'][idx].split('-')[0])
tm_a_end = int(data['TM' + str(tm) + '_active'][idx].split('-')[1])
i_begin = int(data['TM1_inactive'][idx].split('-')[0])
i_end = int(data['TM7_inactive'][idx].split('-')[1])
a_begin = int(data['TM1_active'][idx].split('-')[0])
a_end = int(data['TM7_active'][idx].split('-')[1])

cmd.create('gpcr_a', 'model ' + a_name + ' and chain ' + a_chain + ' and resi ' + str(a_begin) + '-' + str(a_end) + ' and alt A+\"\"')
cmd.create('gpcr_i', 'model ' + i_name + ' and chain ' + i_chain + ' and resi ' + str(i_begin) + '-' + str(i_end) + ' and alt A+\"\"')

cmd.delete(i_name)
cmd.delete(a_name)

if ALIGNMENT:
	cmd.align('gpcr_i', 'gpcr_a', transform=1)

cmd.create('tm' + str(tm) + '_a', 'model gpcr_a and chain ' + a_chain + ' and resi ' + str(tm_a_begin) + '-' + str(tm_a_end) + ' and alt A+\"\"')
cmd.create('tm' + str(tm) + '_i', 'model gpcr_i and chain ' + i_chain + ' and resi ' + str(tm_i_begin) + '-' + str(tm_i_end) + ' and alt A+\"\"')

cmd.delete(i_name)
cmd.delete(a_name)

if not ALIGNMENT:
	cmd.align('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a')

resns_a = set()
resns_i = set()
rms_cur_data = list()

for j in range(tm_i_begin, tm_i_end + 1):

	tm_a_j = None
	if tm_a_begin > 1000 and tm_i_begin < 1000:
		tm_a_j = j + 1000
	elif tm_a_begin < 1000 and tm_i_begin > 1000:
		tm_a_j = j - 1000
	else:
		tm_a_j = j
	cmd.iterate('model ' + 'tm' + str(tm) + '_a and resi ' + str(tm_a_j), 'resns_a.add(resn)')
	cmd.iterate('model ' + 'tm' + str(tm) + '_i and resi ' + str(j), 'resns_i.add(resn)')
	if not resns_a or not resns_i:
		rms_cur_data.append([j, 0])
		continue

	resn_a = list(resns_a)[0]
	resn_i = list(resns_i)[0]
	resns_a.clear()
	resns_i.clear()
	if resn_a != resn_i:
		rms_cur_data.append([j, 0])
		continue

	rms_cur = cmd.rms_cur('model ' + 'tm' + str(tm) + '_i and resi ' + str(j) + ' and alt A+\"\"', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(tm_a_j) + ' and alt A+\"\"', matchmaker = 1)

	if rms_cur == -1:
			rms_cur = cmd.get_distance('model ' + 'tm' + str(tm) + '_i and resi ' + str(j) + ' and name CA and alt A+\"\"', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(tm_a_j) + ' and name CA and alt A+\"\"')
	rms_cur_data.append([j, rms_cur])

X = [x[0] for x in rms_cur_data]
Y = [x[1] for x in rms_cur_data]

data = pd.DataFrame(columns=['resi', 'rms'])

for x, y in zip(X, Y):
	s = pd.Series({'resi': x, 'rms': y})	
	data = data.append(s, ignore_index=True)

data.to_csv('aligned_residues.csv', sep = ',')
print 'done'
python end
