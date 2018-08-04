python

import math
from collections import defaultdict
import time
import pandas as pd

a_name = u'5g53'
i_name = u'3pwh'
tm = 7

data = pd.read_csv('GPCR-TM-table-identity-resis.csv', index_col=None)

def get_row_index(data, a_name, i_name):
	for idx in data.index:
		cur_a_name = data['pdb_active'][idx]
		cur_i_name = data['pdb_inactive'][idx]
		if cur_i_name[:4] == i_name[:4] and cur_a_name[:4] == a_name[:4]:
			return idx	

def clear_rms_array(Y):
    for idx, y in enumerate(Y):
        print 'here'
        if y == -1:
            if (idx - 1) >= 0 and (idx + 1) < len (Y):
                if Y[idx + 1] != -1:
                    Y[idx] = (Y[idx + 1] + Y[idx - 1]) / 2.
                else:
                    Y[idx] = Y[idx - 1]
            elif (idx - 1) > 0:
                Y[idx] = Y[idx -1]
            elif (idx + 1) < len(Y):
                if Y[idx + 1] != -1:
                    Y[idx] = Y[idx + 1]
                else:
                    Y[idx] = np.mean(np.array(Y))
                
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


cmd.create('tm' + str(tm) + '_a', 'model ' + a_name + ' and chain ' + a_chain + ' and resi ' + str(tm_a_begin) + '-' + str(tm_a_end) + ' and alt A+\"\"')
cmd.create('tm' + str(tm) + '_i', 'model ' + i_name + ' and chain ' + i_chain + ' and resi ' + str(tm_a_begin) + '-' + str(tm_a_end) + ' and alt A+\"\"')

cmd.delete(i_name)
cmd.delete(a_name)

cmd.align('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a')

resns_a = set()
resns_i = set()
rms_cur_data = list()

for j in range(tm_a_begin, tm_a_end + 1):
	cmd.iterate('model tm' + str(tm) + '_a and resi ' + str(j), 'resns_a.add(resn)')
	cmd.iterate('model tm' + str(tm) + '_i and resi ' + str(j), 'resns_i.add(resn)')
	if not resns_a or not resns_i:
		continue

	resn_a = list(resns_a)[0]
	resn_i = list(resns_i)[0]
	resns_a.clear()
	resns_i.clear()
	if resn_a != resn_i:
		continue

	rms_cur = cmd.rms_cur('model tm' + str(tm) + '_a and resi ' + str(j), 'model tm' + str(tm) + '_i and resi ' + str(j), matchmaker=-1)
	rms_cur_data.append([j, rms_cur])

X = [x[0] for x in rms_cur_data]
Y = [x[1] for x in rms_cur_data]
#clear_rms_array(Y)

data = pd.DataFrame(columns=['resi', 'rms'])

for x, y in zip(X, Y):
	s = pd.Series({'resi': x, 'rms': y})	
	data = data.append(s, ignore_index=True)

data.to_csv('aligned_residues.csv', sep = ',')
#smoothed_data = smooth_data(Y)
#min_rms = get_min_rms(smoothed_data)
#x, y = get_min_resi(min_rms, X, Y)
#target_resn = set()
#print x, y, get_resn_by_resi(x, target_resn)
#print X[0], X[-1]
#print get_angle(X, x)
print 'done'
python end
