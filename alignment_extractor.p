python

import math
from collections import defaultdict
import time
import pandas as pd
import numpy as np

SMOOTH_STEP = 5

resns_a = set()
resns_i = set()
data = pd.read_csv('GPCR-TM6-table-identity-resis.csv', index_col=None)

def get_angle(X, tg_resi, tm, tm_a_j):
	
	angle_b = cmd.get_angle('model ' + 'tm' + str(tm) + '_i and resi ' + str(X[0]) + ' and name CA', 'model ' + 'tm' + str(tm) + '_i and resi ' + str(tg_resi) + ' and name CA', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(X[0] + tm_a_j) + ' and name CA')
	angle_a = cmd.get_angle('model ' + 'tm' + str(tm) + '_i and resi ' + str(X[-1]) + ' and name CA', 'model ' + 'tm' + str(tm) + '_i and resi ' + str(tg_resi) + ' and name CA', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(X[-1] + tm_a_j) + ' and name CA')
    	return angle_b, angle_a

def clear_rms_array(Y):
	for idx, y in enumerate(Y):

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
                
def smooth_data(Y, step):
	sum_rms = []
	for idx, y in enumerate(Y):
		if (idx + step - 1) < len(Y):
			s = sum(Y[idx:idx + step]) / step
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

def get_str_axis(axis):
	if axis == 'noPRO':
		return 'noPRO'
	elif axis == 'PRO is in the end':
		return 'PRO is in the end'
	elif len(axis) != 3:
		return 'more or less then 3 components'
	else:
		return '(' + str(axis[0]) + ', ' + str(axis[1]) + ', ' + str(axis[2]) + ')'

def get_gpcr_properties(idx):

	i_name = data['pdb_inactive'][idx]
	a_name = data['pdb_active'][idx]
	i_begin = int(data['TM1_inactive'][idx].split('-')[0])
	i_end = int(data['TM7_inactive'][idx].split('-')[1])
	a_begin = int(data['TM1_active'][idx].split('-')[0])
	a_end = int(data['TM7_active'][idx].split('-')[1])
	a_chain = data['chain_a'][idx]
	i_chain = data['chain_i'][idx]

	return i_name, i_begin, i_end, i_chain, a_name, a_begin, a_end, a_chain

def get_tm_borders(idx, tm):

	tm_name = 'TM' + str(tm) + '_'
	tm_i_begin = int(data[tm_name + 'inactive'][idx].split('-')[0])
	tm_i_end = int(data[tm_name + 'inactive'][idx].split('-')[1])
	tm_a_begin = int(data[tm_name + 'active'][idx].split('-')[0])
	tm_a_end = int(data[tm_name + 'active'][idx].split('-')[1])

	return tm_i_begin, tm_i_end, tm_a_begin, tm_a_end

def add_gpcr_properties(data, idx, rmsd_i_a, best_point, best_point_resi, best_axis_b, best_angle_b, best_axis_a, best_angle_a, best_rmsd_t_a, pro_point_resi, pro_axis_b, pro_angle_b, pro_axis_a, pro_angle_a, pro_rmsd_t_a, n_of_pro, tm):

	tm_name = 'TM' + str(tm) + '_'
	data[tm_name + 'RMSD(I_A)'][idx] = rmsd_i_a
	data[tm_name + 'best_rotation_point'][idx] = best_point
	data[tm_name + 'best_point_resi'][idx] = best_point_resi
	data[tm_name + 'best_rotation_axis_b'][idx] = get_str_axis(best_axis_b)
	data[tm_name + 'best_angle_b'][idx] = best_angle_b
	data[tm_name + 'best_rotation_axis_a'][idx] = get_str_axis(best_axis_a)
	data[tm_name + 'best_angle_a'][idx] = best_angle_a
	data[tm_name + 'best_RMSD(T_A)'][idx] = best_rmsd_t_a
	data[tm_name + 'pro_point_resi'][idx] = pro_point_resi
	data[tm_name + 'pro_rotation_axis_b'][idx] = get_str_axis(pro_axis_b)
	data[tm_name + 'pro_angle_b'][idx] = pro_angle_b
	data[tm_name + 'pro_rotation_axis_a'][idx] = get_str_axis(pro_axis_a)
	data[tm_name + 'pro_angle_a'][idx] = pro_angle_a
	data[tm_name + 'pro_RMSD(T_A)'][idx] = pro_rmsd_t_a
	data[tm_name + 'NofPRO'][idx] = n_of_pro

def get_rmsd(obj_1, obj_2, trans):
	rmsd = cmd.rms_cur(obj_1, obj_2, matchmaker = 1)
	return rmsd

def get_resis_resns(model):
	
	cmd.iterate('model ' + model, 'resis[int(resi)] = resn')
	pro_resis = []

	success = False
	for key in resis:
		if resis[key] == 'PRO':
			pro_resis.append(key)
			success = True
	
	if success == False:
		pro_resis = [0]

	return resis, pro_resis

def get_transforms(n_resi, tm_i_begin, tm_i_end, tm):
	cmd.create('tm_i_before_exclusive', 'model tm' + str(tm) + '_i and resi ' + str(tm_i_begin) + '-' + str(n_resi - 1))
	cmd.create('tm_i_after_inclusive', 'model tm' + str(tm) + '_i and resi ' + str(n_resi) + '-' + str(tm_i_end))

	cmd.align('tm_i_before_exclusive', 'tm' + str(tm) + '_a')
	transform_before = cmd.get_object_matrix('tm_i_before_exclusive')
	cmd.align('tm_i_after_inclusive', 'tm' + str(tm) + '_a')
	transform_after = cmd.get_object_matrix('tm_i_after_inclusive')

	cmd.create('tm_i_transformed', 'model tm_i_after_inclusive or model tm_i_before_exclusive')
	cmd.delete('tm_i_before_exclusive')
	cmd.delete('tm_i_after_inclusive')
	
	rmsd = get_rmsd('tm_i_transformed', 'tm' + str(tm) + '_a', 0)
	cmd.delete('tm_i_transformed')

	return transform_before, transform_after, rmsd

def get_axis_angle(matrix):
	cos_angle = (matrix[0] + matrix[5] + matrix[10] - 1) / 2.0
	angle = math.degrees(math.acos(cos_angle))
	if cos_angle == 1:
		axis = ['no rotation']
	else:
		axis = [(matrix[j] - cos_angle) / (1.0 - cos_angle) for j in [0, 5, 10]]

	return axis, angle

def get_best_transforms_old(tm_i_begin, tm_i_end, tm):

	min_rmsd = 100.0
	min_resi = 0.0
	pro_rmsd = 100.0
	pro_resi = 0
	best_transform_b = tuple()
	best_transform_a = tuple()
	pro_transform_b = tuple()
	pro_transform_a = tuple()
	pro_rmsds = []
	pro_transforms_b = []
	pro_transforms_a = []
	pro_nums = []

	resis, pro_resis = get_resis_resns('tm' + str(tm) + '_i')
	pro_status = True

	if len(pro_resis) == 1:
		if pro_resis[0] == 0:
			pro_transform_b = 'noPRO'
			pro_transform_a = 'noPRO'
			pro_status = False
			pro_resi = pro_resis[0]
		elif pro_resis[0] == tm_i_begin or pro_resis[0] == tm_i_end:
			pro_transform_b = 'PRO is in the end'
			pro_transform_a = 'PRO is in the end'
			pro_status = False
			pro_resi = pro_resis[0]
	
	for n_resi in resis:

		if n_resi > tm_i_begin:

			transform_before, transform_after, rmsd = get_transforms(n_resi, tm_i_begin, tm_i_end, tm)
			if rmsd < min_rmsd:
				min_rmsd = rmsd
				min_resi = n_resi
				best_transform_b = transform_before
				best_transform_a = transform_after
			if n_resi in pro_resis and pro_status and len(pro_resis) == 1:
				pro_rmsd = rmsd
				pro_transform_b = transform_before
				pro_transform_a = transform_after
				pro_resi = pro_resis[0]
			elif n_resi in pro_resis and pro_status and len(pro_resis) > 1:
				pro_rmsds.append(rmsd)
				pro_transforms_b.append(transform_before)
				pro_transforms_a.append(transform_after)
				pro_nums.append(n_resi)
				pro_rmsd = min(pro_rmsds)


	if len(pro_resis) > 1:
		for rmsd_it, trans_b_it, trans_a_it, pro_num in zip(pro_rmsds, pro_transforms_b, pro_transforms_a, pro_nums):
			if rmsd_it == pro_rmsd:
				pro_transform_b = trans_b_it
				pro_transform_a = trans_a_it
				pro_resi = pro_num
				

	min_resn = resis[min_resi]

	return min_rmsd, min_resi, min_resn, best_transform_b, best_transform_a, pro_rmsd, pro_resi, pro_transform_b, pro_transform_a, len(pro_resis)

def get_best_transforms(tm_i_begin, tm_i_end, tm_a_begin, tm_a_end, tm):

	min_rmsd = 100.0
	min_resi = 0.0
	pro_rmsd = 100.0
	pro_resi = 0
	best_transform_b = tuple()
	best_transform_a = tuple()
	pro_transform_b = tuple()
	pro_transform_a = tuple()
	pro_rmsds = []
	pro_transforms_b = []
	pro_transforms_a = []
	pro_nums = []

	resis, pro_resis = get_resis_resns('tm' + str(tm) + '_i')
	pro_status = True

	cmd.align('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a')	
	
	if len(pro_resis) == 1:
		if pro_resis[0] == 0:
			pro_transform_b = 'noPRO'
			pro_transform_a = 'noPRO'
			pro_status = False
			pro_resi = pro_resis[0]
		elif pro_resis[0] == tm_i_begin or pro_resis[0] == tm_i_end:
			pro_transform_b = 'PRO is in the end'
			pro_transform_a = 'PRO is in the end'
			pro_status = False
			pro_resi = pro_resis[0]
	
	for n_resi in resis:

		if n_resi > tm_i_begin:
			if n_resi in pro_resis and pro_status and len(pro_resis) == 1:
                		transform_before, transform_after, rmsd = get_transforms(n_resi, tm_i_begin, tm_i_end, tm)
				pro_rmsd = rmsd
				pro_transform_b = transform_before
				pro_transform_a = transform_after
				pro_resi = pro_resis[0]
                
			elif n_resi in pro_resis and pro_status and len(pro_resis) > 1:
                		transform_before, transform_after, rmsd = get_transforms(n_resi, tm_i_begin, tm_i_end, tm)
				pro_rmsds.append(rmsd)
				pro_transforms_b.append(transform_before)
				pro_transforms_a.append(transform_after)
				pro_nums.append(n_resi)
				pro_rmsd = min(pro_rmsds)
                

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
		    continue

		resn_a = list(resns_a)[0]
		resn_i = list(resns_i)[0]
		resns_a.clear()
		resns_i.clear()
		if resn_a != resn_i:
		    continue
		
		rms_cur = cmd.rms_cur('model ' + 'tm' + str(tm) + '_i and resi ' + str(j) + ' and alt A+\"\"', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(tm_a_j) + ' and alt A+\"\"', matchmaker = 1)
		if rms_cur == -1:
			rms_cur = cmd.get_distance('model ' + 'tm' + str(tm) + '_i and resi ' + str(j) + ' and name CA and alt A+\"\"', 'model ' + 'tm' + str(tm) + '_a and resi ' + str(tm_a_j) + ' and name CA and alt A+\"\"')
		
			
		rms_cur_data.append([j, rms_cur])
	
 
        X = [x[0] for x in rms_cur_data]
        Y = [x[1] for x in rms_cur_data]
        clear_rms_array(Y)
        smoothed_data = smooth_data(Y, SMOOTH_STEP)
        min_rms = get_min_rms(smoothed_data)
        x, y = get_min_resi(min_rms, X, Y)

	if x == tm_i_begin:
		x = tm_i_begin + 1
	if x == tm_i_end:
		x == tm_i_end - 1

	
        best_transform_b, best_transform_a, min_rmsd = get_transforms(x, tm_i_begin, tm_i_end, tm)
	
        min_resi = x

	tm_a_j = None
	if tm_a_begin > 1000 and tm_i_begin < 1000:
		tm_a_j = 1000
	elif tm_a_begin < 1000 and tm_i_begin > 1000:
		tm_a_j = -1000
	else:
		tm_a_j = 0

        best_angle_b, best_angle_a = get_angle(X, min_resi, tm, tm_a_j)
        min_resn = resis[min_resi]


	
	if len(pro_resis) > 1:
		for rmsd_it, trans_b_it, trans_a_it, pro_num in zip(pro_rmsds, pro_transforms_b, pro_transforms_a, pro_nums):
			if rmsd_it == pro_rmsd:
				pro_transform_b = trans_b_it
				pro_transform_a = trans_a_it
				pro_resi = pro_num
				
    
	min_resn = resis[min_resi]
    
	return min_rmsd, min_resi, min_resn, best_transform_b, best_transform_a, pro_rmsd, pro_resi, pro_transform_b, pro_transform_a, len(pro_resis), best_angle_b, best_angle_a

error_log = set()
for idx in data['index']:

	i_name, i_begin, i_end, i_chain, a_name, a_begin, a_end, a_chain = get_gpcr_properties(idx)

	cmd.reinitialize()
	cmd.fetch(i_name, async=1)
	cmd.fetch(a_name, async=1)

	time.sleep(2)

	try:
		cmd.create('gpcr_a', 'model ' + a_name + ' and chain ' + a_chain + ' and resi ' + str(a_begin) + '-' + str(a_end) + ' and alt A+\"\"')
		cmd.create('gpcr_i', 'model ' + i_name + ' and chain ' + i_chain + ' and resi ' + str(i_begin) + '-' + str(i_end) + ' and alt A+\"\"')

		cmd.delete(i_name)
		cmd.delete(a_name)

		cmd.align('gpcr_i', 'gpcr_a', transform=1)

	except:
		
		print 'cant create object for row', idx
		break

	for tm in range(1, 8):
		
		try:
			tm_i_begin, tm_i_end, tm_a_begin, tm_a_end = get_tm_borders(idx, tm)

			cmd.create('tm' + str(tm) + '_a', 'model gpcr_a and resi ' + str(tm_a_begin) + '-' + str(tm_a_end))
			cmd.create('tm' + str(tm) + '_i', 'model gpcr_i and resi ' + str(tm_i_begin) + '-' + str(tm_i_end))

			resis = dict()
			min_rmsd, min_resi, min_resn, best_transform_b, best_transform_a, pro_rmsd, pro_resi, pro_transform_b, pro_transform_a, n_of_pro, best_angle_b, best_angle_a = get_best_transforms(tm_i_begin, tm_i_end, tm_a_begin, tm_a_end, tm)
			rmsd_i_a = get_rmsd('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a', 0)
			best_point = min_resn
			best_point_resi = min_resi
			best_axis_b, _ = get_axis_angle(best_transform_b)
			best_axis_a, _ = get_axis_angle(best_transform_a)
			best_rmsd_t_a = min_rmsd
			pro_point_resi = pro_resi
	
			if pro_resi != 0 and pro_transform_b != 'PRO is in the end':
				pro_axis_b, pro_angle_b = get_axis_angle(pro_transform_b)
				pro_axis_a, pro_angle_a = get_axis_angle(pro_transform_a)
				pro_rmsd_t_a = pro_rmsd
			elif pro_transform_b == 'PRO is in the end':
				pro_axis_b, pro_angle_b = 'PRO is in the end', 'PRO is in the end'
				pro_axis_a, pro_angle_a = 'PRO is in the end', 'PRO is in the end'
				pro_rmsd_t_a = 'PRO is in the end'
			else:
				pro_axis_b, pro_angle_b = 'noPRO', 'noPRO'
				pro_axis_a, pro_angle_a = 'noPRO', 'noPRO'
				pro_rmsd_t_a = 'noPRO'
		

			add_gpcr_properties(data, idx, rmsd_i_a, best_point, best_point_resi, best_axis_b, best_angle_b, best_axis_a, best_angle_a, best_rmsd_t_a, pro_point_resi, pro_axis_b, pro_angle_b, pro_axis_a, pro_angle_a, pro_rmsd_t_a, n_of_pro, tm)
			print 'parameters for row ', idx, 'tm', tm, ' were computed'
			cmd.delete('tm' + str(tm) + '_a')
			cmd.delete('tm' + str(tm) + '_i')			
				
		except:
			print 'cant complete on row', idx, 'tm', tm
#			cmd.delete('tm' + str(tm) + '_a')
#			cmd.delete('tm' + str(tm) + '_i')
			error_log.add(idx)

	cmd.delete('gpcr_a')
	cmd.delete('gpcr_i')
	

data.to_csv('GPCR-TM6-table_new_version.csv', index=False)
print 'finish process'
print sorted(list(error_log))

python end
