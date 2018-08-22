python

import pandas as pd
import numpy as np

resis = set()
data = pd.read_csv('GPCR-TM-table-directions.csv', index_col=None)

def get_str_axis(axis):
    return str(axis[0]) + ', ' + str(axis[1]) + ', ' + str(axis[2])

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

def find_borders_resi(model_tm):
    cmd.iterate('model ' + model_tm, 'resis.add(int(resi))')

    list_resis = sorted(list(resis))
    first_resi = list_resis[0]
    last_resi = list_resis[-1]
    resis.clear()

    return first_resi, last_resi

def get_plane_equation(p1, p2, p3, side):
    v1 = np.array(cmd.get_atom_coords('model ' + side + ' and resi ' + str(p1), 1))
    v2 = np.array(cmd.get_atom_coords('model ' + side + ' and resi ' + str(p2), 1))
    v3 = np.array(cmd.get_atom_coords('model ' + side + ' and resi ' + str(p3), 1))

    v_plane_1 = v3 - v1
    v_plane_2 = v2 - v1

    A, B, C = np.cross(v_plane_1, v_plane_2)
    D = (-1.) * np.dot(np.cross(v_plane_1, v_plane_2), v3)

    return A, B, C, D

def get_top_projections(model_tm):
    first_resi, last_resi = find_borders_resi(model_tm)
    first_coord = np.array(cmd.get_atom_coords('model ' + model_tm + ' and resi ' + str(first_resi) + ' and name CA', 1))
    last_coord = np.array(cmd.get_atom_coords('model ' + model_tm + ' and resi ' + str(last_resi) + ' and name CA', 1))

    if C_top != 0.0:
        reference_point = np.array([0.0, 0.0, (-1.) * D_top / C_top])
    elif B_top != 0.0:
        reference_point = np.array([0.0, (-1.) * D_top / B_top, 0.0])
    elif A_top != 0.0:
        reference_point = np.array([(-1.) * D_top / A_top, 0.0, 0.0])
    else:
        reference_point = np.array([0.0, 0.0, 0.0])

    alpha_first = np.dot(first_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_first = first_coord - alpha_first * n_top

    alpha_last = np.dot(last_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_last = last_coord - alpha_last * n_top

    return projection_first, projection_last

def get_bottom_projections(model_tm):
    first_resi, last_resi = find_borders_resi(model_tm)
    first_coord = np.array(cmd.get_atom_coords('model ' + model_tm + ' and resi ' + str(first_resi) + ' and name CA', 1))
    last_coord = np.array(cmd.get_atom_coords('model ' + model_tm + ' and resi ' + str(last_resi) + ' and name CA', 1))

    if C_bottom != 0.0:
        reference_point = np.array([0.0, 0.0, (-1.) * D_bottom / C_bottom])
    elif B_bottom != 0.0:
        reference_point = np.array([0.0, (-1.) * D_bottom / B_bottom, 0.0])
    elif A_bottom != 0.0:
        reference_point = np.array([(-1.) * D_bottom / A_bottom, 0.0, 0.0])
    else:
        reference_point = np.array([0.0, 0.0, 0.0])

    alpha_first = np.dot(first_coord - reference_point, n_bottom) / np.dot(n_bottom, n_bottom)
    projection_first = first_coord - alpha_first * n_bottom

    alpha_last = np.dot(last_coord - reference_point, n_bottom) / np.dot(n_bottom, n_bottom)
    projection_last = last_coord - alpha_last * n_bottom

    return projection_first, projection_last

def add_gpcr_properties(data, idx, tm, active_bottom_first, active_top_first, active_bottom_last, active_top_last, inactive_bottom_first, inactive_top_first, inactive_bottom_last, inactive_top_last):

    tm_name = 'TM' + str(tm) + '_'
    data[tm_name + 'active_bottom_first'][idx] = get_str_axis(active_bottom_first)
    data[tm_name + 'active_top_first'][idx] = get_str_axis(active_top_first)
    data[tm_name + 'active_bottom_last'][idx] = get_str_axis(active_bottom_last)
    data[tm_name + 'active_top_last'][idx] = get_str_axis(active_top_last)
    data[tm_name + 'inactive_bottom_first'][idx] = get_str_axis(inactive_bottom_first)
    data[tm_name + 'inactive_top_first'][idx] = get_str_axis(inactive_top_first)
    data[tm_name + 'inactive_bottom_last'][idx] = get_str_axis(inactive_bottom_last)
    data[tm_name + 'inactive_top_last'][idx] = get_str_axis(inactive_top_last)

cmd.reinitialize()

cmd.load('4eiy_layers.pdb')
cmd.create('reference', 'model 4eiy_layers and resi 2-291 and alt A+\"\"')
cmd.create('top', 'model 4eiy_layers and resn DUM and name O')
cmd.create('bottom', 'model 4eiy_layers and resn DUM and name N')
cmd.delete('4eiy_layers')

A_top, B_top, C_top, D_top =  get_plane_equation(2459, 2460, 2147, 'top')
n_top = np.array([A_top, B_top, C_top])
A_bottom, B_bottom, C_bottom, D_bottom = get_plane_equation(2459, 2460, 2147, 'bottom')
n_bottom = np.array([A_bottom, B_bottom, C_bottom])

cmd.delete('top')
cmd.delete('bottom')

for idx in data['index']:

    i_name, i_begin, i_end, i_chain, a_name, a_begin, a_end, a_chain = get_gpcr_properties(idx)

    cmd.fetch(i_name, async=1)
    cmd.fetch(a_name, async=1)
    time.sleep(2)

    cmd.create('gpcr_a', 'model ' + a_name + ' and chain ' + a_chain + ' and resi ' + str(a_begin) + '-' + str(a_end) + ' and alt A+\"\"')
    cmd.create('gpcr_i', 'model ' + i_name + ' and chain ' + i_chain + ' and resi ' + str(i_begin) + '-' + str(i_end) + ' and alt A+\"\"')

    cmd.delete(i_name)
    cmd.delete(a_name)

    cmd.super('gpcr_a', 'reference')
    cmd.super('gpcr_i', 'gpcr_a')

    for tm in range(1, 8):
        tm_i_begin, tm_i_end, tm_a_begin, tm_a_end = get_tm_borders(idx, tm)
        cmd.create('tm' + str(tm) + '_a', 'model gpcr_a and resi ' + str(tm_a_begin) + '-' + str(tm_a_end) + ' and alt A+\"\"')
        cmd.create('tm' + str(tm) + '_i', 'model gpcr_i and resi ' + str(tm_i_begin) + '-' + str(tm_i_end) + ' and alt A+\"\"')

        active_bottom_first, active_bottom_last = get_bottom_projections('tm' + str(tm) + '_a')
        active_top_first, active_top_last = get_top_projections('tm' + str(tm) + '_a')

        inactive_bottom_first, inactive_bottom_last = get_bottom_projections('tm' + str(tm) + '_i')
        inactive_top_first, inactive_top_last = get_top_projections('tm' + str(tm) + '_i')

        add_gpcr_properties(data, idx, tm, active_bottom_first, active_top_first, active_bottom_last, active_top_last, inactive_bottom_first, inactive_top_first, inactive_bottom_last, inactive_top_last)

        cmd.delete('tm' + str(tm) + '_a')
        cmd.delete('tm' + str(tm) + '_i')

        print 'done for idx', idx, 'tm', tm

    cmd.delete('gpcr_a')
    cmd.delete('gpcr_i')

data.to_csv('GPCR-TM-table-directions.csv', index=False)
print 'finished process'
python end
