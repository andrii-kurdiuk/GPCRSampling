python

import pandas as pd
import numpy as np

resis = set()
data = pd.read_csv('GPCR-TM-table-directions-marked.csv', index_col=None)

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

def get_similar_residues(first_resi_i, last_resi_i, first_resi_a, last_resi_a):

    resi_status_a = 0
    resi_status_i = 0
    if first_resi_a > 1000:
        first_resi_a_fixed = first_resi_a - 1000
        last_resi_a_fixed = last_resi_a - 1000
        resi_status_a = 1
    else:
        first_resi_a_fixed = first_resi_a
        last_resi_a_fixed = last_resi_a

    if first_resi_i > 1000:
        first_resi_i_fixed = first_resi_i - 1000
        last_resi_i_fixed = last_resi_i - 1000
        resi_status_i = 1
    else:
        first_resi_i_fixed = first_resi_i
        last_resi_i_fixed = last_resi_i

    first_resi_uni = 0
    last_resi_uni = 0

    if first_resi_a_fixed > first_resi_i_fixed:
        first_resi_uni = first_resi_i_fixed
    elif first_resi_a_fixed < first_resi_i_fixed:
        first_resi_uni = first_resi_a_fixed
    else:
        first_resi_uni = first_resi_a_fixed

    if last_resi_a_fixed > last_resi_i_fixed:
        last_resi_uni = last_resi_i_fixed
    elif last_resi_a_fixed < last_resi_i_fixed:
        last_resi_uni = last_resi_a_fixed
    else:
        last_resi_uni = last_resi_a_fixed

    return first_resi_uni + resi_status_i * 1000, last_resi_uni + resi_status_i * 1000, first_resi_uni + resi_status_a * 1000, last_resi_uni + resi_status_a * 1000

def get_top_projections(model_tm_i, model_tm_a):
    first_resi_i, last_resi_i = find_borders_resi(model_tm_i)
    first_resi_a, last_resi_a = find_borders_resi(model_tm_a)

    first_i_uni, last_i_uni, first_a_uni, last_a_uni = get_similar_residues(first_resi_i, last_resi_i, first_resi_a, last_resi_a)

    first_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(first_i_uni) + ' and name CA', 1)
    last_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(last_i_uni) + ' and name CA', 1)

    first_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(first_a_uni) + ' and name CA', 1)
    last_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(last_a_uni) + ' and name CA', 1)

    delta_none = 0
    while np.array([first_i_coord is None, last_i_coord is None, first_a_coord is None, last_a_coord is None]).any():
        delta_none += 1
        first_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(first_i_uni + delta_none) + ' and name CA', 1)
        last_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(last_i_uni - delta_none) + ' and name CA', 1)

        first_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(first_a_uni + delta_none) + ' and name CA', 1)
        last_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(last_a_uni - delta_none) + ' and name CA', 1)

    if C_top != 0.0:
        reference_point = np.array([0.0, 0.0, (-1.) * D_top / C_top])
    elif B_top != 0.0:
        reference_point = np.array([0.0, (-1.) * D_top / B_top, 0.0])
    elif A_top != 0.0:
        reference_point = np.array([(-1.) * D_top / A_top, 0.0, 0.0])
    else:
        reference_point = np.array([0.0, 0.0, 0.0])

    alpha_first_i = np.dot(first_i_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_first_i = first_i_coord - alpha_first_i * n_top

    alpha_last_i = np.dot(last_i_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_last_i = last_i_coord - alpha_last_i * n_top

    alpha_first_a = np.dot(first_a_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_first_a = first_a_coord - alpha_first_a * n_top

    alpha_last_a = np.dot(last_a_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_last_a = last_a_coord - alpha_last_a * n_top

    return projection_first_i, projection_last_i, projection_first_a, projection_last_a

def get_bottom_projections(model_tm_i, model_tm_a):
    first_resi_i, last_resi_i = find_borders_resi(model_tm_i)
    first_resi_a, last_resi_a = find_borders_resi(model_tm_a)

    first_i_uni, last_i_uni, first_a_uni, last_a_uni = get_similar_residues(first_resi_i, last_resi_i, first_resi_a, last_resi_a)

    first_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(first_i_uni) + ' and name CA', 1)
    last_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(last_i_uni) + ' and name CA', 1)

    first_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(first_a_uni) + ' and name CA', 1)
    last_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(last_a_uni) + ' and name CA', 1)

    delta_none = 0
    while np.array([first_i_coord is None, last_i_coord is None, first_a_coord is None, last_a_coord is None]).any():
        delta_none += 1
        first_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(first_i_uni + delta_none) + ' and name CA', 1)
        last_i_coord = cmd.get_coords('model ' + model_tm_i + ' and resi ' + str(last_i_uni - delta_none) + ' and name CA', 1)

        first_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(first_a_uni + delta_none) + ' and name CA', 1)
        last_a_coord = cmd.get_coords('model ' + model_tm_a + ' and resi ' + str(last_a_uni - delta_none) + ' and name CA', 1)


    if C_bottom != 0.0:
        reference_point = np.array([0.0, 0.0, (-1.) * D_bottom / C_bottom])
    elif B_bottom != 0.0:
        reference_point = np.array([0.0, (-1.) * D_bottom / B_bottom, 0.0])
    elif A_bottom != 0.0:
        reference_point = np.array([(-1.) * D_bottom / A_bottom, 0.0, 0.0])
    else:
        reference_point = np.array([0.0, 0.0, 0.0])

    alpha_first_i = np.dot(first_i_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_first_i = first_i_coord - alpha_first_i * n_top

    alpha_last_i = np.dot(last_i_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_last_i = last_i_coord - alpha_last_i * n_top

    alpha_first_a = np.dot(first_a_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_first_a = first_a_coord - alpha_first_a * n_top

    alpha_last_a = np.dot(last_a_coord - reference_point, n_top) / np.dot(n_top, n_top)
    projection_last_a = last_a_coord - alpha_last_a * n_top

    return projection_first_i, projection_last_i, projection_first_a, projection_last_a


def add_gpcr_properties(data, idx, tm, active_bottom_first, active_top_first, active_bottom_last, active_top_last, inactive_bottom_first, inactive_top_first, inactive_bottom_last, inactive_top_last):

    tm_name = 'TM' + str(tm) + '_'
    data[tm_name + 'active_bottom_first'][idx] = get_str_axis(active_bottom_first[0])
    data[tm_name + 'active_top_first'][idx] = get_str_axis(active_top_first[0])
    data[tm_name + 'active_bottom_last'][idx] = get_str_axis(active_bottom_last[0])
    data[tm_name + 'active_top_last'][idx] = get_str_axis(active_top_last[0])
    data[tm_name + 'inactive_bottom_first'][idx] = get_str_axis(inactive_bottom_first[0])
    data[tm_name + 'inactive_top_first'][idx] = get_str_axis(inactive_top_first[0])
    data[tm_name + 'inactive_bottom_last'][idx] = get_str_axis(inactive_bottom_last[0])
    data[tm_name + 'inactive_top_last'][idx] = get_str_axis(inactive_top_last[0])

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

        i_bottom_first, i_bottom_last, a_bottom_first, a_bottom_last = get_bottom_projections('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a')
        i_top_first, i_top_last, a_top_first, a_top_last = get_top_projections('tm' + str(tm) + '_i', 'tm' + str(tm) + '_a')

        add_gpcr_properties(data, idx, tm, a_bottom_first, a_top_first, a_bottom_last, a_top_last, i_bottom_first, i_top_first, i_bottom_last, i_top_last)

        cmd.delete('tm' + str(tm) + '_a')
        cmd.delete('tm' + str(tm) + '_i')

        print 'done for idx', idx, 'tm', tm

    cmd.delete('gpcr_a')
    cmd.delete('gpcr_i')

data.to_csv('GPCR-TM-table-directions-fixed.csv', index=False)
print 'finished process'
python end
