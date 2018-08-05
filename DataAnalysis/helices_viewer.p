python

import pandas as pd

a_name = '5g53'
i_name = '3vg9'
tm = 5

data = pd.read_csv('GPCR-TM-table-identity-resis-pair.csv', index_col=None)

def get_row_index(data, a_name, i_name):
	for idx in data.index:
		cur_a_name = data['pdb_active'][idx]
		cur_i_name = data['pdb_inactive'][idx]
		print 
		if cur_i_name[:4] == i_name[:4] and cur_a_name[:4] == a_name[:4]:
			return idx	


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

a_begin = int(data['TM1_active'][idx].split('-')[0])
a_end = int(data['TM7_active'][idx].split('-')[1])
i_begin = int(data['TM1_inactive'][idx].split('-')[0])
i_end = int(data['TM7_inactive'][idx].split('-')[1])

cmd.create('gpcr_a', 'model ' + a_name + ' and chain ' + a_chain + ' and resi ' + str(a_begin) + '-' + str(a_end) + ' and alt A+\"\"')
cmd.create('gpcr_i', 'model ' + i_name + ' and chain ' + i_chain + ' and resi ' + str(i_begin) + '-' + str(i_end) + ' and alt A+\"\"')
cmd.align('gpcr_i', 'gpcr_a')

cmd.create('tm' + str(tm) + '_a', 'model gpcr_a and chain ' + a_chain + ' and resi ' + str(tm_a_begin) + '-' + str(tm_a_end) + ' and alt A+\"\"')
cmd.create('tm' + str(tm) + '_i', 'model gpcr_i and chain ' + i_chain + ' and resi ' + str(tm_i_begin) + '-' + str(tm_i_end) + ' and alt A+\"\"')

cmd.delete(i_name)
cmd.delete(a_name)

print 'done'

python end

hide everything
show cartoon
zoom
