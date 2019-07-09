import numpy as np
import json

def list_append(lst, item):
    lst.append(item)
    return lst

def parse_json(filen):
    with open(filen) as fhandle:
        data = json.load(fhandle)

    for key in data:
        old_data = data[key]
        new_data = []
        for i in range(1,len(old_data)-1):
            if i % 2 == 1:
                new_data.append('${} \pm {}$'.format(np.round(old_data[i], 2),
                                                     np.round(old_data[i+1], 2)))

        data[key] = new_data
    return data

files=[
    'domain_sizes_mean_dma.json',
    'domain_sizes_mean_whitlam.json',
    'domain_sizes_mean_dmoc.json']
names_dict = {}
names_dict[files[0]] = ('dma-c',
                        'dma-as2')
names_dict[files[1]] = ('w. et al.',)
names_dict[files[2]] = ('dmo-c',)

param_dict = {}
param_dict[files[0]] = ('0.5', '0.4')
param_dict[files[1]] = ('0.5',)
param_dict[files[2]] = ('0.5',)

# this is the results table 
table_dict = {}
for fi in files:
    raw_data = parse_json(fi)
    for i, param in enumerate(param_dict[fi]):
        table_dict[names_dict[fi][i]] = raw_data[param]

print(table_dict)
print()
print()

names = [
    '$\sigma_{p}$',
    '$\sigma_{np}$',
    '$\\text{frac}_{np}(\\text{domains})$',
    '$\\text{frac}_{np}(\\text{bonds})$']

# generate latex table 
def generate_table(n_lines,n_col, table_dict, names):
    start_centering='\\begin{center}'
    end_centering='\\end{center}'

    column_format='|l |'
    for i in range(n_col): column_format = column_format + 'l |'
    start_tabular='\\begin{{tabular}}{{{}}}'.format(column_format)
    end_tabular='\\end{tabular}'
    hline='\\hline'
    breakline='  \\\\  '
    line_legend= 'rhombus type  '
    for name in names: line_legend = line_legend + ' & {} '.format(name)
    line_legend = line_legend + breakline + hline 
    col_legend=''
    for key in table_dict.keys(): col_legend = col_legend + '& {}'.format(key)

    # print table
    print(start_centering)
    print(start_tabular)
    print(hline)
    print(line_legend)
    for key in table_dict:
        value_string = ''
        for value in table_dict[key]:
            value_string = value_string + ' & {} '.format(value)
        line = '{} '.format(key) + value_string + breakline + hline
        print(line)
    print(end_tabular)
    print(end_centering)

generate_table(4,4, table_dict, names)
