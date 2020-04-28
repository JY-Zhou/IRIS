import sys, os
import numpy as np

os.system('rm -f ./cov_files/*')

for file in os.listdir('./sto_files'):
    RF_id = file.split('.')[0]
    
    print(RF_id)
    
    # Run R-scape
    os.system('../../lib/R-scape --naive --nofigures --MIr ./sto_files/' + file + ' > /dev/null')
    with open('./sto_files/' + file, 'rb') as f:
        content = f.readlines()
        RF_name = content[3].decode().split(' ')[-1].strip()
        msa_len = len(content[-2].decode().split(' ')[-1].strip())
        
    res_mat = np.zeros((msa_len, msa_len))
    
    with open(RF_id + '_' + RF_name + '.cov') as f:
        for line in f:
            if not '#' in line:
                substr = line.split('\t')
                l = int(substr[1].strip())-1
                r = int(substr[2].strip())-1
                if l >= r:
                    raise ValueError('Illegal position')
                score = float(substr[3].strip())
                res_mat[l, r] = score
                res_mat[r, l] = score
        
    with open('./cov_files/' + RF_id + '.cov', 'w') as f:
        for i in range(msa_len):
            for j in range(msa_len):
                print(res_mat[i, j], file = f, end = '\t')
            print('', file = f)

    os.remove(RF_id + '_' + RF_name + '.cov')
    os.remove(RF_id + '_' + RF_name + '.power')
    os.remove(RF_id + '_' + RF_name + '.sorted.cov')