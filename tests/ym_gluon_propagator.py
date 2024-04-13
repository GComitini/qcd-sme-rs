#!/usr/bin/env python3

def read_file_data(file, proj, f0=-0.876):
    data = []
    norm, renpoint = None, 2
    with open(file, 'r') as fin:
        if proj == 'l':
            idx = 4
        elif proj == 't':
            idx = 6
        else:
            idx = 12
        for line in fin:
            ll = line.strip().split()
            if len(ll) != 16:
                continue
            p = float(ll[0])
            # We need to convert between the numerical data, which are computed
            # for fixed Z = 2.684 and F0 = -0.876, and the paper's data, which
            # use the Z and F0 (pi_0(T)) obtained from the lattice.
            dress_inv = 1/(p**2*float(ll[idx])/2.684)+f0+0.876
            prop = 1/(p**2*dress_inv)
            data += [[p, prop]]
            if data[-1][0] == renpoint:
                norm = data[-1][1]
    assert norm is not None
    norm = norm*(renpoint**2)
    for i in range(len(data)):
        data[i][1] = data[i][1]/norm
    return data


transverse_data = read_file_data('T260/Transverse/m450.txt', 't', -0.42)
longitudinal_data = read_file_data('T260/Longitudinal/m425.txt', 'l', -1.42)

print('*** TRANSVERSE PROJECTION ***')
for d in transverse_data:
    print(f'{d[0]}\t{d[1]}')

print('\n*** LONGITUDINAL PROJECTION ***')
for d in longitudinal_data:
    print(f'{d[0]}\t{d[1]}')
