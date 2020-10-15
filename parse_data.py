# --------------------------------------------------------------
# TERMINAL VERSION: CODE FOR PARSING DATA FROM THE OUTPUT FILES
# --------------------------------------------------------------

import sys

import sys

# Stores the file name into a variable and then opens that file.
filename = sys.argv[1]
with open('{0}'.format(filename), 'r') as inp:
    
    anharm_raw = []
    copy = False
    for line in inp:
        line = line.strip()
        if line == 'Integrated intensity (I) in km.mol^-1':
            copy = True
        elif line == 'Dipole strengths (DS) in 10^-40 esu^2.cm^2':
            copy = False
        elif copy:
            anharm_raw.append(line)
            
    anharm_raw = list(filter(None, anharm_raw))
    anharm_raw = list(filter(lambda line: not line.startswith('-'), anharm_raw))
    
    # Gets fundamentals info.
    fundamentals = []
    for row in anharm_raw:
        row = row.strip()
        if row == 'Fundamental Bands':
            copy = True
        elif row == 'Overtones':
            copy = False
        elif copy:
            fundamentals.append(row)

    # Breaks each component in the fundamentals list to get the separated info from the modes, harmonic and anharmonic
    # freqs and harmonic and anharmonic intensities.
    modes_fund = []
    modes_fund_kind = []
    freqs_fund = []
    ints_fund = []
    for line in fundamentals:
        splitted_line = line.split()
        if 'Mode(n,l)' in splitted_line:
            modes_fund_index = splitted_line.index('Mode(n,l)')
        elif 'Mode(n)' in splitted_line:
            modes_fund_index = splitted_line.index('Mode(n)')
        else:
            modes_fund.append(line.split()[modes_fund_index])
            modes_fund.append(line.split()[modes_fund_index])
            
        if 'E(harm)' in splitted_line:
            freqs_fund_index = splitted_line.index('E(harm)')
        else:
            freqs_fund.append(line.split()[freqs_fund_index])
            modes_fund_kind.append('FundHarmonic')
            
        if 'E(anharm)' in splitted_line:
            anfreqs_fund_index = splitted_line.index('E(anharm)')
        else:
            freqs_fund.append(line.split()[anfreqs_fund_index])
            modes_fund_kind.append('FundAnharm')
            
        if 'I(harm)' in splitted_line:
            ints_fund_index = splitted_line.index('I(harm)')
        else:
            ints_fund.append(line.split()[ints_fund_index])
            
        if 'I(anharm)' in splitted_line:
            anints_fund_index = splitted_line.index('I(anharm)')
        else:
            ints_fund.append(line.split()[anints_fund_index])

    # Gets scaled fundamentals (scaling factor for wB97X-D/def2-SVPD is 0.9542 ref. 15KeBrMa)
    modes_fund_kind_scaled = []
    modes_fund_scaled = []
    freqs_fund_scaled = []
    ints_fund_scaled = []
    for kind_s,mode_s,harm_freq_s,harm_int_s in zip(modes_fund_kind,modes_fund,freqs_fund,ints_fund):
        if kind_s == 'FundHarmonic':
            modes_fund_kind_scaled.append('FundScaled')
            modes_fund_scaled.append(mode_s)
            freqs_fund_scaled.append(str(float(harm_freq_s)*0.9542))
            ints_fund_scaled.append(harm_int_s)
        
    # Gets overtones info.
    overtones = []
    for row in anharm_raw:
        if row == 'Overtones':
            copy = True
        elif row == 'Combination Bands':
            copy = False
        elif copy:
            overtones.append(row)
    
    # Breaks overtones info to get the modes and anharmonic frequencies and intensities.
    modes_overt = []
    modes_overt_kind = []
    freqs_overt = []
    ints_overt = []
    for line in overtones:
        splitted_line = line.split()
        if 'Mode(n,l)' in splitted_line:
            modes_overt_index = splitted_line.index('Mode(n,l)')
        elif 'Mode(n)' in splitted_line:
            modes_overt_index = splitted_line.index('Mode(n)')
        else:
            modes_overt.append(line.split()[modes_overt_index])
            modes_overt_kind.append('Overtone')
            
        if 'E(anharm)' in splitted_line:
            freqs_overt_index = splitted_line.index('E(anharm)')
        else:
            freqs_overt.append(line.split()[freqs_overt_index])
            
        if 'I(anharm)' in splitted_line:
            ints_overt_index = splitted_line.index('I(anharm)')
        else:
            ints_overt.append(line.split()[ints_overt_index])
            
    # Gets combination bands info.
    comb_bands = []
    for row in anharm_raw:
        if row == 'Combination Bands':
            copy = True
        elif row == 'Units: Transition energies (E) in cm^-1':
            copy = False
        elif copy:
            comb_bands.append(row)
        
    # Adds a new heading to the combination bands info to make it readable. The change in the heading is only
    # on the second term that is changed from Mode(n,l) to ModeC(n,l).
    comb_bands = list(filter(lambda line: not line.startswith('M'), comb_bands))
    comb_bands.insert(0, 'Mode(n,l)   ModeC(n,l)    E(harm)   E(anharm)                      I(anharm)')
    
    # Gets the modes and anharmonic frequencies and intensities info from the combination bands list.
    modes_cb_1 = []
    modes_cb_2 = []
    modes_cb_tot = []
    modes_cb_kind = []
    freqs_cb = []
    ints_cb = []
    for line in comb_bands:
        splitted_line = line.split()
        if 'Mode(n,l)' in splitted_line:
            modes_cb_1_index = splitted_line.index('Mode(n,l)')
        elif 'Mode(n)' in splitted_line:
            modes_cb_1_index = splitted_line.index('Mode(n)')
        else:
            modes_cb_1.append(line.split()[modes_cb_1_index])
            
        if 'ModeC(n,l)' in splitted_line:
            modes_cb_2_index = splitted_line.index('ModeC(n,l)')
        else:
            modes_cb_2.append(line.split()[modes_cb_2_index])
            
        if 'E(anharm)' in splitted_line:
            freqs_cb_index = splitted_line.index('E(anharm)')
        else:
            freqs_cb.append(line.split()[freqs_cb_index])
            
        if 'I(anharm)' in splitted_line:
            ints_cb_index = splitted_line.index('I(anharm)')
        else:
            ints_cb.append(line.split()[ints_cb_index])
            
    for i,j in zip(modes_cb_1, modes_cb_2):
        modes_cb_tot.append(i+j)
        
    for element in modes_cb_tot:
        modes_cb_kind.append('CombBand')
    
# Gets all the information together
all_kinds = modes_fund_kind + modes_fund_kind_scaled + modes_overt_kind + modes_cb_kind
all_modes = modes_fund + modes_fund_scaled + modes_overt + modes_cb_tot
all_freqs = freqs_fund + freqs_fund_scaled + freqs_overt + freqs_cb
all_ints = ints_fund + ints_fund_scaled + ints_overt + ints_cb

# Removes the comma form the modes notation and creates a list with the molecular formula for all
# frequencies and intensities being considered.
list_filename = filename.split('.')
all_modes_nocomma = []
mol_id = []
for i in all_modes:
    all_modes_nocomma.append(i.replace(',', '_'))
    mol_id.append(list_filename[0])
    
# Prints a CSV file sotring all information.
for mol,kind,mode,freqs,ints in zip(mol_id,all_kinds,all_modes_nocomma,all_freqs,all_ints):
    print(mol+','+kind+','+mode+','+freqs+','+ints, file = open('pmol_vibinfo.csv','a'))
