import numpy as np
import matplotlib.pyplot as plt

target_PR_total = 550e3/900

target_PR_total /= np.power(1.1, 6)

HX_PD_per_K = np.linspace(0, 0.425/100, 50)
HX_PD_fixed = 0.01
targ_n_blocks = np.array(list(range(6, 50)))

HX_PD_total = np.linspace(0.01, 0.201, 50)


ESM = []

eta = 0.85
T_in_compr = 250
gamma = 1.31
cp = 870
massflow = 0.05


blocknum = np.zeros((HX_PD_per_K.size, targ_n_blocks.size))
totalweight = np.zeros((HX_PD_per_K.size, targ_n_blocks.size))
PR_final = np.zeros((HX_PD_per_K.size, targ_n_blocks.size))
PD_final = np.zeros((HX_PD_per_K.size, targ_n_blocks.size))

for PD_i, PD_var in enumerate(HX_PD_total):
    for targ_i, target_n_blocks in enumerate(targ_n_blocks):
        n_blocks = target_n_blocks - 1

        PD_crude = PD_var
        required_PR_block = np.power(target_PR_total,1/n_blocks)

        PR_compr = required_PR_block / (1 - PD_crude)

        TR_compr = np.power(PR_compr, (gamma-1)/gamma)
        work_compr = T_in_compr * ((TR_compr - 1)/eta) * cp * massflow
        PD = PD_var #HX_PD_fixed + (T_in_compr * ((TR_compr - 1)/eta) * PD_var)

        PR_compr = 3

        while PR_compr > 2.8: 
            n_blocks += 1
            required_PR_block = np.power(target_PR_total,1/n_blocks)

            PR_compr = required_PR_block / (1 - PD)

            TR_compr = np.power(PR_compr, (gamma-1)/gamma)
            work_compr = T_in_compr * ((TR_compr - 1)/eta) * cp * massflow
            PD = PD_var#HX_PD_fixed + (T_in_compr * ((TR_compr - 1)/eta) * PD_var)
        
        #print(PD_var, target_n_blocks, PR_compr, PD, PR_compr * (1 - PD))

        work_total = work_compr * n_blocks

        ESM_work = work_total * .149
        ESM_cooling = work_total * .121
        ESM_HX_weight = work_total * 0.0815

        blocknum[PD_i, targ_i] = n_blocks
        totalweight[PD_i, targ_i] = ESM_work + ESM_cooling + ESM_HX_weight
        PR_final[PD_i, targ_i] = PR_compr
        PD_final[PD_i, targ_i] = PD

fig, ax = plt.subplots()
plt.contourf(targ_n_blocks, HX_PD_total*100, totalweight, levels=15)
plt.xlabel("Target number of blocks in compressor")
plt.ylabel("Pressure drop in HX, % of static ")
plt.title("Approximate system ESM weight (kg)")
plt.colorbar()

fig, ax = plt.subplots()
plt.contourf(targ_n_blocks, HX_PD_total*100, blocknum, levels=targ_n_blocks[::4])
plt.xlabel("Target number of blocks in compressor")
plt.ylabel("Pressure drop in HX, % of static ")
plt.title("Actual number of blocks in system \n 1% fixed pressure drop plus variable in HX")
plt.colorbar()

fig, ax = plt.subplots()
plt.contourf(targ_n_blocks, HX_PD_total*100, PR_final)
plt.xlabel("Target number of blocks in compressor")
plt.ylabel("Pressure drop in HX, % of static ")
plt.title("Pressure ratio in block compressor \n 1% fixed pressure drop plus variable in HX")
plt.colorbar()

fig, ax = plt.subplots()
plt.contourf(targ_n_blocks, HX_PD_total*100, PD_final*100)
plt.xlabel("Target number of blocks in compressor")
plt.ylabel("Pressure drop in HX, % of static ")
plt.title("Actual pressure drop in HX, % of static pressure")
plt.colorbar()


plt.show()