from scipy.stats import chi2_contingency
import math

file = open("./Kinase_psite_relation.txt",'r')
title = file.readline()
lines = file.readlines()
kinase_site_dic = {}
for line in lines:
    items = line.strip().split("\t")
    sites = items[-1].split()
    sites_new = []
    for item in sites:
        temp = item.strip(",").strip("\'")
        temp_items = temp.split(",")
        this_site = ",".join(temp_items[0:3])
        #print(this_site)
        sites_new.append(this_site)
    kinase_site_dic[items[0]] = sites_new

file_DRP = open("./Phos_DRP_T_N__Matrix.txt",'r')
title = file_DRP.readline()
title_items = title.strip().split("\t")
lines = file_DRP.readlines()
samp_site_exp = {}
for i in range(5,len(title_items)):
    samp_site_exp[title_items[i]] = {}
for line in lines:
    items = line.strip().split("\t")
    new_key = ",".join([items[1],items[0],items[4]])
    for i in range(5,len(items)):
        samp_site_exp[title_items[i]][new_key] = float(items[i])

#samp kinase score
kinase_KNS_samp = {}
for kinase in kinase_site_dic.keys():
    kinase_KNS_samp[kinase] = {}
    sites = kinase_site_dic[kinase]
    for samp in samp_site_exp.keys():
        temp_score = 0
        for site in sites:
            if site in samp_site_exp[samp].keys():
                temp_score += samp_site_exp[samp][site]
        kinase_KNS_samp[kinase][samp] = temp_score

#total chi2
file_experiment = open("./kinase_result.txt",'r')
lines = file_experiment.readlines()
kinase_2_name = {}
for line in lines:
    items = line.strip().split("\t")
    kinase_2_name[items[0]] = items[1]

file_out = open("./Table S5C.txt",'w')
file_out.write("\t".join(["Uniprot ID","Kinase","Activity"])+"\n")
total_T_score = 0
total_N_score = 0
kinase_T_N = {}
for kinase in kinase_KNS_samp.keys():
    kinase_T_N[kinase] = {"T":0,"N":0}
    for samp in kinase_KNS_samp[kinase].keys():
        if "-T" in samp:
            total_T_score += kinase_KNS_samp[kinase][samp]
            kinase_T_N[kinase]["T"] += kinase_KNS_samp[kinase][samp]
        elif "-N" in samp:
            total_N_score += kinase_KNS_samp[kinase][samp]
            kinase_T_N[kinase]["N"] += kinase_KNS_samp[kinase][samp]

for kinase in kinase_T_N.keys():
    table = [[kinase_T_N[kinase]["N"],total_N_score-kinase_T_N[kinase]["N"]],
             [kinase_T_N[kinase]["T"],total_T_score-kinase_T_N[kinase]["T"]]]
    kf = chi2_contingency(table)
    if kinase in kinase_2_name.keys():
        file_out.write("\t".join([kinase,kinase_2_name[kinase],str(kf[1])])+"\n")
    
file_out.close()

#samp chi2
file_class = open(r"./Proteome_Subtype/k_means_out.txt",'r')
lines = file_class.readlines()
samp_class = {}
SI_samp = []
SII_samp = []
SIII_samp = []
for line in lines:
    items = line.strip().split("\t")
    samp_class[items[0]] = items[1]
    if items[1] == "0":
        SII_samp.append(items[0])
    elif items[1] == '1':
        SI_samp.append(items[0])
    elif items[1] == '2':
        SIII_samp.append(items[0])
SI_samp.sort()
SII_samp.sort()
SIII_samp.sort()

file_out = open("./Subtype_Kinase_pvalue.txt",'w')
for samp in samp_site_exp.keys():
    if "-T" in samp:
        samp_T_score = 0
        for kinase in kinase_KNS_samp.keys():
            samp_T_score += kinase_KNS_samp[kinase][samp]
        samp_n = samp.split("-")[0]+"-N"
        samp_N_score = 0
        for kinase in kinase_KNS_samp.keys():
            samp_N_score += kinase_KNS_samp[kinase][samp_n]
        for kinase in kinase_KNS_samp.keys():
            table = [[kinase_KNS_samp[kinase][samp_n],samp_N_score-kinase_KNS_samp[kinase][samp_n]],
                 [kinase_KNS_samp[kinase][samp],samp_T_score-kinase_KNS_samp[kinase][samp]]]
            kf = chi2_contingency(table)
            if samp.split("-")[0] in SI_samp:
                subtype = "SI"
            elif samp.split("-")[0] in SII_samp:
                subtype = "SII"
            elif samp.split("-")[0] in SIII_samp:
                subtype = "SIII"
            file_out.write("\t".join([subtype,samp.split("-")[0],kinase,str(kf[1])])+"\n")
file_out.close()

file_in = open("./Subtype_Kinase_pvalue.txt",'r')
kinase_result = {}
lines = file_in.readlines()
for line in lines:
    items = line.strip().split("\t")
    if items[2] in kinase_2_name.keys():
        if items[2] not in kinase_result.keys():
            kinase_result[items[2]] = []
        kinase_result[items[2]].append("\t".join([items[0],items[1],kinase_2_name[items[2]],items[3]])+"\n")
for kinase in kinase_result.keys():
    kinase_name = kinase_2_name[kinase]
    file_out = open("./Kinase_acitivity_in_Subtype/{0}.txt".format(kinase_name),'w')
    file_out.write("\t".join(["Subtype","Sample","Kinase","Activity"])+"\n")
    for l in kinase_result[kinase]:
        file_out.write(l)
    file_out.close()

