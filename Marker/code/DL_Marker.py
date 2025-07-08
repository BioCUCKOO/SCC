from sklearn.metrics import roc_curve,auc,roc_auc_score,accuracy_score
import numpy as np
import math
import os

def Combination_AUC(c_type):
    #Tumor specific check
    f_list = os.listdir("../Model_Save/{0}".format(c_type))
    ori_combination = []
    for f in f_list:
        if "_data" in f:
            ori_combination.append(f.split("_")[0])

    file_p = open("./{0}_AUC.txt".format(c_type),'w')
    for c in ori_combination:
        data_re_in = open("../Model_Save/{1}/{0}_data.txt".format(c,c_type),'r')
        title = data_re_in.readline()
        lines = data_re_in.readlines()
        matrix_check = []
        label_check = []
        for line in lines:
            items = line.strip().split("\t")
            l = int(items[0])
            values = [float(x) for x in items[1:]]
            matrix_check.append(values)
            label_check.append(l)

        matrix_re_in = open("../Model_Save/{1}/{0}_model.txt".format(c,c_type),'r')
        lines = matrix_re_in.readlines()
        value_items = lines[1].strip().split(",")
        w = [float(x) for x in value_items]
        b = float(lines[-1].strip())
        model_array = [w,b]

        matrix_check_np = np.array(matrix_check)
        re_score_check = []
        for k in range(len(matrix_check_np)):
            score =  np.dot(np.array(model_array[0]).T,np.array(matrix_check_np[k]))+model_array[1]
            sig_score = 1/(1+math.e**(-score))
            re_score_check.append(sig_score)
            
        fpr,tpr,thresholds = roc_curve(label_check,re_score_check,pos_label=1)
        roc_auc_check = auc(fpr,tpr)
        file_p.write("\t".join([c,str(roc_auc_check)])+"\n")
    file_p.close()

Combination_AUC("Tumor_specific")
Combination_AUC("SI_specific")
Combination_AUC("SII_specific")
Combination_AUC("SIII_specific")
