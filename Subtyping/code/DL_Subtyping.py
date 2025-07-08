from  scipy.stats import chi2_contingency
from tensorflow.keras.models import Model
import numpy as np
import keras as K

def Load_Model():
    CF_model = K.models.load_model("../CF_Model/Ori_Model/DNN_Best.model")
    Proteome_model = K.models.load_model("../Proteome_Model/Ori_Model/DNN_Best.model")
    Phos_model = K.models.load_model("../Phosphoproteome_Model/Ori_Model/DNN_Best.model")
    RNA_model = K.models.load_model("../Transcriptome_Model/Ori_Model/DNN_Best.model")

    layer_name = "encoder"
    CF_encoder = Model(inputs=CF_model.input,outputs=CF_model.get_layer(layer_name).output)
    Proteome_encoder = Model(inputs=Proteome_model.input,outputs=Proteome_model.get_layer(layer_name).output)
    Phos_encoder = Model(inputs=Phos_model.input,outputs=Phos_model.get_layer(layer_name).output)
    RNA_encoder = Model(inputs=RNA_model.input,outputs=RNA_model.get_layer(layer_name).output)

    return(CF_encoder,Proteome_encoder,Phos_encoder,RNA_encoder)

def load_data(path):
    file = open(path)
    title = file.readline()
    lines = file.readlines()
    samp_list = []
    matrix = []
    for line in lines:
        items = line.strip().split("\t")
        samp_list.append(items[0])
        matrix.append([float(x) for x in items[1:]])
    return(samp_list,matrix)

def cluster_in():
    CF_encoder,Proteome_encoder,Phos_encoder,RNA_encoder = Load_Model()
    #CF data
    cf_samp,cf_matrix = load_data("../CF_Model/CF_keras_in.txt")
    cf_10 = CF_encoder.predict(np.array(cf_matrix))
    file_out = open("./CF_Encoder_Feature.txt",'w')
    for i in range(len(cf_samp)):
        wrt_str = [cf_samp[i],",".join([str(x) for x in cf_10[i]])]
        file_out.write("\t".join(wrt_str)+"\n")
    file_out.close()

    #Proteome data
    proteome_samp,proteome_matrix = load_data("../Proteome_Model/Proteome_keras_in.txt")
    proteome_10 = Proteome_encoder.predict(np.array(proteome_matrix))
    file_out = open("./Proteome_Encoder_Feature.txt",'w')
    for i in range(len(proteome_samp)):
        wrt_str = [proteome_samp[i],",".join([str(x) for x in proteome_10[i]])]
        file_out.write("\t".join(wrt_str)+"\n")
    file_out.close()

    #Phosphoproteome data
    phos_samp,phos_matrix = load_data("../Phosphoproteome_Model/Phos_keras_in.txt")
    phos_10 = Phos_encoder.predict(np.array(phos_matrix))
    file_out = open("./Phosphoproteome_Encoder_Feature.txt",'w')
    for i in range(len(phos_samp)):
        wrt_str = [phos_samp[i],",".join([str(x) for x in phos_10[i]])]
        file_out.write("\t".join(wrt_str)+"\n")
    file_out.close()

    #RNA data
    RNA_samp,RNA_matrix = load_data("../Transcriptome_Model/RNA_keras_in.txt")
    rna_10 = RNA_encoder.predict(np.array(RNA_matrix))
    file_out = open("./Transcriptome_Encoder_Feature.txt",'w')
    for i in range(len(RNA_samp)):
        wrt_str = [RNA_samp[i],",".join([str(x) for x in rna_10[i]])]
        file_out.write("\t".join(wrt_str)+"\n")
    file_out.close()

def CF_Chi2_Check():
    #load cf model
    file_cf = open(r"../CF_Model/Ori_Subtype/k_means_out.txt",'r')
    lines = file_cf.readlines()
    subtype_dic = {}
    for line in lines:
        items = line.strip().split("\t")
        subtype_dic[items[0]] = items[1]

    #juge protein group and now subtype relation
    file_protein = open(r"../Proteome_Model/Ori_Subtype/k_means_out.txt",'r')
    lines = file_protein.readlines()
    protein_subtype = {}
    for line in lines:
        items = line.strip().split("\t")
        protein_subtype[items[0]] = items[1]
    
    overlap_dic = {'0':[0,0,0],"1":[0,0,0],"2":[0,0,0]}
    for samp in subtype_dic.keys():
        if samp in protein_subtype.keys():
            pre_subtype_index = int(protein_subtype[samp])
            overlap_dic[subtype_dic[samp]][pre_subtype_index] += 1
    
    chi_matrix = np.array([overlap_dic["0"],overlap_dic["1"],overlap_dic["2"]])
    kf = chi2_contingency(chi_matrix)
    pvalue1 = kf[1]

    #juge phos group and now subtype relation
    file_phos = open(r"../Phosphoproteome_Model/Ori_Subtype/k_means_out.txt",'r')
    lines = file_phos.readlines()
    phos_subtype = {}
    for line in lines:
        items = line.strip().split("\t")
        phos_subtype[items[0]] = items[1]
    
    overlap_dic = {'0':[0,0,0],"1":[0,0,0],"2":[0,0,0]}
    for samp in subtype_dic.keys():
        if samp in phos_subtype.keys():
            pre_subtype_index = int(phos_subtype[samp])
            overlap_dic[subtype_dic[samp]][pre_subtype_index] += 1
    
    chi_matrix = np.array([overlap_dic["0"],overlap_dic["1"],overlap_dic["2"]])
    kf = chi2_contingency(chi_matrix)
    pvalue2 = kf[1]

    #juge RNA group and now subtype relation
    file_rna = open(r"../Transcriptome_Model/Ori_Subtype/k_means_out.txt",'r')
    lines = file_rna.readlines()
    RNA_subtype = {}
    for line in lines:
        items = line.strip().split("\t")
        RNA_subtype[items[0]] = items[1]
    
    overlap_dic = {'0':[0,0,0],"1":[0,0,0],"2":[0,0,0]}
    for samp in subtype_dic.keys():
        if samp in RNA_subtype.keys():
            pre_subtype_index = int(RNA_subtype[samp])
            overlap_dic[subtype_dic[samp]][pre_subtype_index] += 1
    
    chi_matrix = np.array([overlap_dic["0"],overlap_dic["1"],overlap_dic["2"]])
    kf = chi2_contingency(chi_matrix)
    pvalue3 = kf[1]

    file_p_out = open("./P_Calculate.txt",'w')
    file_p_out.write("\t".join(["Proteome",str(pvalue1)])+"\n")
    file_p_out.write("\t".join(["Phos",str(pvalue2)])+"\n")
    file_p_out.write("\t".join(["RNA",str(pvalue3)])+"\n")
    file_p_out.close()

cluster_in()
CF_Chi2_Check()
