# -*- coding: utf-8 -*-
# @Time    : 2023/4/10  14:03
# @Author  : Gou Yujie
# @File    : classify.py
from tensorflow.keras.models import Model
from keras.utils import to_categorical
from sklearn import metrics
from keras.callbacks import Callback,ModelCheckpoint
import cv2
import numpy as np
import os
import matplotlib.pyplot as plt
from tensorflow.keras.metrics import AUC
from keras.models import load_model
from PIL import Image

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

def load_test():
    image=np.load('recognize_image.npy')
    label=np.load('recognize_types.npy')
    return image,label

images,label=load_test()
cellmodel=load_model("F:/cervical_cancer/pic_process/models/classify_filt_canimu_1_5.model")
x_pred=cellmodel.predict(images, batch_size=100)
x_pred=np.array(x_pred)
y_test = np.array(to_categorical(types))
for i in range(3):
    AUC = roc_auc_score(np.array(y_test[:,i]).reshape(-1, 1), np.array(x_pred[:,i]).reshape(-1, 1))
    print(AUC)
