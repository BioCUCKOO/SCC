# -*- coding: utf-8 -*-
# @Time    : 2023/4/10  17:26
# @Author  : Gou Yujie
# @File    : unet_sample.py
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam,RMSprop
from tensorflow.keras.layers import Conv2D,Dense,Flatten,Dropout,MaxPool2D,BatchNormalization,Activation,Input,Conv2DTranspose,Concatenate,concatenate,UpSampling2D,AveragePooling2D,SeparableConv2D,MaxPooling2D
from sklearn.model_selection import train_test_split
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
from sklearn.model_selection import KFold,StratifiedKFold
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

image=np.load('recognize_image.npy')
unet=load_model("get_unet_half_add4.h5")
out=unet.predict(image)
for pred in list(out):
    pred=pred[:,:,1].reshape(64,64)
    pred=np.uint8(pred*255)
    binary = cv2.adaptiveThreshold(pred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 31, -1)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    img_erode = cv2.erode(binary, kernel, iterations=3)
    si = cv2.medianBlur(img_erode, 7)
    si2 = cv2.medianBlur(si, 3)
    img_dilated = cv2.dilate(si2, kernel, iterations=3)
    contours, hierarchy = cv2.findContours(img_dilated, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    for k in range(len(contours)):
        M2 = cv2.moments(contours[k])
        if M2["m00"] != 0:
            center_x2 = int(M2["m10"] / M2["m00"])
            center_y2 = int(M2["m01"] / M2["m00"])
            if (contours[k].size > 15):
                cv2.circle(image, (center_x2, center_y2), 2, (0, 255, 0), -1)


