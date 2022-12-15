from __future__ import print_function



import os

# ran without problems on fermilab gpu server singularity environment
# (tensorflow 2)

import tensorflow as tf
if tf.__version__.startswith("2."):
  tf = tf.compat.v1
tf.disable_eager_execution()

import keras


import math

import sys

import argparse


from keras.callbacks import ModelCheckpoint

from  matplotlib import pyplot as plt

import pylab

import glob

from tensorflow.python.framework import ops

from tensorflow.python.ops import clip_ops

from tensorflow.python.ops import math_ops

from tensorflow.python.ops import nn

from numpy import concatenate as concatenatenp

import random



from keras.layers import AlphaDropout



from keras import backend as K

from keras.models import load_model

# adjust input file (h5) path here
weight_file_path = '/storage/local/data1/gpuscratch/hichemb/DeepCore_git/DeepCore/Convert_h5_to_pb/DeepCore_model_1017.h5'

# defining loss functions

# Valerio did not define epsilon in his script despite being used in the loss functions
_Epsilon = 1e-7 #value needed for the loss functione valuation
def epsilon():
    return _Epsilon

def _to_tensor(x, dtype):
    return ops.convert_to_tensor(x, dtype=dtype)

def loss_mse_select_clipped(y_true, y_pred) :
    wei = y_true[:,:,:,:,-1:]
    pred = y_pred[:,:,:,:,:-1]
    true =  y_true[:,:,:,:,:-1]
    out =K.square(tf.clip_by_value(pred-true,-5,5))*wei
  # Valerio here multplied denominator by 4 and does not add 0.00001, which is what he did during the training, so the loss function here is not defined correctly
  # return tf.reduce_sum(out, axis=None)/(tf.reduce_sum(wei,axis=None)*4) #4=numPar
  # Fixed
  return tf.reduce_sum(out, axis=None)/(tf.reduce_sum(wei,axis=None)*5 + 0.00001) #4=numPar

def loss_ROI_crossentropy(target, output):
    epsilon_ = _to_tensor(keras.backend.epsilon(), output.dtype.base_dtype)
    # epsilon_ = keras.backend.epsilon()
    output = clip_ops.clip_by_value(output, epsilon_, 1 - epsilon_)
    # print("target shape=", target.shape)
    wei = target[:,:,:,:,-1:]
    target = target[:,:,:,:,:-1]
    # print("target shape=", target.shape, ", wei shape=", wei.shape, "output shape ", output.shape)
    output = output[:,:,:,:,:-1]
    output = math_ops.log(output / (1 - output))
    retval = nn.weighted_cross_entropy_with_logits(labels=target, logits=output, pos_weight=10)#900=works #2900=200x200, 125=30x30
    retval = retval*wei
    # print("output shape ", retval.shape)
    # ROI here will diverge since he did not add 0.00001 in denominator like he did in the training
    # return tf.reduce_sum(retval, axis=None)/(tf.reduce_sum(wei,axis=None))
    # Fixed  
    return tf.reduce_sum(retval, axis=None)/(tf.reduce_sum(wei,axis=None) + 0.00001)

def loss_ROIsoft_crossentropy(target, output):
    epsilon_ = _to_tensor(keras.backend.epsilon(), output.dtype.base_dtype)
    # epsilon_ = keras.backend.epsilon()
    output = clip_ops.clip_by_value(output, epsilon_, 1 - epsilon_)
    # print("target shape=", target.shape)
    wei = target[:,:,:,:,-1:]
    target = target[:,:,:,:,:-1]
    # print("target shape=", target.shape, ", wei shape=", wei.shape, "output shape ", output.shape)
    output = output[:,:,:,:,:-1]
    output = math_ops.log(output / (1 - output))
    retval = nn.weighted_cross_entropy_with_logits(labels=target, logits=output, pos_weight=10)#900=works #2900=200x200, 125=30x30
    retval = retval*(wei+0.01)
    ## ROIsoft with Valerio's definition: does not work since denominator goes to 0 sometimes 
    ## return tf.reduce_sum(retval, axis=None)/(tf.reduce_sum(wei,axis=None))
    ## ROIsoft fix 1
    ## return tf.reduce_sum(retval, axis=None)/(tf.reduce_sum(wei,axis=None)+0.00001)
    ## ROIsoft fix 2
    return tf.reduce_sum(retval,axis=None)/(tf.reduce_sum((wei+0.01),axis=None))

# Loading our model
# This is what Valerio wrote, it's loading ROI instead of ROIsoft so it's
# wrong, also it's not loading epsilon function
# net_model = load_model(weight_file_path,custom_objects={'loss_mse_select_clipped':loss_mse_select_clipped,'loss_ROIsoft_crossentropy':loss_ROI_crossentropy, '_to_tensor':_to_tensor})
# Fixed: if last loss function used is ROIsoft, use this
net_model = load_model(weight_file_path,custom_objects={'loss_mse_select_clipped':loss_mse_select_clipped,'loss_ROIsoft_crossentropy':loss_ROIsoft_crossentropy, '_to_tensor':_to_tensor, 'epsilon': epsilon})
# If last loss function used is ROI, use this
# net_model = load_model(weight_file_path,custom_objects={'loss_mse_select_clipped':loss_mse_select_clipped,'loss_ROI_crossentropy':loss_ROI_crossentropy, '_to_tensor':_to_tensor, 'epsilon':epsilon})

# renaming output nodes
num_output = 2
pred = [None]*num_output
pred_node_names = [None]*num_output
for i in range(num_output):
      pred_node_names[i] = "output_node"+str(i)
      pred[i] = tf.identity(net_model.outputs[i], name=pred_node_names[i])
print('output nodes names are: ', pred_node_names)

sess = tf.Session()
sess.run(tf.global_variables_initializer())

outputs = pred_node_names 

#IMPORTANT NOTE:
# For some reason output nodes are flipped, so you need to flip the order when running DeepCoreSeeGenerator.cc in CMSSW

# conversion here
constant_graph = tf.graph_util.convert_variables_to_constants(sess, sess.graph.as_graph_def(), outputs)

# adjust output file (pb) directory and name here, you should rename it after
tf.train.write_graph(constant_graph,"/storage/local/data1/gpuscratch/hichemb/DeepCore_git/DeepCore/Convert_h5_to_pb/", "constantgraph.pb", as_text=False)


print("saved pb file")
