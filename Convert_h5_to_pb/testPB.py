import numpy as np
import tensorflow.compat.v1 as tf
import sys
tf.disable_eager_execution()

sess = tf.Session()

with tf.gfile.GFile(sys.argv[1],"rb") as f:
    graph_def = tf.GraphDef()
    graph_def.ParseFromString(f.read())
    graph = tf.Graph().as_default()
    sess.graph.as_default()

    tf.import_graph_def(graph_def,name='')

tensor_input_1 = sess.graph.get_tensor_by_name('input_1:0')
tensor_input_2 = sess.graph.get_tensor_by_name('input_2:0')
tensor_input_3 = sess.graph.get_tensor_by_name('input_3:0')
tensor_output_1 = sess.graph.get_tensor_by_name('output_node0:0')
tensor_output_2 = sess.graph.get_tensor_by_name('output_node1:0')

tensor_output = [tensor_output_1,tensor_output_2]
tfPar,tfProb = sess.run(tensor_output, {tensor_input_1: [[0]], tensor_input_2:[[0]], tensor_input_3: [30*[30*[[0]*4]]]})
print(tfProb)

#Alternatively, can test h5 model with these lines
#import keras
#keras = tf.keras
#weight_file = "../Training0628/DeepCore_model_0628.h5"
#kmodel = keras.models.load_model(weight_file,custom_objects={'loss_mse_select_clipped': (lambda x,x1 : x),'loss_ROIsoft_crossentropy': (lambda y,y1 : y), '_to_tensor': (lambda z : z), 'epsilon': 0})
#a=kmodel.predict([np.array([30*[30*[[0]*4]]],dtype=float),np.array([0],dtype=float),np.array([0],dtype=float)])
#print(a)
