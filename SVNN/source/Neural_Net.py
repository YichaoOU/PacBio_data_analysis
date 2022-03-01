#!/usr/bin/env python3

import pandas as pd
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle


def drawPlot(y_pred, y_test):
    cnf_matrix = metrics.confusion_matrix(y_test, y_pred)

    class_names=[0,1] # name  of classes
    fig, ax = plt.subplots()
    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names)
    plt.yticks(tick_marks, class_names)
    # create heatmap
    sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
    ax.xaxis.set_label_position("top")
    plt.tight_layout()
    plt.title('Confusion matrix', y=1.1)
    plt.ylabel('Actual label')
    plt.xlabel('Predicted label')
    
    print("total number of points to predict: {}".format(len(y_pred)))
    print ("total number of reads detected as 1: {}".format(cnf_matrix[0][1] + cnf_matrix[1][1]))
    print ("detected 1 to all reads: {}".format((cnf_matrix[0][1] + cnf_matrix[1][1]) / len(y_pred)))
    print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
    print("Precision:",metrics.precision_score(y_test, y_pred))
    print("Recall:",metrics.recall_score(y_test, y_pred))
    print ("================================")    

            
colNames = ['sequenceName', 'read_start', 'read_end','s1', 's2', 'cm', 'length', 'biggest_del', 'biggest_ins', 'big_del_tol', 'big_ins_tol','insertion' \
            , 'deletion' , 'mismatch', 'AS' ,'de', 'most_dif', 'second_diff' ,'third_diff', 'NM', 'ms','label']
SV_test = pd.read_csv('temp_test.csv', header = None, names = colNames)
detected_NN = open("temp_detected_reads.txt", "w+")
featureCols = ['read_start', 'read_end','s1', 's2', 'cm', 'length','biggest_del', 'biggest_ins', 'big_del_tol', 'big_ins_tol','insertion' , 'deletion' \
               , 'mismatch', 'AS' ,'most_dif', 'second_diff' , 'third_diff','NM', 'ms']


X_test = SV_test[featureCols]
#y_test = SV_test.label
#print (featureCols)
#print (len(y_test))

filename = '/home/yli11/Programs/SVNN/model/finalized_model.sav'
loaded_model = pickle.load(open(filename, 'rb'))
y_pred_NN = loaded_model.predict(X_test)
#print(y_pred_NN)
#drawPlot(y_pred_NN, y_test)


for i in range(len(y_pred_NN)):
    if y_pred_NN[i] == 1:
        detected_NN.write(SV_test.iloc[i,0] + '\n')
    
detected_NN.close()






