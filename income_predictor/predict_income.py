# -*- coding: utf-8 -*-
"""Untitled3.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1F2BRT2VGHUlF2MMiUk3Sj-Op05jm4Jjo
"""

# Commented out IPython magic to ensure Python compatibility.
# %tensorflow_version 2.x
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import clear_output
import tensorflow as tf

from google.colab import files 
traindata = files.upload()
testdata = files.upload()

#naming columns - unnecessary columns labeled as r0 through r3 
columns = ['age', 'workclass', 'r0', 'education', 'r1', 'marital_status', 'occupation', 'relationship', 'ethnicity', 'gender', 'r2', 'r3', 'hrs_per_week', 'native_country', 'income'] 


#read csv files into panda dataframe
from io import StringIO 
traindata = pd.read_csv(StringIO(traindata['adult.data'].decode('utf-8')), header = None, names = columns) 
testdata = pd.read_csv(StringIO(testdata['adult.test'].decode('utf-8')), header = None, names = columns) 

#removing 1st row of testdata which contains title 
testdata.drop(0, inplace = True) 
testdata = testdata.reset_index()

print(traindata.head(), testdata.head())
print(traindata.describe(), testdata.describe())

'''
Data Cleaning 
'''

#concatenating training and testing data as one for data cleaning, will be separated later 
fulldata = pd.concat([traindata, testdata]) 


#dropping repetitive columns
fulldata = fulldata.drop(['r0', 'r1', 'r2', 'r3'], axis = 1)


#labels columns should not be categorical, replace with numerical 
#I spent hours trying to debug this !! Especially as some elements had a period after them 
fulldata['income'] = fulldata['income'].replace({' <=50K': 0, ' >50K': 1})
fulldata['income'] = fulldata['income'].replace({' <=50K.': 0, ' >50K.': 1})


#checking that each column contains elements of the same dtype 
for x in fulldata.columns: 
  print('dtypes of each element: \n',x, pd.Series([type(y) for y in fulldata[x]]).value_counts()) 

#as age contains a mix of int and str, convert all to int 
fulldata['age'] = pd.to_numeric(fulldata['age']) 
  


#checking for empty values 
print('\n\nEmpty values of each column: \n') 
for x in fulldata.columns : 
  print(pd.isna(fulldata[x]).value_counts()) 


#separting columns into categorical and numerical 
numerical_var = ['age', 'hrs_per_week'] 
categorical_var = ['workclass', 'education', 'marital_status', 'occupation', 'relationship', 'ethnicity', 'gender', 'native_country']

#filling in empty values - mean value for numerical, most occurring for categorical 
for x in numerical_var:
  fulldata[x].fillna(fulldata[x].mean()) 

for x in categorical_var: 
  fulldata[x].fillna(fulldata[x].value_counts().index[0]) 


#separating dataset into training and testing data again 
traindata = fulldata[0: len(traindata)] 
testdata = fulldata[len(traindata):len(fulldata)] 

#separating features and labels
traineval = traindata.pop('income')
testeval = testdata.pop('income')

#creating feature columns 
feature_columns = [] 

for x in categorical_var: 
  vocabulary = []
  vocabulary = traindata[x].unique()
  feature_columns.append(tf.feature_column.categorical_column_with_vocabulary_list(key = x, vocabulary_list= vocabulary))

for x in numerical_var: 
  feature_columns.append(tf.feature_column.numeric_column(key = x, dtype = tf.float32)) 

print(feature_columns)

#creating input function 
def input_fn(features, label, epochs = 20, training = True, batch = 32): 
  def input_fn_inner(): 
    input = tf.data.Dataset.from_tensor_slices((dict(features), label))
    if training: 
      input = input.shuffle(1000)
    input = input.batch(batch).repeat(epochs) 
    return input
  return input_fn_inner 

train_input_fn = input_fn(traindata, traineval) 
test_input_fn = input_fn(testdata, testeval, training = False)

#training and evaluating the model
model = tf.estimator.LinearClassifier(feature_columns = feature_columns)
model.train(test_input_fn) 
print(model.evaluate(test_input_fn))

#takes user input of candidate number to predict if salary is above/below 50K 
results_list = list(model.predict(test_input_fn)) 
subject_num = int(input('Please enter participant number you want to predict for: '))
classid = results_list[subject_num]['class_ids']
print(results_list[subject_num]['probabilities'][classid]) 

print(traindata.iloc[subject_num])
if classid == 0: 
  print("\nThis person's income is predicted to be below 50K. ") 
if classid == 1: 
  print("This person's income is predicted to be above 50K. ") 

if testeval[subject_num] == 0: 
  print('Actual income is below 50K. ') 
if testeval[subject_num] == 1: 
  print('Actual income is above 50K. ') 

if classid == testeval[subject_num]: 
  print('Prediction is correct! ')
if classid != testeval[subject_num]: 
  print('Maybe next time... ')


#accuracy of model 
misprediction = 0
total = 0
for x in range(len(testeval)): 
  total = total + 1
  if testeval[x] != results_list[x]['class_ids']: 
    misprediction = misprediction + 1
  

print('\nAccuracy: ', (total - misprediction) / total)