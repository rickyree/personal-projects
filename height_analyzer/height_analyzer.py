# -*- coding: utf-8 -*-
"""height_analyzer.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1eTWOswpDo0TCg9017EyfgzskwcUqtQrs
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from google.colab import files 
uploaded = files.upload()


raw = pd.read_csv('main.csv')

#too much data. Data is filtered to only include 18 yr old boys 
temp = []
for x in range(raw.shape[0]):
   if raw.iloc[x,3] == 18: 
     if raw.iloc[x,1] == 'Boys': 
      temp.append(x) 


mainyear = raw.take(temp)

temp = []
for x in range(raw.shape[0]):
   if raw.iloc[x,2] == 2019: 
     if raw.iloc[x,1] == 'Boys': 
      temp.append(x) 

mainage = raw.take(temp) 

print(mainyear.info(), mainage.info()) 

#drops null values 
mainyear = mainyear.dropna() 
mainage = mainage.dropna()

print(mainyear.head())
print(mainage.head())

countries = []

#input for number of countries to compare 
invalidinput = True  
while invalidinput == True: 
  invalidinput = False
  number = input('How many countries do you want to compare? ')
  try: 
    int(number)
  except ValueError: 
    print('Please enter a number. ') 
    invalidinput = True

number = int(number) 
 
#input for each country 
for x in range(number): 
  invalidinput = True
  while invalidinput == True: 
    invalidinput = False 
    countryinput = input('Enter a country: ')
    countryinput = countryinput.title()
    if countryinput == 'Usa': 
      countryinput = 'United States of America' 
    if countryinput == 'Us': 
      countryinput = 'United States of America' 
    if countryinput == 'Uk': 
      countryinput = 'United Kingdom' 
    if countryinput == 'Korea': 
      countryinput = 'South Korea' 
    
    if any(countryinput == x for x in mainyear['Country'].unique()): 
      invalidinput = False
      countries.append(countryinput) 

    else: 
      invalidinput = True 
      print('Not a valid country! Try again: ')
  
print(countries)


#sorting data appropriately into input 
heightdict1 = {}

for x in range(len(countries)): 
  heightdict1[countries[x]] = []

heightdict2 = {}

for x in range(len(countries)): 
  heightdict2[countries[x]] = []

yeardict = {}

for x in range(len(countries)): 
  yeardict[countries[x]] = [] 

agedict = {}

for x in range(len(countries)): 
  agedict[countries[x]] = []


for x in range(len(countries)): 
  for y in range(mainyear.shape[0]): 
    if mainyear.iloc[y,0] == countries[x]: 
      heightdict1[countries[x]].append(mainyear.iloc[y,4])
      yeardict[countries[x]].append(mainyear.iloc[y,2]) 

  for y in range(mainage.shape[0]): 
    if mainage.iloc[y,0] == countries[x]: 
      heightdict2[countries[x]].append(mainage.iloc[y,4])
      agedict[countries[x]].append(mainage.iloc[y,3]) 



#creating array for boxplot
cheights = []

for x, y in enumerate(heightdict1.values()): 
  cheights.append(y) 



#creating arrays for lineplot
yearheights = [[] for x in range(len(countries))] 
ageheights = [[] for x in range(len(countries))]
years = [[] for x in range(len(countries))]
age = [[] for x in range(len(countries))] 


for x, y in zip(range(len(countries)), heightdict1):
  yearheights[x] = np.append(yearheights[x], heightdict1[y])
  years[x] = np.append(years[x], yeardict[y]) 
for x, y in zip(range(len(countries)), heightdict2): 
  ageheights[x] = np.append(ageheights[x], heightdict2[y])
  age[x] = np.append(age[x], agedict[y])



#plot boxplot to compare between countries 
plt.xticks([*range(1,len(countries)+1)], countries)
plt.xlabel('Countries', weight = 'bold')
plt.ylabel('Height (cm)', weight = 'bold')
plt.xticks(rotation=45) 
plt.boxplot(cheights)
plt.show()

print('\n\n')
#plot lineplot to show changes in height over the years   
for x in range(len(countries)): 
  plt.plot(years[x], yearheights[x])
plt.legend(countries)
plt.xlabel('Year', weight = 'bold')
plt.ylabel('Height (cm)', weight = 'bold') 
plt.xticks(rotation=45) 
plt.title('Height of each country by year')
plt.show()


print('\n\n')
#plot lineplot to show changes in height over age 
plt.figure().set_figheight(10) 
for x in range(len(countries)): 
  plt.plot(age[x], ageheights[x]) 
plt.legend(countries)
plt.xlabel('Age', weight = 'bold')
plt.ylabel('Height (cm)', weight = 'bold') 
plt.xticks(rotation=45) 
plt.title('Height of each country by age')
plt.show()