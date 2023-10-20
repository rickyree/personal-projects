## Introduction 
Hi, I'm Rick. 

I'm a uni student studying biomed who codes for fun. 
Here are a few projects that I have been working on.


## Projects 

### - Website scraper - 
Two projects involving scraping data from the online forum Reddit. 

#### Files 
redditscraper_notifier.py : notifies the user when there is a new post in r/minecraft, then the program terminates (r/minecraft is chosen as there is a lot of new posts per time). 

redditscraper_top3posts.py : returns the title, content, votes and comments of today's top 3 posts of a reddit thread chosen by the user. For example, enter 'minecraft' to obtain a summary of the top 3 posts of r/minecraft. 


### - Income Predictor - 
Linear classification to predict whether a person's income is above or below 50K/year (uses tensorflow). Data obtained from the UCI machine learning database (https://archive.ics.uci.edu/ml/machine-learning-databases/adult/), which was extracted from the American census bureau database. 

#### Files 

adult.data : database that contains various features such as occupation, ethnicity, gender etc., and if a person makes above or below 50K/year. This is the dataset used for training the model. 

adult.test : database used for testing the model. 

predict_income.py : linear classification model to predict a person's income. Uses adult.data to train and adult.test to evaulate. 


### - Height Analyzer - 
Analyzes heights of 18 yr olds in different countries (source: Our World in Data). Takes number of countries to compare and countries as input, and outputs a line plot showing height changes of each country over time, as well as a boxplot comparing people's heights in different countries. 

#### Files 

main.csv: database containing mean height for each country, year (1985 ~ 2019) and age. 

height_analyzer.py : python code that processes main.csv to produce a boxplot and a line plot according to user's input of countries. 


### - R Stat Project -
A markdown notebook in R that I created for the STAT module in my course. Explores a fictional anti-viral drug VK001's therapeutic effects on different types of viruses. Demonstrates my skills in utiliizing R as a tool in data analysis and presentation. 

#### Files 

STAT_ICA1_21.rdata : database containing viral load of different types of viruses according to VK001 treatment. My job is to figure out for which virus VK001 is effective against. 

hl2320_STAT.Rmd : Contains the data analysis process in a markdown file in Rmd and html format. Analyses STAT_ICA1_21.rdata. 


### - Coronavirus Data Analyzer -
An algorithm that analyzes worldwide Covid-19 cases. 

The databank from Johns Hopkins Coronavirus Resource Center is analyzed. An example database file is included, but any file from any date in this link can also be used: 
https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_daily_reports  


#### Files 

Corona-12-10-2020.csv : database containing Covid-19 cases of different countries. 

Coronavirus data analyzer.Rmd : Outputs summary statistics comparing Covid-19 cases of various countries according to user input. Can use Corona-12-10-2020.csv or any other files from the link above. Run this file and answer questions in the console for customized analysis of cases in chosen countries. 



