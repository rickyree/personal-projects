# -*- coding: utf-8 -*-
"""redditscraper_top3posts.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1PRogt8RIo1TbfbzCaTfLAgXl5R3X817F
"""

import lxml 
from bs4 import BeautifulSoup 
import request 


userinput = str(input('Which thread would you like to look up? ')) 

#continually makes requests to reddit.com until data is obtained (because sometimes the scraping request returns blank values)
success = False 
while success == False: 
  hi = requests.get('https://www.reddit.com/r/' + userinput + '/top/') 
  data = hi.text 
  soup = BeautifulSoup(data, 'lxml')
  items = [a for a in soup.select('._1poyrkZ7g36PawDueRza-J')] 

  
  #uses css class selectors to extract relevant data 
  post = [] 
  content = [] 
  votes = [] 
  comment = []
  for x in range(len(items)): 
    post.append(items[x].select('._eYtD2XCVieq6emjKBH3m')) 
    content.append(items[x].select('._292iotee39Lmt0MkQZ2hPV')) 
    votes.append(items[x].select('._1rZYMD_4xY3gRcSS3p8ODO')) 
    comment.append(items[x].select('.FHCV02u6Cp2zYL0fhQPsO')) 

  #if the request returns blank data try again  
  if not post: 
    print('Loading data...')
  else:  
    success = True 

  
#output of the top 3 posts 
print('\nThese are the top 3 posts of r/', userinput,'as of today: \n') 
for x in range(3): 
  print('Post #', x+1)
  print('Title: \n', post[x][0].text, '\n') 

  if not not content[x]: 
    print('Content: \n', content[x][0].text, '\n\n\n') 
  
  #if the post has no written content and only contains image 
  else: 
    print('Content: \n', 'No written content for this post...', '\n') 

  print(votes[x][0].text, 'votes, ', comment[x][0].text, '\n\n\n')