import json
import requests
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import numpy

response = requests.get("https://covidtracking.com/api/v1/us/daily.json")
todos = json.loads(response.text)

counter = 0
dates=[]
pos=[]
tot=[]
for todo in todos:
    counter = counter + 1
    dates.append(todo["date"])
    pos.append(todo["positive"])
    tot.append(todo["negative"])
    if counter > 70: break # get at most 70 data points

date0 = dates[0]
for i in range(len(dates)):
    date = dates[i]
    dates[i] = datetime(date/10000,date/100%100,date%100)

logpos = [0]*len(pos)
logtot = [0]*len(pos)
datefit = [i for i in range(len(pos))]
for i in range(len(pos)-1):
    pos[i] = pos[i] - pos[i+1]
    tot[i] = tot[i] - tot[i+1] + pos[i]
    logpos[i] = numpy.log(pos[i])
    logtot[i] = numpy.log(tot[i])

lfit = 16;
cpos=numpy.polyfit(datefit[0:lfit],logpos[0:lfit],1)
ctot=numpy.polyfit(datefit[0:lfit],logtot[0:lfit],1)
posfit = [0]*len(pos)
totfit = [0]*len(pos)
for i in range(len(datefit)):
    posfit[i]= numpy.exp(cpos[0]*datefit[i] + cpos[1])
    totfit[i]= numpy.exp(ctot[0]*datefit[i] + ctot[1])
plt.semilogy(dates[:-2],pos[:-2])
plt.semilogy(dates[:-2],tot[:-2])
plt.semilogy(dates[0:lfit],totfit[0:lfit])
plt.semilogy(dates[0:lfit],posfit[0:lfit])
plt.text(dates[15],posfit[0]/1.8,str(format(-cpos[0]*100, '.2f'))+"% per day")
plt.text(dates[15],totfit[0]/1.8,str(format(-ctot[0]*100, '.2f'))+"% per day")
plt.xlim([datetime(2020,03,20),datetime(2020,date0/100%100,date0%100)])
plt.ylim([5000,500000])
plt.xticks([datetime(2020,3,16),datetime(2020,4,1),datetime(2020,4,16),datetime(2020,5,1),datetime(2020,5,16)],['2020/3/16','2020/4/1','2020/4/16','2020/5/1','2020/5/16'])
plt.xlabel('Date')
plt.ylabel('Daily New Cases')
plt.grid(which='both')
plt.savefig('DailyNewCases')
