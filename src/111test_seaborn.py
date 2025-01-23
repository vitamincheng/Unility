import numpy as np
a=np.array([1,4,5,6,20])
b=np.array([True,False,True,False,False])
print(a[b])


import os 
os.environ['OMP_NUM_THREADS']='8'
import seaborn as sns
import matplotlib.pyplot as plt
sns.get_dataset_names()
tips=sns.load_dataset("tips")
iris=sns.load_dataset("iris")
titanic=sns.load_dataset("titanic")
planets=sns.load_dataset("planets")
#sns.scatterplot(x="tip",y="total_bill",data=tips,hue="day",size="size",palette="YlGnBu")
#plt.show()
print(tips)
print(tips['tip'])
import numpy as np
#sns.histplot(tips,x="tip",hue="sex",bins=20)
#sns.barplot(x="sex",y="tip",data=tips)
sns.boxplot(x="day",y="tip",data=tips,hue="sex")

plt.show()
print(titanic)
titanic.corr()
#from icecream import ic
#ic(tips)
