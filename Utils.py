import matplotlib.pyplot as plt
import matplotlib as mat 
import seaborn as sns
import pandas as pd
import importlib
"""
df = pd.read_csv("dataset_luminex.csv" , sep=';')

df['class'] = df['class'].astype('category')

y = df['class'].cat.codes

print(df.shape)

fig, axes = plt.subplots(5, 10, figsize=(60, 30))
cpt_i = 0
cpt_j = 0

for i in df:
	
	if i != "class" and i != "ID":
		cpt_i +=1
		if cpt_i == 5:
			cpt_i=0
			cpt_j+=1
		sns.boxplot(ax=axes[cpt_i, cpt_j], data=df, x='class', y=i)
		
		plt.savefig("all.png")

"""

df2 = pd.read_csv("data_facs.csv" , sep=';')

df2['class'] = df2['class'].astype('category')

fig, axes = plt.subplots(3, 4, figsize=(30, 15))
cpt_i = 0
cpt_j = 0

for i in df2:
	
	if i != "class" and i != "ID":
		cpt_i +=1
		if cpt_i == 3:
			cpt_i=0
			cpt_j+=1
		sns.boxplot(ax=axes[cpt_i, cpt_j], data=df2, x='class', y=i)
		
		plt.savefig("all_facs.png")
	
	
df3 = pd.read_csv("Dataset_clinic.csv" , sep=';')

df3['class'] = df3['class'].astype('category')
fig, axes = plt.subplots(4, 10, figsize=(60, 30))
cpt_i = 0
cpt_j = 0

for i in df3:
	
	if i != "class" and i != "ID":
		cpt_i +=1
		if cpt_i == 4:
			cpt_i=0
			cpt_j+=1
		sns.boxplot(ax=axes[cpt_i, cpt_j], data=df3, x='class', y=i)
		
		plt.savefig("all_clinic.png")

print(df2.shape)
print(df3.shape)



