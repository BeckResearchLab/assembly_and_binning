# modified from /Users/janet/Dropbox/meta4_bins_data_and_files/170118_read_mappings_by_sample/plot_frac_mapped.py  
# coding: utf-8

print('import packages...')

# In[1]:

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

print('done importing packages...')


# In[2]:

info = pd.read_csv('./data/total_num_reads_across_samples_with_sample_info.tsv',
                  sep='\t')


# In[3]:

sample_info = info[['sample id', 'week', 'oxygen', 'replicate']].drop_duplicates()
sample_info.shape


# In[4]:

info.head()


# In[5]:

df = info.copy() 
df['oxygen, rep'] = df['oxygen'].map(str) + ' oxygen, rep ' + df['replicate'].map(str)  # 


# In[6]:

df.head()


# In[7]:

df['frac mapped to contigs (upper-bound)'] = df['reads mapped (includes multiply mapped)']/df['number of reads']


# In[8]:
print('plotting time')

fig, axs = plt.subplots(2, 1, figsize=(10, 4))
print('fig, axs initialized')
axs_dict = {'low':axs[0], 'high':axs[1]}
axs[0].set_title('low oxygen')
axs[1].set_title('high oxygen')
color_dict = {1:'#bdc9e1', 2:'#74a9cf', 3:'#2b8cbe', 4:'#045a8d'}

print('loop through groupby dataframes')
for (o2, rep), sdf in df.groupby(['oxygen', 'replicate']):
    #print(sdf.head(1))
    #print(rep)
    ax = axs_dict[o2]
    sdf.sort_values('week', ascending=True)
    color = color_dict[rep]
    ax.plot(sdf['week'], sdf['frac mapped to contigs (upper-bound)'],
                linestyle='-', marker='o', color=color)
    
labels = ['rep {}'.format(n) for n in [1, 2, 3, 4]]
axs[0].legend(labels, bbox_to_anchor=(1.15, 1.05))
lgd = axs[1].legend(labels, bbox_to_anchor=(1.15, 1.05))
plt.gcf().suptitle('170118 upper bound frac of reads aligned to contigs', size=15)
fig.subplots_adjust(top=0.85)

print('save fig')
plt.savefig('170118_approx_frac_reads_mapping_to_contigs.pdf', 
            bbox_extra_artists=(lgd,), bbox_inches='tight')

