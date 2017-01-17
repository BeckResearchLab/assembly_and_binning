# Developed on Janet's laptop: /Users/janet/Dropbox/meta4_bins_data_and_files/170116_compare_metabat_runs 
# coding: utf-8

# In[1]:

import pandas as pd
import seaborn as sns

import matplotlib
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt


# In[2]:

nones = pd.read_csv('./abundances_of_none.tsv', sep='\t')


# In[3]:

nones.head(2)


# In[4]:

sample_info = nones[['sample id', 'week', 'oxygen', 'replicate']].drop_duplicates()


# In[5]:

sample_info.shape


# In[6]:

pivoted = nones.pivot(index='sample id', columns='metabat run', values='abundance')
pivoted.reset_index(inplace=True)
pivoted.index.rename('index', inplace=True)
print(pivoted.columns)
#del pivoted['metabat run']
pivoted.head()


# In[7]:

results = pd.merge(pivoted, sample_info, on='sample id')
results.set_index('sample id', inplace=True)
results.head()


# In[8]:

results.sort_values(by=['week', 'oxygen', 'replicate'], inplace=True, ascending=False)


# In[9]:

results.head()


# In[10]:

results[['default settings', 'specific', 'specific 2500']].plot()


# In[11]:

results.set_index('week')[['default settings', 'specific', 'specific 2500']].plot()


# In[12]:

results.set_index('week')[['default settings', 'specific', 'specific 2500']].head()


# In[13]:

results.head(2)


# In[14]:

df = results[(results['replicate'] == 1) & (results['oxygen'] == 'low')]


# In[15]:

df.head()


# In[16]:

df.set_index('week')[['default settings', 'specific', 'specific 2500']].head()


# In[17]:

sns.heatmap(df.set_index('week')[['default settings', 'specific', 'specific 2500']].T)


# In[18]:

import matplotlib.pyplot as plt


# In[19]:

def plot_abundance(dataframe, replicate, oxygen):
    df = dataframe.sort_values('week', ascending=True)
    df = df.set_index('week')[['default settings', 'specific', 'specific 2500']].T
    #print(df.head(3))
    fig, ax = plt.subplots(1,1, figsize=(8,3.5))
    sns.heatmap(df, ax=ax, vmin=0, vmax=1, annot=True)
    title = 'abundance of contigs not included in bins: replicate {}, {} oxygen'.format(replicate, oxygen)
    fig.suptitle(title)
    fig.savefig(filename = '{}_oxygen_replicate_{}'.format(oxygen, replicate) + '.pdf')
    return fig
    


# In[20]:

p = plot_abundance(df, 'a', 'x')
#p.suptitle('abc')


# In[21]:

df.head(3)


# In[22]:

for tup, sub_df in results.groupby(['oxygen', 'replicate']):
    print(tup)
    plot_abundance(sub_df, tup[1], tup[0])


# In[23]:

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


# In[24]:

fig, axs = plt.subplots(2, 1, figsize=(7, 5))
o2_dict = {'low': axs[0], 'high': axs[1]}
colors = {1:'#66c2a5', 2:'#fc8d62', 3:'#8da0cb', 4:'#e78ac3'}
for (o2, rep), df in results.groupby(['oxygen', 'replicate']):
    #print(o2)
    #print(rep)
    ax = o2_dict[o2]
    #print(sdf.head(2))
    color = colors[rep]
    #df.plot.line(x='week', y='default settings', ax=ax, c=color) #, label='replicate {}'.format(rep))
    #df.plot.scatter(x='week', y='default settings', ax=ax, c=color, label='replicate {}'.format(rep))
    ax.plot(df['week'], df['default settings'], linestyle='-', marker='o', color=color)

labels = ['rep {}'.format(n) for n in [1, 2, 3, 4]]
axs[0].set_title('low oxygen')
axs[0].set_ylabel('abundance: not binned')
axs[0].legend(labels, bbox_to_anchor=(1.15, 1.05))

axs[1].set_title('high oxygen')
axs[1].set_ylabel('abundance: not binned')
lgd = axs[1].legend(labels, bbox_to_anchor=(1.15, 1.05))

# instead of plt.tight_layout(), 
# http://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box
plt.savefig('170117_abundance_of_unbinned_contigs_by_series.pdf', 
            bbox_extra_artists=(lgd,), bbox_inches='tight')

