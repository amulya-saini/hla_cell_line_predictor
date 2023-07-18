#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Ignoring warnings
import warnings
warnings.filterwarnings("ignore")


# In[2]:


#The HLA haplotype list given by the client
haplotype_list = {"haplotype":['A*01:01', 'A*02:01', 'A*03:01', 'A*11:01', 'A*23:01', 'A*24:02', 'A*26:01', 'A*29:02', 
'A*30:01', 'A*30:02', 'A*31:01', 'A*32:01', 'A*33:03', 'A*68:01', 'A*68:02', 'B*35:01', 
'B*07:02', 'B*40:01', 'B*08:01', 'B*44:03', 'B*51:01', 'B*53:01', 'B*44:02', 'B*15:01', 
'B*18:01', 'B*58:01', 'B*14:02', 'B*27:05', 'B*15:03', 'B*48:01', 'B*52:01', 'B*42:01', 
'B*49:01', 'B*57:01', 'B*13:02', 'B*39:01', 'B*40:02', 'B*58:02', 'B*46:01', 'B*35:03', 
'B*57:03', 'B*38:01', 'B*15:10', 'B*50:01', 'B*55:01', 'B*15:02', 'B*37:01', 'B*54:01', 
'B*38:02', 'B*45:01', 'B*81:01', 'B*07:05', 'B*14:01', 'B*35:02']}


# In[3]:


#importing pandas
import pandas as pd

#Converting the list to a pandas dataframe
haplotype_df = pd.DataFrame(haplotype_list)


# In[4]:


haplotype_df.head()


# In[5]:


#spliting the haplotype into two columns: locus and alleles
haplotype_df[['locus','allele']] = pd.DataFrame(haplotype_df["haplotype"].str.split('*',1).tolist())


# In[6]:


haplotype_df.head()


# In[7]:


#loading the dataset
df = pd.read_csv("hla_data_clean.csv")


# In[8]:


df.head()


# In[9]:


#Since the required HLA haplotypes list contains only A and B, The other haplotypes are being dropped from the given dataset
df.drop(['HLA-C 1', 'HLA-C 2','HLA-DQB1 1', 'HLA-DQB1 2','HLA-DRB1 1','HLA-DRB1 2'], axis=1, inplace=True)

#setting column "SampleID" as dataframe Index
df.set_index('SampleID', inplace=True)


# In[10]:


df.head()


# In[11]:


#copying the dataframe based on the column "locus"
haplotype_A_df = haplotype_df[haplotype_df.locus == "A"]


# In[12]:


haplotype_A_df.head()


# In[13]:


#Creating a list with data from "allele" column, that contains the position of HLA A haplotype 
allele_A_list = haplotype_A_df.allele.unique()


# In[14]:


allele_A_list


# In[15]:


#Creating an empty two column dataframe to store the allele location and their significant SampleID 
results_A = pd.DataFrame(columns=['allele', 'SampleID'])

#creating a for loop to iterate through every haplotype location in the list and the rows and get the index (SampleID)
for i in allele_A_list:
    index_list = df.index[(df['HLA-A 1'] == i) | (df['HLA-A 2'] == i) ].tolist()
    if len(index_list) > 0:
        new_row = {'allele': i, 'SampleID': index_list}
        results_A = results_A.append(new_row, ignore_index=True)


# In[16]:


results_A.head()


# In[18]:


#Creating a dictionary with SampleID as values and the allele as index
alleles_dict_A = pd.Series(results_A.SampleID.values,index = results_A.allele).to_dict()


# In[19]:


#defining a function to find the minimum of number of Smaples required to buy
def find_min_sample_ids(alleles_dict):
    
    # creating an empty list to hold the sample IDs
    selected_samples = []
    
    # Creating a for loop to iterate through each allele in the alleles_dict dictionary
    for allele in alleles_dict:
        
        # finding the sample ids where the allele is present
        allele_samples = set(alleles_dict[allele])
        
        # Check if any of these sample IDs are already in the selected_samples list
        common_samples = allele_samples.intersection(selected_samples)
        
        if len(common_samples) > 0:
            continue
        
        # Selecting the sample ID that covers the most alleles
        best_sample = None
        best_coverage = 0
        for sample in allele_samples:
            
            # Finding the set of alleles covered by this sample ID
            sample_alleles = set([a for a in alleles_dict if sample in alleles_dict[a]])
            
            # Counting the number of remaining alleles covered by this sample ID
            coverage = len(sample_alleles - set(selected_samples))
            
            # Updating the best sample by sample if this one covers more alleles
            if coverage > best_coverage:
                best_sample = sample
                best_coverage = coverage
        
        # Adding the best sample to the selected_samples list
        selected_samples.append(best_sample)
    
    # Return the selected_samples list
    return selected_samples


# In[20]:


# Calling the find_min_sample_ids function
min_sample_ids_A = find_min_sample_ids(alleles_dict_A)

# Printing the minimum list of sample IDs that covers all of required alleles
print('No of min samples required to cover "A" samples is : ', len(min_sample_ids_A))
print(min_sample_ids_A)


# In[21]:


haplotype_B_df = haplotype_df[haplotype_df.locus == "B"]


# In[22]:


allele_B_list = haplotype_B_df.allele.unique()


# In[23]:


results_B = pd.DataFrame(columns=['allele', 'SampleID'])

for i in allele_B_list:
    index_list = df.index[(df['HLA-B 1'] == i) | (df['HLA-B 2'] == i) ].tolist()
    if len(index_list) > 0:
        new_row = {'allele': i, 'SampleID': index_list}
        results_B = results_B.append(new_row, ignore_index=True)


# In[24]:


results_B.shape


# In[25]:


#Creating a dictionary with SampleID as values and the allele as index
alleles_dict_B = pd.Series(results_B.SampleID.values,index = results_B.allele).to_dict()


# In[26]:


# Calling the find_min_sample_ids function
min_sample_ids_B = find_min_sample_ids(alleles_dict_B)

# Printing the minimum list of sample IDs that covers all of required alleles
print('No of min samples required to cover "B" samples is : ', len(min_sample_ids_B))
print(min_sample_ids_B)


# In[27]:


#Combining the two lists of SampleID and removing duplicates to obtain minimum samples required to get all the haplotypes from the given list
min_sample_IDs = list(set(min_sample_ids_A + min_sample_ids_B))
len(min_sample_IDs), min_sample_IDs

