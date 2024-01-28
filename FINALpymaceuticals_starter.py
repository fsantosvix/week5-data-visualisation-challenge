#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ## *Analysis*
# 
# 
# #### <u>Observation 1</u>
# - The mice submitted to the Capomulin regimen demonstrated a high correlation between weight and tumour volume, meaning that the heavier the mouse, the larger the tumour.
# 
# #### <u>Observation 2</u>
# - Even though the tumour in mouse l509 increased its size at the beginning of the treatment with the Capomulin drug, it rapidly decreased after day 20, which indicates a positive reaction to the medicine. However, around day 35, the tumour started increasing again. Further measurements are required to assess more accurate results.
# 
# #### <u>Observation 3</u>
# - Only one potential outlier was identified in the dataset for the four more promising treatment regimens. This represents that the data we are dealing with is of good quality, and the results can be considered reliable.
#  

# In[37]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import sem
import numpy as np


# In[38]:


# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"


# In[39]:


# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)


# In[40]:


# Display the mouse_metadata dataframe to identify potential keys for merging
mouse_metadata.head(2)


# In[41]:


# Display the study_results dataframe to identify potential keys for merging
study_results.head(2)


# In[42]:


# Combine the dataframes into a single dataframe
merged_df = pd.merge(study_results, mouse_metadata, on='Mouse ID')

# Display the data table for preview
merged_df.head()


# In[43]:


# Display the number of unique mice.
len(merged_df['Mouse ID'].unique())


# In[44]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_ID = merged_df[merged_df.duplicated(['Mouse ID','Timepoint'])]
duplicate_ID['Mouse ID'].unique()


# In[45]:


# Optional: Get all the data for the duplicate mouse ID. 
g989_records = merged_df.loc[merged_df['Mouse ID']== 'g989']
g989_records


# In[46]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_merged_df = merged_df.loc[merged_df['Mouse ID'] != 'g989']

# Checking the number of mice in the clean DataFrame.
len(clean_merged_df['Mouse ID'].unique())


# ## Summary Statistics

# In[47]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume.

#  1st step - group the clean dataframe by Drug Regimen
grouped_by_regimen = clean_merged_df.groupby('Drug Regimen')

# Identify the average tumour volume by Drug Regiment and store in a series.
tumor_volume_mean = grouped_by_regimen['Tumor Volume (mm3)'].mean()

# Identify the median for tumour volumes by Drug Regiment and store in a series.
tumor_volume_median = grouped_by_regimen['Tumor Volume (mm3)'].median()

# Identify the variance for tumour volumes by Drug Regiment and store in a series.
tumor_volume_variance = grouped_by_regimen['Tumor Volume (mm3)'].var()

# Identify the standard deviation for tumour volumes by Drug Regiment and store in a series.
tumor_volume_std = grouped_by_regimen['Tumor Volume (mm3)'].std()

# Identify the standard error of the mean for the tumour volumes by Drug Regiment and store in a series.
tumor_volume_sem = grouped_by_regimen['Tumor Volume (mm3)'].sem()


# In[48]:


# Assemble the resulting series into a single summary DataFrame.
summary_statistics_df = pd.DataFrame({'Mean Tumor Volume':tumor_volume_mean,
                                    'Median Tumor Volume':tumor_volume_median,
                                    'Tumor Volume Variance':tumor_volume_variance,
                                    'Tumor Volume Std Deviation':tumor_volume_std,
                                    'Tumor Volume Std Error':tumor_volume_sem,
                                    })

# Display the summary statistics dataframe
summary_statistics_df


# ## Bar and Pie Charts

# In[49]:


# Identify the number of records for each drug regimen
regimen_count = clean_merged_df['Drug Regimen'].value_counts()
regimen_count


# In[50]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.

# Specify the series to be plotted
regimen_count.plot(kind='bar')

# Define the label for the x axis
plt.xlabel('Drug Regimen')

# Define label for the y axis
plt.ylabel('# of observed mouse timepoints')

# Display the graph
plt.show()


# In[51]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.

# Pass the index of the dataframe as a list to be used in the x axis
x_axis = list(regimen_count.index)


# In[52]:


# Specify the type of chart and the axis
plt.bar(x_axis,regimen_count)

# Define the label for the x axis
plt.xlabel('Drug Regimen')

# Define the label for the y axis
plt.ylabel('# of observed mouse timepoints')

# Rotate the x axis labels so it can be read
plt.xticks(rotation = 'vertical')

# Display the graph
plt.show()


# In[53]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas

# Count how many of records of each sex type are in the dataset
sex_distribution = clean_merged_df['Sex'].value_counts()

# Use the index as labels for the graph
labels = list(sex_distribution.index)


# In[54]:


# Plot the graph using Pandas, identifying the type and displaying the percentage
sex_distribution.plot(kind='pie', autopct="%.1f")

# Display the graph
plt.show()


# In[55]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot

# Plot the chart, identyfing the dataset, labels and percentage format
plt.pie(sex_distribution, 
        labels = labels, 
        autopct="%.1f"
       )

# Include a label to the vertical axis
plt.ylabel('Sex')

# Display the graph
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[56]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
grouped_by_ID = clean_merged_df.groupby(['Mouse ID'])['Timepoint'].max()

# Merge the grouped df with the original DataFrame to get the tumor volume at the last timepoint
final_tumor_volume = merged_df.merge(grouped_by_ID, on=['Mouse ID', 'Timepoint'], how='right')

# Display the merged dataframe for conference
final_tumor_volume.head()


# In[57]:


# Put treatments into a list for for loop (and later for plot labels)
treatments_list = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']


# Create empty list to fill with tumor vol data (for plotting)
container_tumor_data = []


# Loop through the merged dataframe to populate the each empty list
for item in treatments_list:
    # Filter the tumour volumes by each item in the treatment list
    final_volume = final_tumor_volume.loc[final_tumor_volume['Drug Regimen'] == item, 'Tumor Volume (mm3)']
    # append the filtered data into the container list
    container_tumor_data.append(final_volume)
    # make the calculations to find the lower and upper quartiles
    quartile = final_volume.quantile([.25, .5, .75])
    lowerq = quartile[0.25]
    upperq = quartile[0.75]
    # calculate iqr and identify the upper and lower whiskers
    iqr = upperq - lowerq
    lower_bound = lowerq - 1.5 * iqr
    upper_bound = upperq + 1.5 * iqr
    # filter the outliers based on the calculations made
    outlier = final_volume.loc[((final_volume > upper_bound)
                                     | (final_volume < lower_bound))]
    # print the outliers
    print(f'The outliers for {item} are: {outlier}')


# In[58]:


# create the boxplot figure based on the container list with the tumor data from the 4 drug regimens
plt.boxplot(container_tumor_data, 
            labels = treatments_list, # use the list with the 4 drugs as the x axis labels
            sym='^g') # modify the colour and style of the outlier symbol
# add a label to the vertical axis
plt.ylabel('Final Tumor Volume (mm3)')
# add a title
plt.title('Most Promissing Drug Regimens')

#  display the graph
plt.show()


# ## Line and Scatter Plots

# In[59]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin


# Filter data related to the Capomulin treatment from the cleaned merged df
capomulin_df = clean_merged_df.loc[clean_merged_df['Drug Regimen'] == 'Capomulin']

# check existing Mouse IDs within the Capomulin database that can be used as an example
capomulin_df['Mouse ID'].unique()


# In[60]:


# LINE PLOT

# Filter data for mouse i509
l509_data = capomulin_df.loc[capomulin_df['Mouse ID'] == 'l509']

# Slice the Timepoint column to be used as the x axis info
x_axis_line = l509_data['Timepoint']

# Define the Tumour Volume as the y axis information
y_axis_line = l509_data['Tumor Volume (mm3)']

# Plot the line graph identifying the x axis data, y axis data, the line color and transparency 
plt.plot(x_axis_line, y_axis_line, color = 'red', alpha = 0.5)

#  Add a title
plt.title('Capomulin results - Mouse l509')

# Add a label to x and y axis
plt.xlabel('Timepoint')
plt.ylabel('Tumor Volume (mm3)')

# Display the chart
plt.show()


# In[61]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen

# group the data by Mouse ID and calculate the mean
grouped_cap_id = capomulin_df.groupby(['Mouse ID']).mean()

# Specify the Weight column as the x axis
x_scatter = grouped_cap_id['Weight (g)']

# Specify the Tumour Volume column as the y axis
y_scatter = grouped_cap_id['Tumor Volume (mm3)']

# Create list to be used as the x axis values
array = list(x_scatter.unique())
array.sort()
array

# Plot the scatter plot
plt.scatter(x_scatter,y_scatter, color = 'black', alpha=0.3)

# Add a grid for better visualisation
plt.grid()

# Define the created list as the x axis legend
plt.xticks(array)

# Add labels to the axis
plt.xlabel('Weight(g)')
plt.ylabel('Average Tumor Volume (mm3)')

# Add title
plt.title('Capomulin regimen data')

# Display the chart
plt.show()


# ## Correlation and Regression

# In[62]:


# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen


# Calculate the correlation between the variables
correlation = st.pearsonr(x_scatter,y_scatter)
# Print the result
print(f'The correlation between the mouse weight and the average observed tumor volume for the entire Capomulin regimen is {round(correlation[0],2)}.')

# Find the relevant information to build the linear regression line
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(x_scatter,y_scatter)

# Calculate values for the regression line
regression_line = x_scatter * slope + intercept

# Recycle the code to bring the recently created scatter plot
plt.scatter(x_scatter,y_scatter, color = 'black', alpha=0.3)
plt.xticks(array)
plt.xlabel('Weight(g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.title('Capomulin regimen data')

# Plot the linear graph together with the code from the scatter plot.
plt.plot(x_scatter, regression_line)

# Display the graphs together
plt.show()


# In[ ]:




