import pandas as pd
import numpy as np

# Read the data
data = pd.read_csv('merged_for_dim_reduction2.csv', index_col=0)
data.index = data['SAMPLE_ID']

# Get the data for the scores
scores = data.iloc[: , ~data.columns.isin(['SAMPLE_ID', 'MTHFD2', 'UQCR11', 'X19pDEL', 'SummaryMTHFD2']) ]
scores_columns = scores.columns # Save the column names

# Transform and scale the data
from sklearn.preprocessing import FunctionTransformer
transformer = FunctionTransformer(np.cbrt, validate=True)
scores = transformer.transform(scores) # Perform squate root transformation
from sklearn.preprocessing import StandardScaler
scores = StandardScaler().fit_transform(scores) # Scale the scores


# Run PHATE
import phate
embedding = phate.PHATE()
scores_transformed = embedding.fit_transform(scores)
isomapDf = pd.DataFrame(data = scores_transformed,
                        columns = ['PHATE 1', 'PHATE 2'])
isomapDf.index = data.index
isomapfinalDf = pd.concat([isomapDf,  data[[ 'MTHFD2', 'UQCR11', 'X19pDEL', 'SummaryMTHFD2']]], axis = 1)

# Visualise
import matplotlib.pyplot as plt
marker_size = 15
fig = plt.figure(figsize = (15,15))
ax = fig.add_subplot(2,2,1)
ax.set_xlabel('PHATE 1', fontsize = 15)
ax.set_ylabel('PHATE 2', fontsize = 15)
ax.set_title('19p13 Loss', fontsize = 20)
targets = [0, 1]
colors = ['r', 'b']
ax.scatter(isomapfinalDf.loc[~isomapfinalDf['X19pDEL'].isin(targets), 'PHATE 1']
           , isomapfinalDf.loc[~isomapfinalDf['X19pDEL'].isin(targets), 'PHATE 2']
           , c='grey'
           , s=marker_size)
for target, color in zip(targets,colors):
    indicesToKeep = isomapfinalDf['X19pDEL'] == target
    ax.scatter(isomapfinalDf.loc[indicesToKeep, 'PHATE 1']
               , isomapfinalDf.loc[indicesToKeep, 'PHATE 2']
               , c = color
               , s = marker_size)
ax.legend(['Status Unknown', 'No Deletion', '19p13 Deletion'])
#ax.grid()

ax = fig.add_subplot(2,2,2)
ax.set_xlabel('PHATE 1', fontsize = 15)
ax.set_ylabel('PHATE 2', fontsize = 15)
ax.set_title('UQCR11 Expression', fontsize = 20)
plt.scatter(isomapfinalDf.loc[isomapfinalDf['UQCR11'].isna(), 'PHATE 1']
            , isomapfinalDf.loc[isomapfinalDf['UQCR11'].isna(), 'PHATE 2']
            , c='grey'
            , s=marker_size)
ax.legend(['Expression Unknown'])
plt.scatter(isomapfinalDf.loc[:, 'PHATE 1'],
            isomapfinalDf.loc[:, 'PHATE 2'],
            c = isomapfinalDf['UQCR11'],
            cmap = "RdYlGn",
            s = marker_size)
plt.clim(max(isomapfinalDf['UQCR11'].min(), -3) , min(isomapfinalDf['UQCR11'].max(), 3) )
plt.colorbar()

ax = fig.add_subplot(2,2,3)
ax.set_xlabel('PHATE 1', fontsize = 15)
ax.set_ylabel('PHATE 2', fontsize = 15)
ax.set_title('MTHFD2 Expression', fontsize = 20)
plt.scatter(isomapfinalDf.loc[isomapfinalDf['MTHFD2'].isna(), 'PHATE 1']
            , isomapfinalDf.loc[isomapfinalDf['MTHFD2'].isna(), 'PHATE 2']
            , c='grey'
            , s=marker_size)
ax.legend(['Expression Unknown'])
plt.scatter(isomapfinalDf.loc[:, 'PHATE 1'],
            isomapfinalDf.loc[:, 'PHATE 2'],
            c = isomapfinalDf['MTHFD2'],
            cmap = "RdYlGn",
            s = marker_size)
plt.clim(max(isomapfinalDf['MTHFD2'].min(), -3) , min(isomapfinalDf['MTHFD2'].max(), 3) )
plt.colorbar()


ax = fig.add_subplot(2,2,4)
ax.set_xlabel('PHATE 1', fontsize = 15)
ax.set_ylabel('PHATE 2', fontsize = 15)
ax.set_title('MTHFD2 Response', fontsize = 20)
targets = [0, 1]
colors = ['r', 'b']
ax.scatter(isomapfinalDf.loc[~isomapfinalDf['SummaryMTHFD2'].isin(targets), 'PHATE 1']
           , isomapfinalDf.loc[~isomapfinalDf['SummaryMTHFD2'].isin(targets), 'PHATE 2']
           , c='grey'
           , s=marker_size)
for target, color in zip(targets,colors):
    indicesToKeep = isomapfinalDf['SummaryMTHFD2'] == target
    ax.scatter(isomapfinalDf.loc[indicesToKeep, 'PHATE 1']
               , isomapfinalDf.loc[indicesToKeep, 'PHATE 2']
               , c = color
               , s = marker_size)
ax.legend(['Response unknown', 'Predicted Unresponsive', 'Predicted Responsive'])
#ax.grid()

# Plot the stromal scores for each cell type and save it to pdf plots
scoresDf = pd.DataFrame(data = scores,
                        columns = scores_columns, index= data.index)
scoresDf = scoresDf.merge(isomapfinalDf, right_index= True, left_index=True, how = 'inner')

for i in range(0, len(scores_columns)):
    plt.figure(i+2)
    plt.scatter(scoresDf.loc[:, 'PHATE 1'],
                scoresDf.loc[:, 'PHATE 2'],
                c = scoresDf[scores_columns[1]],
                cmap = "RdYlGn",
                s = marker_size)
    plt.xlabel('ISOMAP 1', fontsize = 15)
    plt.ylabel('ISOMAP 2', fontsize = 15)
    plt.title(scores_columns[i] + " Score", fontsize = 20)
    plt.clim(max(scoresDf[scores_columns[i]].min(), -3) , min(scoresDf[scores_columns[i]].max(), 3) )
    plt.colorbar()
    plt.savefig( scores_columns[i] + ".pdf")
    plt.close(i+2)


# Make correlation plots between MTHFD2 or UQCr11 expression and the ISOMAP components
from scipy import stats

plt.figure(51)
x = isomapfinalDf.loc[:, 'PHATE 1']
y = isomapfinalDf.loc[:, 'MTHFD2']
idx = np.isfinite(x) & np.isfinite(y)
spearman = stats.spearmanr(x[idx], y[idx])
plt.scatter(x[idx], y[idx] , s = marker_size, c = 'darkred')
plt.title("Spearman Correlation = " + str(round(spearman.correlation,2)) + " \n" + f"p-value = {spearman.pvalue:.2e}")
plt.xlabel("IS0MAP 1")
plt.ylabel("MTHFD2 Expression")
m, b = np.polyfit(x[idx], y[idx], 1)
plt.plot(x[idx], m*x[idx]+b, c = 'black', ls = '-', linewidth = 1)
plt.savefig( "MTHFD2_ISOMAP1.pdf")


plt.figure(52)
x = isomapfinalDf.loc[:, 'PHATE 2']
y = isomapfinalDf.loc[:, 'MTHFD2']
idx = np.isfinite(x) & np.isfinite(y)
spearman = stats.spearmanr(x[idx], y[idx])
plt.scatter(x[idx], y[idx] , s = marker_size, c = 'darkred')
plt.title("Spearman Correlation = " + str(round(spearman.correlation,2)) + " \n" + f"p-value = {spearman.pvalue:.2e}")
plt.xlabel("PHATE 2")
plt.ylabel("MTHFD2 Expression")
m, b = np.polyfit(x[idx], y[idx], 1)
plt.plot(x[idx], m*x[idx]+b, c = 'black', ls = '-', linewidth = 1)
plt.savefig( "MTHFD2_PHATE2.pdf")


plt.figure(53)
x = isomapfinalDf.loc[:, 'PHATE 1']
y = isomapfinalDf.loc[:, 'UQCR11']
idx = np.isfinite(x) & np.isfinite(y)
spearman = stats.spearmanr(x[idx], y[idx])
plt.scatter(x[idx], y[idx] , s = marker_size, c = 'darkred')
plt.title("Spearman Correlation = " + str(round(spearman.correlation,2)) + " \n" + f"p-value = {spearman.pvalue:.2e}")
plt.xlabel("IS0MAP 1")
plt.ylabel("UQCR11 Expression")
m, b = np.polyfit(x[idx], y[idx], 1)
plt.plot(x[idx], m*x[idx]+b, c = 'black', ls = '-', linewidth = 1)
plt.savefig( "UQCR11_ISOMAP1.pdf")



plt.figure(54)
x = isomapfinalDf.loc[:, 'ISOMAP 2']
y = isomapfinalDf.loc[:, 'UQCR11']
idx = np.isfinite(x) & np.isfinite(y)
spearman = stats.spearmanr(x[idx], y[idx])
plt.scatter(x[idx], y[idx] , s = marker_size, c = 'darkred')
plt.title("Spearman Correlation = " + str(round(spearman.correlation,2)) + " \n" + f"p-value = {spearman.pvalue:.2e}")
plt.xlabel("IS0MAP 2")
plt.ylabel("UQCR11 Expression")
m, b = np.polyfit(x[idx], y[idx], 1)
plt.plot(x[idx], m*x[idx]+b, c = 'black', ls = '-', linewidth = 1)
plt.savefig( "UQCR11_ISOMAP2.pdf")