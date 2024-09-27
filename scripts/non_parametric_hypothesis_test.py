# # Non-parametric Hypothesis Testing

# Our data:
#  - Genetic difference between two samples in the dataset, for samples with OXA-48.
#  - Population 1: All pairs but excluding pairs in Population 2
#  - Population 2: Pairwise samples from the _same_ patients

# Our hypotheses:
#  - $H_0$: Population 1 and 2 have the same distribution of SNPs
#  - $H_1$: Population 1 and 2 have different distributions of SNPs

# Our (non-parametric) tests:
#  - Wilcoxon Signed-Rank Test
#  - Mann-Whitney Test





# Import packages:
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
sns.set_theme(style="white", font="Arial")


# Load data
df = pd.read_csv("data/genetic_distance_table.csv")
print(df.head())

df.loc[:,'same_sample'] = df.pl1_masked == df.pl2_masked
df = df[(df.pl1_masked != 'Reference') & (df.pl2_masked != 'Reference') & (df.same_sample != 1)].rename(columns = {'within_patient':'population'}).copy()


# Mean and standard deviation of each population:
print("Population 1 (not pairs)")
print(f"Mean: {df.loc[df.population == 0, 'SNPs'].mean():.2f}")
print(f"Standard Deviation: {df.loc[df.population == 0, 'SNPs'].std():.2f}")

print("Population 2 (pairs)")
print(f"Mean: {df.loc[df.population == 1, 'SNPs'].mean():.2f}")
print(f"Standard Deviation: {df.loc[df.population == 1, 'SNPs'].std():.2f}")


# Data visualisation:
fig, ax = plt.subplots(1,2)
fig.set_size_inches(10, 5)
sns.histplot(data=df.loc[df.population == 0, ['SNPs']], x='SNPs', binwidth=1, ax=ax[0])
sns.histplot(data=df.loc[df.population == 1, ['SNPs']], x='SNPs', binwidth=1, ax=ax[1])
ax[0].set_title("Population 1")
ax[1].set_title("Population 2")
plt.show()
fig.savefig('figures/distributions_of_genetic_distance.png', dpi=fig.dpi)

fig, ax = plt.subplots(1,2)
fig.set_size_inches(10, 5)
sns.histplot(data=df.loc[(df.population == 0) & (df.SNPs < 8), ['SNPs']], x='SNPs', binwidth=1, ax=ax[0])
sns.histplot(data=df.loc[(df.population == 1) & (df.SNPs < 8), ['SNPs']], x='SNPs', binwidth=1, ax=ax[1])
ax[0].set_title("Population 1 (8 SNP cutoff)")
ax[1].set_title("Population 2 (8 SNP cutoff)")
plt.show()
fig.savefig('figures/distributions_of_genetic_distance_up_to_8_SNPs.png', dpi=fig.dpi)



# Mann-Whitney Test:
mw_statistic, mw_p_value = stats.mannwhitneyu(df.loc[df.population == 0, 'SNPs'], df.loc[df.population == 1, 'SNPs'])
print("MANN-WHITNEY TEST")
# Display the results
print(f'Mann-Whitney U statistic: {mw_statistic}')
print(f'P-value: {mw_p_value}')


# Wilcoxon Signed-Rank Test:
print("WILCOXON SIGNED-RANK TEST")
# Perform Wilcoxon signed-rank test
wsr_statistic, wsr_p_value, wsr_median, wsr_table = stats.median_test(df.loc[df.population == 0, 'SNPs'], df.loc[df.population == 1, 'SNPs'])
wsr_table_df = pd.DataFrame(wsr_table, 
                            columns=pd.MultiIndex.from_tuples([('Population', 1), ('Population', 2)]),
                            index=['Above Median', 'Below Median'])
# Display the results
print(f'Wilcoxon Signed-Rank statistic: {wsr_statistic}')
print(f'P-value: {wsr_p_value}')
print(f'Median: {wsr_median}')
print(f"Contingency Table:")
print(wsr_table_df)


# Repeating analysis just for the N West
region_data = pd.read_csv("data/finalEpiDataRef.csv")[["MOLIS", "Region"]]

NW_df = (
    pd.read_csv("data/genetic_distance_table.csv")
    .merge(region_data.rename({"MOLIS": "pl1", "Region": "pl1_region"}, axis=1), on="pl1", how="left")
    .merge(region_data.rename({"MOLIS": "pl2", "Region": "pl2_region"}, axis=1), on="pl2", how="left")
)
NW_df.loc[:,'same_sample'] = NW_df.pl1 == NW_df.pl2
NW_df = NW_df.loc[(NW_df.pl1_region == "N WEST") & (NW_df.pl2_region == "N WEST")]
NW_df = NW_df.loc[(NW_df.pl1_pat != 'Reference') & (NW_df.pl2_pat != 'Reference') & (NW_df.same_sample != 1)].rename(columns = {'within_patient':'population'}).copy()

print("FOR NORTH WEST ONLY")

# Mann-Whitney Test:
NW_mw_statistic, NW_mw_p_value = stats.mannwhitneyu(NW_df.loc[NW_df.population == 0, 'SNPs'], NW_df.loc[NW_df.population == 1, 'SNPs'])
# Display the results
print("MANN-WHITNEY TEST (NW)")
print(f'Mann-Whitney U statistic: {NW_mw_statistic}')
print(f'P-value: {NW_mw_p_value}')




# Wilcoxon Signed-Rank Test:
print("WILCOXON SIGNED-RANK TEST (NW)")
# Perform Wilcoxon signed-rank test
NW_wsr_statistic, NW_wsr_p_value, NW_wsr_median, NW_wsr_table = stats.median_test(NW_df.loc[NW_df.population == 0, 'SNPs'], NW_df.loc[NW_df.population == 1, 'SNPs'])
NW_wsr_table_df = pd.DataFrame(NW_wsr_table, 
                            columns=pd.MultiIndex.from_tuples([('Population', 1), ('Population', 2)]),
                            index=['Above Median', 'Below Median'])
# Display the results
print(f'Wilcoxon Signed-Rank statistic: {NW_wsr_statistic}')
print(f'P-value: {NW_wsr_p_value}')
print(f'Median: {NW_wsr_median}')
print(f"Contingency Table:")
NW_wsr_table_df
