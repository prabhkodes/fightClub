import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import io

# Data from previous turn
csv_data = """Config,Computation,Communication
N1_R1_T32,18.220009,0.000116
N2_R2_T8,17.255191,0.726682
N2_R2_T16,11.797811,0.432212
N4_R4_T16,6.014215,0.316104
N4_R4_T32,5.184021,0.573707
N4_R64_T2,2.302217,0.688571
N8_R8_T32,2.871148,0.728856
N8_R128_T2,1.503566,0.921988
N16_R16_T32,1.710424,0.971230
"""

df = pd.read_csv(io.StringIO(csv_data), index_col='Config')

# Plot
fig, ax = plt.subplots(figsize=(12, 7))

# Plot stacked bars
# Colors: Blue for Compute, Red/Orange for Comm (Common distinction)
df.plot(kind='bar', stacked=True, ax=ax, color=['#4c72b0', '#c44e52'], width=0.7)

# Set Log Scale
ax.set_yscale('log')

# Format Y-axis to be cleaner (powers of 10 or scalar)
# User liked "10^2 scale", so we use LogFormatterMathtext
ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
# Or we can just set simple scalar labels if the range is small (0.1 to 100)
# Let's stick to the 10^x style as requested previously.

# Labels and Title
ax.set_title('Computation vs Communication Time (Log Scale)\nFiltered for Scalability Analysis', fontsize=14)
ax.set_xlabel('System Configuration', fontsize=12)
ax.set_ylabel('Time (Seconds)', fontsize=12)

# Rotate x-labels
plt.xticks(rotation=45, ha='right')

# Grid
ax.grid(axis='y', which='major', linestyle='-', alpha=0.5)
ax.grid(axis='y', which='minor', linestyle=':', alpha=0.3)

# Legend
plt.legend(title='Category')

plt.tight_layout()
plt.savefig('stacked_histogram_clean.png')

print("Plot created.")