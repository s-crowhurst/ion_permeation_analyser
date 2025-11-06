"""
Example usage of the ion_analysis package.
Reads multiple trajectories, segments ion trajectories,
detects membrane boundaries, and classifies crossing events.
"""

from pathlib import Path

from ion_analysis.io import read_multiple_positions
from ion_analysis.segmentation import segment_trajectories
from ion_analysis.membrane import determine_membrane_bounds
from ion_analysis.crossings import classify_crossings


# List simulations
sims = [
    ("top1.tpr", "traj1.xtc"),
    ("top2.tpr", "traj2.xtc"),
    ("top3.tpr", "traj3.xtc"),
    ("top4.tpr", "traj4.xtc"),
    ("top5.tpr", "traj5.xtc"),
]


# Load positions for all simulations
print("Loading ion positions...")
positions_df = read_multiple_positions(sims, ion_resname="K")
print("Positions loaded:")
print(positions_df.head())


# Segment trajectories to remove PBC jumps
print("\nSegmenting trajectories...")
seg_df = segment_trajectories(positions_df)
print("Segments:")
print(seg_df.head())


# Filter out tiny segments for efficiency
filtered_df = seg_df.groupby(["simulation", "ion_id", "segment_id"]).filter(lambda x: len(x) >= 5)
print("\nFiltered segments:")
print(filtered_df.head())

# Compute membrane boundaries for each simulation
print("\nDetermining membrane bounds...")
membrane_bounds = determine_membrane_bounds(filtered_df, axis='z')
print("Membrane bounds:")
print(membrane_bounds)


# Classify membrane crossing events
print("\nClassifying crossings...")
crossings_df = classify_crossings(filtered_df, membrane_bounds)
print(crossings_df.head(30))

# Summary: count crossings per simulation
crossings_per_sim = crossings_df.groupby('simulation').size()
print("\nCrossings per simulation:")
print(crossings_per_sim)
