import pandas as pd
import numpy as np

def segment_trajectories(positions_df, axis='z', jump_threshold=None, std_thresh=2):
    """
    Split ion trajectories at PBC jumps along a given axis.
    """
    # Compute displacement between consecutive frames for each ion
    disp = positions_df[axis].groupby(level=["simulation", "ion_id"]).diff().fillna(0)

    segments = []

    # Process each ion independently
    for (sim, ion_id), ion_data in positions_df.groupby(level=["simulation", "ion_id"]):
        ion_disp = disp.loc[ion_data.index]

        # Use the jump_threshold or calculate it from the given standard deviation threshold
        thresh = jump_threshold or (np.abs(ion_disp).mean() + std_thresh * np.abs(ion_disp).std())

        jumps = np.abs(ion_disp.values) > thresh
        segment_ids = np.cumsum(jumps)

        seg_data = ion_data.copy()
        seg_data['segment_id'] = segment_ids

        segments.append(seg_data)

    # Combine all ions
    seg_df = pd.concat(segments)
    seg_df.set_index('segment_id', append=True, inplace=True)
    seg_df = seg_df.reorder_levels(['simulation', 'ion_id', 'segment_id', 'time'])

    return seg_df

