import pandas as pd
import numpy as np


def classify_crossings(df, bounds):
    """
    Classifies ion membrane penetration events.
    """
    penetrations = []

    for sim_id, sim_group in df.groupby(level='simulation'):
        lower, upper = bounds[sim_id]

        for (ion_id, seg_id), seg_group in sim_group.groupby(['ion_id', 'segment_id']):
                times = seg_group.index.get_level_values('time').to_numpy()
                zs = seg_group['z'].to_numpy()

                region = np.where(zs < lower, -1, np.where(zs > upper, 1, 0))

                lower_cross_time = None
                upper_cross_time = None
                last = region[0]

                for t, curr in enumerate(region[1:], start=1):

                    if last == -1 and curr == 0:
                        lower_cross_time = times[t]

                    elif last == 1 and curr == 0:
                        upper_cross_time = times[t]

                    if lower_cross_time and curr == 1:
                        penetrations.append({
                            "simulation": sim_id,
                            "ion_id": ion_id,
                            "segment_id": seg_id,
                            "start_time": lower_cross_time,
                            "end_time": times[t],
                            "direction": "upward"
                        })
                        lower_cross_time = None

                    elif upper_cross_time and curr == -1:
                        penetrations.append({
                            "simulation": sim_id,
                            "ion_id": ion_id,
                            "segment_id": seg_id,
                            "start_time": upper_cross_time,
                            "end_time": times[t],
                            "direction": "downward"
                        })
                        upper_cross_time = None

                    last = curr

    return pd.DataFrame(penetrations)
