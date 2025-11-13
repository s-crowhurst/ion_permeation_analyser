import pandas as pd
import numpy as np
from pathlib import Path
import MDAnalysis as mda


def read_positions(top, traj, ion_resname, simulation_id, stride=1, chunk_size=100):
    """Loads ion positions from an MD trajectory into a MultiIndex DataFrame."""
    
    u = mda.Universe(top, traj)
    ions = u.select_atoms(f"resname {ion_resname}")

    if len(ions) == 0:
        raise ValueError(f"No ions with resname '{ion_resname}' found in {traj}")

    chunks = []
    frame_data = []

    for i, ts in enumerate(u.trajectory[::stride]):
        
        pos = ions.positions.copy()
        ion_ids = ions.ids

        df = pd.DataFrame({
            "simulation": simulation_id,
            "time": ts.time,
            "ion_id": ion_ids,
            "x": pos[:, 0],
            "y": pos[:, 1],
            "z": pos[:, 2]
        })
        frame_data.append(df)

        # periodically flush to reduce memory
        if (i + 1) % chunk_size == 0:
            chunks.append(pd.concat(frame_data, ignore_index=True))
            frame_data.clear()

    if frame_data:
        chunks.append(pd.concat(frame_data, ignore_index=True))

    df = pd.concat(chunks, ignore_index=True)
    df.set_index(["simulation", "time", "ion_id"], inplace=True)

    return df


def read_multiple_positions(simulations, ion_resname, stride=1, chunk_size=100):
    """
    Combine multiple simulation outputs safely.
    """
    dfs = []
    for top, traj in simulations:
        sim_id = f"{Path(top).stem}_{Path(traj).stem}"
        df = read_positions(top, traj, ion_resname, sim_id, stride=stride, chunk_size=chunk_size)
        dfs.append(df)
    return pd.concat(dfs)
