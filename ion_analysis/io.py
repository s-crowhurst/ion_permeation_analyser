import pandas as pd
import numpy as np
from pathlib import Path
import MDAnalysis as mda


def read_positions(top, traj, ion_resname, simulation_id):
    """Loads ion positions from an MD trajectory into a MultiIndex DataFrame."""

    u = mda.Universe(top, traj)
    ions = u.select_atoms(f"resname {ion_resname}")

    records = []
    for ts in u.trajectory:
        for ion, pos in zip(ions, ions.positions):
            records.append({
                "simulation": simulation_id,
                "time": ts.time,
                "ion_id": ion.id,
                "x": pos[0],
                "y": pos[1],
                "z": pos[2]
            })

    df = pd.DataFrame.from_records(records)
    df.set_index(['simulation', 'time', 'ion_id'], inplace=True)

    return df


def read_multiple_positions(simulations, ion_resname):
    dfs = []
    for top, traj in simulations:
        sim_id = f"{Path(top).stem}_{Path(traj).stem}"
        df = read_positions(top, traj, ion_resname, simulation_id=sim_id)
        dfs.append(df)
    return pd.concat(dfs)