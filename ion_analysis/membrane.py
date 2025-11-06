import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

def determine_membrane_bounds(positions_df, axis='z', grid_points=1000):
    """
    Identifies membrane/channel boundaries as the first inward
    zero-crossings of the second derivative of the ion-density KDE
    from each solvent peak.
    """

    def find_zero_crossings(x, y):
        """Return x-positions where y crosses zero (linear interpolation)."""
        sign_changes = np.where(np.diff(np.sign(y)) != 0)[0]
        z = []
        for i in sign_changes:
            x0, x1 = x[i], x[i+1]
            y0, y1 = y[i], y[i+1]
            z.append(x0 - y0 * (x1 - x0) / (y1 - y0))
        return np.array(z)

    bounds = {}

    for sim, sim_df in positions_df.groupby(level="simulation"):
        coords = sim_df[axis].values

        # KDE smoothing
        kde = gaussian_kde(coords)
        x = np.linspace(coords.min(), coords.max(), grid_points)
        density = kde(x)

        # Derivatives
        peaks, _ = find_peaks(density)
        if len(peaks) < 2:
            raise ValueError(f"Simulation {sim}: expected 2 peaks, found {len(peaks)}.")

        largest_two = np.argsort(density[peaks])[-2:]
        peak_positions = np.sort(x[peaks[largest_two]])
        left_peak, right_peak = peak_positions

        first_deriv = np.gradient(density, x)
        second_deriv = np.gradient(first_deriv, x)

        # Zero crossings of second derivative
        zeros = find_zero_crossings(x, second_deriv)
        if len(zeros) < 2:
            raise ValueError(f"Simulation {sim}: could not find zero crossings.")

        # First inward zero crossing from each peak
        left_candidates = zeros[zeros > left_peak]
        right_candidates = zeros[zeros < right_peak]

        if len(left_candidates) == 0 or len(right_candidates) == 0:
            raise ValueError(f"Simulation {sim}: boundary not found.")

        left_boundary = left_candidates.min()
        right_boundary = right_candidates.max()

        bounds[sim] = (left_boundary, right_boundary)

    return bounds

