import time
from collections import Counter

import numpy as np
from imblearn.over_sampling import SMOTE

from ACDC import PhenoModel as Model

if __name__ == "__main__":
    model = Model()

    times, n_samples = [], []

    for n in [len(model.y0)] + [200000, 400000, 800000, 1600000]:
        if n > len(model.y0):
            sampling_strategy = dict(Counter(np.random.choice(model.y0, n)))
            sm = SMOTE(sampling_strategy=sampling_strategy, random_state=42)
            model.X0, model.y0 = sm.fit_resample(model.X0, model.y0)

        n_samples.append(n)

        start = time.perf_counter()
        model()
        times.append(time.perf_counter() - start)

        print("n:", n_samples)
        print("times:", times)
