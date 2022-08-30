import time
from collections import Counter

import numpy as np
from imblearn.over_sampling import SMOTE

from MP import MpModel

PATH_DATA = "ACDC/AML_benchmark/"

if __name__ == "__main__":
    model = MpModel(PATH_DATA)

    times, n_samples = [], []

    for n in [len(model.Y)] + [200000, 400000, 800000, 1600000]:
        if n > len(model.Y):
            sampling_strategy = dict(Counter(np.random.choice(model.Y, n)))
            sm = SMOTE(sampling_strategy=sampling_strategy, random_state=42)
            model.data, model.Y = sm.fit_resample(model.data, model.Y)

        n_samples.append(n)

        start = time.perf_counter()
        model()
        times.append(time.perf_counter() - start)

        print("n:", n_samples)
        print("times:", times)
