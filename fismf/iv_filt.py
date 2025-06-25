#!/usr/bin/env python3

from scipy.signal import butter, filtfilt

def butter_filter(data, fcc, fss, tp, order=4):
    nyq = 0.5 * fss  # Nyquist frequency
    normal_cutoff = fcc / nyq
    b, a = butter(order, normal_cutoff, btype=tp, analog=False)
    y = filtfilt(b, a, data)
    return y

