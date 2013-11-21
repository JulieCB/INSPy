from scipy import io

class EEG(object):
    def __init__(self, struct):

        self.struct = struct

        # pull out most stuff automatically
        for k in struct[0][0].dtype.fields.iterkeys():
            setattr(self, k, struct[0][0][k])

        # hand tuning
        self.srate = float(self.srate)
        self.labels = [s[0] for s in self.chanlocs['labels'][:, 0]]
