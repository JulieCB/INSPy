
import sys
import os
import os.path
import numpy

from IPython.utils.traitlets import CUnicode
import pandas

from ins import core
from ins.viz import tsplot

class Ades(core.Data):
    """
    This imports an "Anywave descriptor" file as well sa the
    associated data file. 

    When reading the header, we assume that samplingRate
    and numberOfSamples are specified before the channel types

    """

    path = CUnicode()

    def __init__(self, path=None):

        filepath = path

        self.dir = os.path.dirname(filepath)
        self.adesname = os.path.basename(filepath)
        self.datname = os.path.join(self.dir, 
            self.adesname[:-5] + '.dat')

        self.read_config()
        self.load_data()

    def read_config(self):
        with open(os.path.join(self.dir, self.adesname)) as fd:
            ades = fd.readlines()

        fs, nsamp, label, type = False, False, [], []
        for line in ades:
            if line[0] is not '#' and '=' in line:
                key, value = [p.strip() for p in line.strip().split('=')]
                if not (fs and nsamp):
                    if key == 'samplingRate':
                        fs = float(value)
                    elif key == 'numberOfSamples':
                        nsamp = int(value)
                    else:
                        msg = 'Unexpected: %r %r' % (key, value)
                        print msg
                else:
                    label.append(key)
                    type.append(value)

        self.fs = fs
        self._nsamp = nsamp
        self.chan = pandas.DataFrame({'label': label, 'modality': type})

    def load_data(self):
        """
        CGB's MATLAB code stores using

            fwrite(h, reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials), 'float32')

        """

        nchan = self.chan.shape[0]
        msg = "reading nchan=%d x nsamp=%d from %s..."
        msg %= (nchan, self._nsamp, self.datname)
        sys.stdout.write(msg)
        sys.stdout.flush()
        data = numpy.fromfile(self.datname, 'f').reshape((self._nsamp, nchan))
        self.data = pandas.DataFrame(data)
        sys.stdout.write(' done!\n')
   
    def to_tsplot(self):
        labs = [lab for lab, _ in self.chan]
        return tsplot.TSPlot(self.data.T, labs, self.fs)



