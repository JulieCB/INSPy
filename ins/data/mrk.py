"""
Anywave marker files

- comment lines begin with "\\"
- each line is "{tag}\t{trigger}\t{time}"
    where trigger is an integer and time is a float

We try to do some parsing of the tag as well for

- point or osc
- channel identifier as [A-Z]+\d+

"""

import re
import numpy
import pandas

from ins.util import parse_marker_tag

class Mrk(object):

    def __init__(self, mrkfile):

        self.mrkfile = mrkfile
        self.items = []
        self.tags = []
        self.triggers = []
        self.times = []
        self.parses = []

        with open(mrkfile, 'r') as fd:
            self.text = fd.read()

        for line in self.text.split('\n'):
            line = line.strip()
            if not line.startswith('//'):
                parts = line.split('\t')
                if len(parts) == 3:
                    tag, s_trigger, s_time = line.split('\t')
                    self.tags.append(tag)
                    self.triggers.append(int(s_trigger))
                    self.times.append(float(s_time))
                    parse = parse_marker_tag(tag)
                    self.parses.append(parse)
                else:
                    print 'line %r ignored' % (line,)

        self.triggers = numpy.array(self.triggers, numpy.int32)
        self.times = numpy.array(self.times)
