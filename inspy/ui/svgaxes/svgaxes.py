# coding: utf-8
from xml.dom import minidom

def read_axes(filename):

    doc = minidom.parse(filename)
    get = lambda k: doc.getElementsByTagName(k)

    svg = get('svg')[0]
    sw, sh = [float(svg.getAttribute(k)) for k in ['width', 'height']]
    axes = {}

    for r in doc.getElementsByTagName('rect'):
        id = r.getAttribute('id')
        x, y, w, h = [float(r.getAttribute(k)) for k in ['x', 'y', 'width', 'height']]
        axes[id] = (x/sw, y/sh, w/sw, h/sh)

    return sw, sh, axes
        
