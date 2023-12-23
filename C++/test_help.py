import sys
import lagmodel as lm
from lagmodel import geometries as geo
from lagmodel import internal_class as lm_intr

try:
    help(getattr(lm, sys.argv[1]))
except AttributeError:
    try:
        help(getattr(lm_intr, sys.argv[1]))
    except AttributeError:
        try:
            help(getattr(geo, sys.argv[1]))
        except AttributeError:
            exit("Try different name")
