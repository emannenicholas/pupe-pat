from collections.__init__ import deque
from numbers import Number
from typing import Set, Mapping
import sys
from unittest import TestCase


class PupilSynthTestCase(TestCase):
    def setUp(self):
        print(self.__class__.__name__ + '::' + self._testMethodName)


def getsize(obj_0):
    """Recursively iterate to sum size of object & members.

    This is being used to test the notion of decal and cutout
    not keeping their numpy arrays around but just generating
    them as needed.

    see https://stackoverflow.com/questions/449560
    """
    _seen_ids = set()

    def inner(obj):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)

        zero_depth_bases = (str, bytes, Number, range, bytearray)
        if isinstance(obj, zero_depth_bases):
            pass  # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, Set, deque)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, Mapping) or hasattr(obj, 'items'):
            size += sum(inner(k) + inner(v) for k, v in getattr(obj, 'items')())
        # Check for custom object instances - may subclass above too
        if hasattr(obj, '__dict__'):
            size += inner(vars(obj))
        if hasattr(obj, '__slots__'):  # can have __slots__ with __dict__
            size += sum(inner(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))
        return size
    return inner(obj_0)
