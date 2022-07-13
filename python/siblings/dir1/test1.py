#!/usr/bin/env python3

from dir2.test2 import func2
import dir2
from dir2 import test2

func2(10)
dir2.test2.func2(10)
test2.func2(10)
