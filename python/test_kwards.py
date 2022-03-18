#!/usr/bin/env python

def myfunc(arg1, arg2, arg3=True, arg4=4, **kw):
    '''Testing the kwargs function'''

    print(f"arg1 = {arg1}")
    print(f"arg2 = {arg2}")
    print(f"arg3 = {arg3}")
    print(f"arg4 = {arg4}")

    print(kw)

    for k, v in kw.items():
        print(f"input argument {k} with value {v}")

if __name__ == '__main__':
    myfunc(1, 2,  arg3=False, arg4=5, arg6='/path/to/make', arg7=True, arg8=8.12)

