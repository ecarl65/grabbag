
def mandel(data):
    for ii in range(1000):
        im = (ii - 500)*.002
        for jj in range(1500):
            re = (jj - 1000)*.002
            xx = re
            yy = im
            cc = 0
            while 1:
                tx = xx*xx - yy*yy + re
                ty = 2*xx*yy + im
                r2 = tx*tx + ty*ty
                xx = tx
                yy = ty
                cc += 1
                if r2 > 4.0 or cc == 255:
                    data[ii*1500 + jj] = cc
                    break


if __name__ == '__main__':
    from os import write
    data = bytearray(1500*1000)
    mandel(data)
    write(1, b"P5 1500 1000 255\n")
    write(1, data)
