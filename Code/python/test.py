__author__ = 'koorosh'
def answer(x):
    # your code here
    out = 1
    for ix in range(1, x+1):
        out = out + 7**ix
    return out

print(answer(2))