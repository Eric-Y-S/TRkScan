def say_hello():
    print("Hello, World!")
    
# 通过 Cython 优化计算
def fibonacci(int n):
    cdef int a = 0, b = 1, temp
    cdef int i
    for i in range(n):
        temp = a
        a = b
        b = temp + b
    return a