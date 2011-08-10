#!/usr/bin/env python

def f(func):
    print "do stuff"
    return func

@f
def foo():
    print "do foo()"
    return 1

# f(foo)
foo()
