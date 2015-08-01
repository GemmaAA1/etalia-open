def attr_getter(attr):
    def get_any(self):
        return getattr(self, attr)
    return get_any

def attr_setter(attr, attr2):
    def set_any(self, value):
        if value < 0:
            setattr(self, attr, value)
        elif value > 1000:
            setattr(self, attr, value)
        else:
            setattr(self, attr, value)
        setattr(self.q, attr2, value)
    return set_any


class Q:
    def __init__(self, x, y):
        self.x = x
        self.y2 = y

class P:

    def __init__(self, x, y):
        self.q = Q(x, y)
        self._x = x
        self._y = y

    x = property(fget=attr_getter('_x'), fset=attr_setter('_x', 'x'))
    y = property(fget=attr_getter('_y'), fset=attr_setter('_y', 'y2'))

