__all__ = (
    "isect_segments",
    "isect_polygon",
    "isect_segments_include_segments",
    "isect_polygon_include_segments",
    "isect_segments__naive",
    "isect_polygon__naive",
    )
USE_IGNORE_SEGMENT_ENDINGS = True
USE_DEBUG = True
USE_VERBOSE = False
USE_PARANOID = False
USE_VERTICAL = True
X, Y = 0, 1
NUMBER_TYPE = 'native'

if NUMBER_TYPE == 'native':
    Real = float
    NUM_EPS = Real("1e-10")
    NUM_INF = Real(float("inf"))
elif NUMBER_TYPE == 'decimal':
    import decimal
    Real = decimal.Decimal
    decimal.getcontext().prec = 80
    NUM_EPS = Real("1e-10")
    NUM_INF = Real(float("inf"))
elif NUMBER_TYPE == 'numpy':
    import numpy
    Real = numpy.float64
    del numpy
    NUM_EPS = Real("1e-10")
    NUM_INF = Real(float("inf"))
elif NUMBER_TYPE == 'gmpy2':
    import gmpy2
    gmpy2.set_context(gmpy2.ieee(128))
    Real = gmpy2.mpz
    NUM_EPS = Real(float("1e-10"))
    NUM_INF = gmpy2.get_emax_max()
    del gmpy2
else:
    raise Exception("Type not found")

NUM_EPS_SQ = NUM_EPS * NUM_EPS
NUM_ZERO = Real(0.0)
NUM_ONE = Real(1.0)

class Event:
    __slots__ = (
        "type",
        "point",
        "segment",
        "slope",
        "span",
        ) + (() if not USE_DEBUG else (
        "other",
        "in_sweep",
        ))

    class Type:
        END = 0
        INTERSECTION = 1
        START = 2
        if USE_VERTICAL:
            START_VERTICAL = 3

    def __init__(self, type, point, segment, slope):
        assert(isinstance(point, tuple))
        self.type = type
        self.point = point
        self.segment = segment
        self.slope = slope
        if segment is not None:
            self.span = segment[1][X] - segment[0][X]

        if USE_DEBUG:
            self.other = None
            self.in_sweep = False
    def __hash__(self):
        return hash(self.point)

    def is_vertical(self):
        return self.span == NUM_ZERO

    def y_intercept_x(self, x: Real):
        if USE_VERTICAL:
            if self.is_vertical():
                return None

        if x <= self.segment[0][X]:
            return self.segment[0][Y]
        elif x >= self.segment[1][X]:
            return self.segment[1][Y]

        delta_x0 = x - self.segment[0][X]
        delta_x1 = self.segment[1][X] - x
        if delta_x0 > delta_x1:
            ifac = delta_x0 / self.span
            fac = NUM_ONE - ifac
        else:
            fac = delta_x1 / self.span
            ifac = NUM_ONE - fac
        assert(fac <= NUM_ONE)
        return (self.segment[0][Y] * fac) + (self.segment[1][Y] * ifac)

    @staticmethod
    def Compare(sweep_line, this, that):
        if this is that:
            return 0
        if USE_DEBUG:
            if this.other is that:
                return 0
        current_point_x = sweep_line._current_event_point_x
        this_y = this.y_intercept_x(current_point_x)
        that_y = that.y_intercept_x(current_point_x)
        # print(this_y, that_y)
        if USE_VERTICAL:
            if this_y is None:
                this_y = this.point[Y]
            if that_y is None:
                that_y = that.point[Y]

        delta_y = this_y - that_y

        assert((delta_y < NUM_ZERO) == (this_y < that_y))
       
        if abs(delta_y) > NUM_EPS:
            return -1 if (delta_y < NUM_ZERO) else 1
        else:
            this_slope = this.slope
            that_slope = that.slope
            if this_slope != that_slope:
                if sweep_line._before:
                    return -1 if (this_slope > that_slope) else 1
                else:
                    return 1 if (this_slope > that_slope) else -1

        delta_x_p1 = this.segment[0][X] - that.segment[0][X]
        if delta_x_p1 != NUM_ZERO:
            return -1 if (delta_x_p1 < NUM_ZERO) else 1

        delta_x_p2 = this.segment[1][X] - that.segment[1][X]
        if delta_x_p2 != NUM_ZERO:
            return -1 if (delta_x_p2 < NUM_ZERO) else 1

        return 0

    def __repr__(self):
        return ("Event(0x%x, s0=%r, s1=%r, p=%r, type=%d, slope=%r)" % (
            id(self),
            self.segment[0], self.segment[1],
            self.point,
            self.type,
            self.slope,
            ))


class SweepLine:
    __slots__ = (
        
        "intersections",
        "queue",
        "_events_current_sweep",
        "_current_event_point_x",
        "_before",
        )

    def __init__(self):
        self.intersections = {}

        self._current_event_point_x = None
        self._events_current_sweep = RBTree(cmp=Event.Compare, cmp_data=self)
        self._before = True

    def get_intersections(self):
        if Real is float:
            return list(self.intersections.keys())
        else:
            return [(float(p[0]), float(p[1])) for p in self.intersections.keys()]
    def get_intersections_with_segments(self):
        if Real is float:
            return [
                (p, [event.segment for event in event_set])
                for p, event_set in self.intersections.items()
            ]
        else:
            return [
                (
                    (float(p[0]), float(p[1])),
                    [((float(event.segment[0][0]), float(event.segment[0][1])),
                      (float(event.segment[1][0]), float(event.segment[1][1])))
                     for event in event_set],
                )
                for p, event_set in self.intersections.items()
            ]
    def _check_intersection(self, a: Event, b: Event):
        if ((a is None or b is None) or
                (a.type == Event.Type.INTERSECTION) or
                (b.type == Event.Type.INTERSECTION)):

            return

        if a is b:
            return
        
        p = isect_seg_seg_v2_point(
                a.segment[0], a.segment[1],
                b.segment[0], b.segment[1])
       
        if p is None:
            return
        if USE_IGNORE_SEGMENT_ENDINGS:
            if ((len_squared_v2v2(p, a.segment[0]) < NUM_EPS_SQ or
                 len_squared_v2v2(p, a.segment[1]) < NUM_EPS_SQ) and
                (len_squared_v2v2(p, b.segment[0]) < NUM_EPS_SQ or
                 len_squared_v2v2(p, b.segment[1]) < NUM_EPS_SQ)):

                return
        events_for_point = self.intersections.pop(p, set())
        is_new = len(events_for_point) == 0
        events_for_point.add(a)
        events_for_point.add(b)
        self.intersections[p] = events_for_point
        
        if is_new and p[X] >= self._current_event_point_x:
            event_isect = Event(Event.Type.INTERSECTION, p, None, None)
            self.queue.offer(p, event_isect)

    def _sweep_to(self, p):
        if p[X] == self._current_event_point_x:
            return

        self._current_event_point_x = p[X]

    def insert(self, event):
        assert(event not in self._events_current_sweep)
        assert(not USE_VERTICAL or event.type != Event.Type.START_VERTICAL)
        if USE_DEBUG:
            assert(event.in_sweep == False)
            assert(event.other.in_sweep == False)

        self._events_current_sweep.insert(event, None)

        if USE_DEBUG:
            event.in_sweep = True
            event.other.in_sweep = True

    def remove(self, event):
        try:
            self._events_current_sweep.remove(event)
            if USE_DEBUG:
                assert(event.in_sweep == True)
                assert(event.other.in_sweep == True)
                event.in_sweep = False
                event.other.in_sweep = False
            return True
        except KeyError:
            if USE_DEBUG:
                assert(event.in_sweep == False)
                assert(event.other.in_sweep == False)
            return False

    def above(self, event):
        return self._events_current_sweep.succ_key(event, None)

    def below(self, event):
        return self._events_current_sweep.prev_key(event, None)
    
    def above_all(self, event):
        return self._events_current_sweep.key_slice(event, None, reverse=False)

    def handle(self, p, events_current):
        if len(events_current) == 0:
            return
        assert(p[0] == self._current_event_point_x)

        if not USE_IGNORE_SEGMENT_ENDINGS:
            if len(events_current) > 1:
                for i in range(0, len(events_current) - 1):
                    for j in range(i + 1, len(events_current)):
                        self._check_intersection(
                                events_current[i], events_current[j])

        for e in events_current:
            self.handle_event(e)

    def handle_event(self, event):
        t = event.type
        if t == Event.Type.START:
            # print("  START")
            self._before = False
            self.insert(event)

            e_above = self.above(event)
            e_below = self.below(event)

            self._check_intersection(event, e_above)
            self._check_intersection(event, e_below)
            if USE_PARANOID:
                self._check_intersection(e_above, e_below)

        elif t == Event.Type.END:
            # print("  END")
            self._before = True

            e_above = self.above(event)
            e_below = self.below(event)

            self.remove(event)

            self._check_intersection(e_above, e_below)
            if USE_PARANOID:
                self._check_intersection(event, e_above)
                self._check_intersection(event, e_below)

        elif t == Event.Type.INTERSECTION:
            
            self._before = True
            event_set = self.intersections[event.point]
            
            reinsert_stack = []  # Stack
            for e in event_set:
                
                if self.remove(e):
                    reinsert_stack.append(e)
            self._before = False

            
            while reinsert_stack:
                e = reinsert_stack.pop()

                self.insert(e)

                e_above = self.above(e)
                e_below = self.below(e)

                self._check_intersection(e, e_above)
                self._check_intersection(e, e_below)
                if USE_PARANOID:
                    self._check_intersection(e_above, e_below)
        elif (USE_VERTICAL and
                (t == Event.Type.START_VERTICAL)):
         
            assert(event.segment[0][X] == event.segment[1][X])
            assert(event.segment[0][Y] <= event.segment[1][Y])
           
            y_above_max = event.segment[1][Y]
            for e_above in self.above_all(event):
                if e_above.type == Event.Type.START_VERTICAL:
                    continue
                y_above = e_above.y_intercept_x(
                        self._current_event_point_x)
                if USE_IGNORE_SEGMENT_ENDINGS:
                    if y_above >= y_above_max - NUM_EPS:
                        break
                else:
                    if y_above > y_above_max:
                        break
                self._check_intersection(event, e_above)
class EventQueue:
    __slots__ = (
        "events_scan",
        )

    def __init__(self, segments, line: SweepLine):
        self.events_scan = RBTree()
        for s in segments:
            assert(s[0][X] <= s[1][X])

            slope = slope_v2v2(*s)

            if s[0] == s[1]:
                pass
            elif USE_VERTICAL and (s[0][X] == s[1][X]):
                e_start = Event(Event.Type.START_VERTICAL, s[0], s, slope)

                if USE_DEBUG:
                    e_start.other = e_start

                self.offer(s[0], e_start)
            else:
                e_start = Event(Event.Type.START, s[0], s, slope)
                e_end   = Event(Event.Type.END,   s[1], s, slope)

                if USE_DEBUG:
                    e_start.other = e_end
                    e_end.other = e_start

                self.offer(s[0], e_start)
                self.offer(s[1], e_end)

        line.queue = self

    def offer(self, p, e: Event):
        existing = self.events_scan.setdefault(
                p, ([], [], [], []) if USE_VERTICAL else
                   ([], [], []))
        existing[e.type].append(e)
    def poll(self):
        assert(len(self.events_scan) != 0)
        p, events_current = self.events_scan.pop_min()
        return p, events_current


def isect_segments_impl(segments, include_segments=False) -> list:
    if Real is float:
        segments = [
            (s[0], s[1]) if (s[0] <= s[1]) else
            (s[1], s[0])
            for s in segments]
    else:
        segments = [
            (
                (Real(s[0][0]), Real(s[0][1])),
                (Real(s[1][0]), Real(s[1][1])),
            ) if (s[0] <= s[1]) else
            (
                (Real(s[1][0]), Real(s[1][1])),
                (Real(s[0][0]), Real(s[0][1])),
            )
            for s in segments]

    sweep_line = SweepLine()
    queue = EventQueue(segments, sweep_line)

    while len(queue.events_scan) > 0:
        if USE_VERBOSE:
            print(len(queue.events_scan), sweep_line._current_event_point_x)
        p, e_ls = queue.poll()
        for events_current in e_ls:
            if events_current:
                sweep_line._sweep_to(p)
                sweep_line.handle(p, events_current)

    if include_segments is False:
        return sweep_line.get_intersections()
    else:
        return sweep_line.get_intersections_with_segments()


def isect_polygon_impl(points, include_segments=False) -> list:
    n = len(points)
    segments = [
        (tuple(points[i]), tuple(points[(i + 1) % n]))
        for i in range(n)]
    return isect_segments_impl(segments, include_segments=include_segments)


def isect_segments(segments) -> list:
    return isect_segments_impl(segments, include_segments=False)


def isect_polygon(segments) -> list:
    return isect_polygon_impl(segments, include_segments=False)


def isect_segments_include_segments(segments) -> list:
    return isect_segments_impl(segments, include_segments=True)


def isect_polygon_include_segments(segments) -> list:
    return isect_polygon_impl(segments, include_segments=True)

def slope_v2v2(p1, p2):
    if p1[X] == p2[X]:
        if p1[Y] < p2[Y]:
            return NUM_INF
        else:
            return -NUM_INF
    else:
        return (p2[Y] - p1[Y]) / (p2[X] - p1[X])


def sub_v2v2(a, b):
    return (
        a[0] - b[0],
        a[1] - b[1])


def dot_v2v2(a, b):
    return (
        (a[0] * b[0]) +
        (a[1] * b[1]))


def len_squared_v2v2(a, b):
    c = sub_v2v2(a, b)
    return dot_v2v2(c, c)


def line_point_factor_v2(p, l1, l2, default=NUM_ZERO):
    u = sub_v2v2(l2, l1)
    h = sub_v2v2(p, l1)
    dot = dot_v2v2(u, u)
    return (dot_v2v2(u, h) / dot) if dot != NUM_ZERO else default


def isect_seg_seg_v2_point(v1, v2, v3, v4, bias=NUM_ZERO):
    if v1 > v2:
        v1, v2 = v2, v1
    if v3 > v4:
        v3, v4 = v4, v3

    if (v1, v2) > (v3, v4):
        v1, v2, v3, v4 = v3, v4, v1, v2

    div = (v2[0] - v1[0]) * (v4[1] - v3[1]) - (v2[1] - v1[1]) * (v4[0] - v3[0])
    if div == NUM_ZERO:
        return None

    vi = (((v3[0] - v4[0]) *
           (v1[0] * v2[1] - v1[1] * v2[0]) - (v1[0] - v2[0]) *
           (v3[0] * v4[1] - v3[1] * v4[0])) / div,
          ((v3[1] - v4[1]) *
           (v1[0] * v2[1] - v1[1] * v2[0]) - (v1[1] - v2[1]) *
           (v3[0] * v4[1] - v3[1] * v4[0])) / div,
          )

    fac = line_point_factor_v2(vi, v1, v2, default=-NUM_ONE)
    if fac < NUM_ZERO - bias or fac > NUM_ONE + bias:
        return None

    fac = line_point_factor_v2(vi, v3, v4, default=-NUM_ONE)
    if fac < NUM_ZERO - bias or fac > NUM_ONE + bias:
        return None
    return vi

def isect_segments__naive(segments) -> list:
    isect = []
    
    if Real is float:
        segments = [
            (s[0], s[1]) if s[0][X] <= s[1][X] else
            (s[1], s[0])
            for s in segments]
    else:
        segments = [
            (
                (Real(s[0][0]), Real(s[0][1])),
                (Real(s[1][0]), Real(s[1][1])),
            ) if (s[0] <= s[1]) else
            (
                (Real(s[1][0]), Real(s[1][1])),
                (Real(s[0][0]), Real(s[0][1])),
            )
            for s in segments]

    n = len(segments)

    for i in range(n):
        a0, a1 = segments[i]
        for j in range(i + 1, n):
            b0, b1 = segments[j]
            if a0 not in (b0, b1) and a1 not in (b0, b1):
                ix = isect_seg_seg_v2_point(a0, a1, b0, b1)
                if ix is not None:
                    isect.append(ix)

    return isect


def isect_polygon__naive(points) -> list:
    isect = []

    n = len(points)

    if Real is float:
        pass
    else:
        points = [(Real(p[0]), Real(p[1])) for p in points]


    for i in range(n):
        a0, a1 = points[i], points[(i + 1) % n]
        for j in range(i + 1, n):
            b0, b1 = points[j], points[(j + 1) % n]
            if a0 not in (b0, b1) and a1 not in (b0, b1):
                ix = isect_seg_seg_v2_point(a0, a1, b0, b1)
                if ix is not None:

                    if USE_IGNORE_SEGMENT_ENDINGS:
                        if ((len_squared_v2v2(ix, a0) < NUM_EPS_SQ or
                             len_squared_v2v2(ix, a1) < NUM_EPS_SQ) and
                            (len_squared_v2v2(ix, b0) < NUM_EPS_SQ or
                             len_squared_v2v2(ix, b1) < NUM_EPS_SQ)):
                            continue

                    isect.append(ix)

    return isect


from operator import attrgetter
_sentinel = object()


class _ABCTree(object):
    def __init__(self, items=None, cmp=None, cmp_data=None):
        self._root = None
        self._count = 0
        if cmp is None:
            def cmp(cmp_data, a, b):
                if a < b:
                    return -1
                elif a > b:
                    return 1
                else:
                    return 0
        self._cmp = cmp
        self._cmp_data = cmp_data
        if items is not None:
            self.update(items)

    def clear(self):
        def _clear(node):
            if node is not None:
                _clear(node.left)
                _clear(node.right)
                node.free()
        _clear(self._root)
        self._count = 0
        self._root = None

    @property
    def count(self):
        return self._count

    def get_value(self, key):
        node = self._root
        while node is not None:
            cmp = self._cmp(self._cmp_data, key, node.key)
            if cmp == 0:
                return node.value
            elif cmp < 0:
                node = node.left
            else:
                node = node.right
        raise KeyError(str(key))

    def pop_item(self):
        if self.is_empty():
            raise KeyError("pop_item(): tree is empty")
        node = self._root
        while True:
            if node.left is not None:
                node = node.left
            elif node.right is not None:
                node = node.right
            else:
                break
        key = node.key
        value = node.value
        self.remove(key)
        return key, value
    popitem = pop_item

    def min_item(self):
        if self.is_empty():
            raise ValueError("Tree is empty")
        node = self._root
        while node.left is not None:
            node = node.left
        return node.key, node.value

    def max_item(self):
        if self.is_empty():
            raise ValueError("Tree is empty")
        node = self._root
        while node.right is not None:
            node = node.right
        return node.key, node.value

    def succ_item(self, key, default=_sentinel):
        node = self._root
        succ_node = None
        while node is not None:
            cmp = self._cmp(self._cmp_data, key, node.key)
            if cmp == 0:
                break
            elif cmp < 0:
                if (succ_node is None) or self._cmp(self._cmp_data, node.key, succ_node.key) < 0:
                    succ_node = node
                node = node.left
            else:
                node = node.right

        if node is None:
            if default is _sentinel:
                raise KeyError(str(key))
            return default
       
        if node.right is not None:
            node = node.right
            while node.left is not None:
                node = node.left
            if succ_node is None:
                succ_node = node
            elif self._cmp(self._cmp_data, node.key, succ_node.key) < 0:
                succ_node = node
        elif succ_node is None:
            if default is _sentinel:
                raise KeyError(str(key))
            return default
        return succ_node.key, succ_node.value

    def prev_item(self, key, default=_sentinel):
        node = self._root
        prev_node = None

        while node is not None:
            cmp = self._cmp(self._cmp_data, key, node.key)
            if cmp == 0:
                break
            elif cmp < 0:
                node = node.left
            else:
                if (prev_node is None) or self._cmp(self._cmp_data, prev_node.key, node.key) < 0:
                    prev_node = node
                node = node.right

        if node is None:
            if default is _sentinel:
                raise KeyError(str(key))
            return default
        if node.left is not None:
            node = node.left
            while node.right is not None:
                node = node.right
            if prev_node is None:
                prev_node = node
            elif self._cmp(self._cmp_data, prev_node.key, node.key) < 0:
                prev_node = node
        elif prev_node is None:
            if default is _sentinel:
                raise KeyError(str(key))
            return default
        return prev_node.key, prev_node.value

    def __repr__(self):

        tpl = "%s({%s})" % (self.__class__.__name__, '%s')
        return tpl % ", ".join(("%r: %r" % item for item in self.items()))

    def __contains__(self, key):
       
        try:
            self.get_value(key)
            return True
        except KeyError:
            return False

    def __len__(self):

        return self.count

    def is_empty(self):
        
        return self.count == 0

    def set_default(self, key, default=None):
       
        try:
            return self.get_value(key)
        except KeyError:
            self.insert(key, default)
            return default
    setdefault = set_default

    def get(self, key, default=None):
     
        try:
            return self.get_value(key)
        except KeyError:
            return default

    def pop(self, key, *args):
        
        if len(args) > 1:
            raise TypeError("pop expected at most 2 arguments, got %d" % (1 + len(args)))
        try:
            value = self.get_value(key)
            self.remove(key)
            return value
        except KeyError:
            if len(args) == 0:
                raise
            else:
                return args[0]

    def prev_key(self, key, default=_sentinel):
        
        item = self.prev_item(key, default)
        return default if item is default else item[0]

    def succ_key(self, key, default=_sentinel):
        
        item = self.succ_item(key, default)
        return default if item is default else item[0]

    def pop_min(self):
        
        item = self.min_item()
        self.remove(item[0])
        return item

    def pop_max(self):
       
        item = self.max_item()
        self.remove(item[0])
        return item

    def min_key(self):
    
        return self.min_item()[0]

    def max_key(self):
        return self.max_item()[0]

    def key_slice(self, start_key, end_key, reverse=False):
        
        return (k for k, v in self.iter_items(start_key, end_key, reverse=reverse))

    def iter_items(self,  start_key=None, end_key=None, reverse=False):

        if self.is_empty():
            return []
        if reverse:
            return self._iter_items_backward(start_key, end_key)
        else:
            return self._iter_items_forward(start_key, end_key)

    def _iter_items_forward(self, start_key=None, end_key=None):
        for item in self._iter_items(left=attrgetter("left"), right=attrgetter("right"),
                                     start_key=start_key, end_key=end_key):
            yield item

    def _iter_items_backward(self, start_key=None, end_key=None):
        for item in self._iter_items(left=attrgetter("right"), right=attrgetter("left"),
                                     start_key=start_key, end_key=end_key):
            yield item

    def _iter_items(self, left=attrgetter("left"), right=attrgetter("right"), start_key=None, end_key=None):
        node = self._root
        stack = []
        go_left = True
        in_range = self._get_in_range_func(start_key, end_key)

        while True:
            if left(node) is not None and go_left:
                stack.append(node)
                node = left(node)
            else:
                if in_range(node.key):
                    yield node.key, node.value
                if right(node) is not None:
                    node = right(node)
                    go_left = True
                else:
                    if not len(stack):
                        return  # all done
                    node = stack.pop()
                    go_left = False

    def _get_in_range_func(self, start_key, end_key):
        if start_key is None and end_key is None:
            return lambda x: True
        else:
            if start_key is None:
                start_key = self.min_key()
            if end_key is None:
                return (lambda x: self._cmp(self._cmp_data, start_key, x) <= 0)
            else:
                return (lambda x: self._cmp(self._cmp_data, start_key, x) <= 0 and
                        self._cmp(self._cmp_data, x, end_key) < 0)


class Node(object):
    
    __slots__ = ['key', 'value', 'red', 'left', 'right']

    def __init__(self, key=None, value=None):
        self.key = key
        self.value = value
        self.red = True
        self.left = None
        self.right = None

    def free(self):
        self.left = None
        self.right = None
        self.key = None
        self.value = None

    def __getitem__(self, key):
        
        return self.left if key == 0 else self.right

    def __setitem__(self, key, value):
    
        if key == 0:
            self.left = value
        else:
            self.right = value


class RBTree(_ABCTree):
    
    @staticmethod
    def is_red(node):
        if (node is not None) and node.red:
            return True
        else:
            return False

    @staticmethod
    def jsw_single(root, direction):
        other_side = 1 - direction
        save = root[other_side]
        root[other_side] = save[direction]
        save[direction] = root
        root.red = True
        save.red = False
        return save

    @staticmethod
    def jsw_double(root, direction):
        other_side = 1 - direction
        root[other_side] = RBTree.jsw_single(root[other_side], other_side)
        return RBTree.jsw_single(root, direction)

    def _new_node(self, key, value):
      
        self._count += 1
        return Node(key, value)

    def insert(self, key, value):

        if self._root is None:
            self._root = self._new_node(key, value)
            self._root.red = False
            return

        head = Node()
        grand_parent = None
        grand_grand_parent = head
        parent = None
        direction = 0
        last = 0

        
        grand_grand_parent.right = self._root
        node = grand_grand_parent.right
        while True:
            if node is None:
                node = self._new_node(key, value)
                parent[direction] = node
            elif RBTree.is_red(node.left) and RBTree.is_red(node.right):
                node.red = True
                node.left.red = False
                node.right.red = False

            
            if RBTree.is_red(node) and RBTree.is_red(parent):
                direction2 = 1 if grand_grand_parent.right is grand_parent else 0
                if node is parent[last]:
                    grand_grand_parent[direction2] = RBTree.jsw_single(grand_parent, 1 - last)
                else:
                    grand_grand_parent[direction2] = RBTree.jsw_double(grand_parent, 1 - last)

            
            if self._cmp(self._cmp_data, key, node.key) == 0:
                node.value = value
                break

            last = direction
            direction = 0 if (self._cmp(self._cmp_data, key, node.key) < 0) else 1
            
            if grand_parent is not None:
                grand_grand_parent = grand_parent
            grand_parent = parent
            parent = node
            node = node[direction]

        self._root = head.right
        self._root.red = False

    def remove(self, key):
        
        if self._root is None:
            raise KeyError(str(key))
        head = Node()
        node = head
        node.right = self._root
        parent = None
        grand_parent = None
        found = None
        direction = 1

       
        while node[direction] is not None:
            last = direction

           
            grand_parent = parent
            parent = node
            node = node[direction]

            direction = 1 if (self._cmp(self._cmp_data, node.key, key) < 0) else 0

            
            if self._cmp(self._cmp_data, key, node.key) == 0:
                found = node

           
            if not RBTree.is_red(node) and not RBTree.is_red(node[direction]):
                if RBTree.is_red(node[1 - direction]):
                    parent[last] = RBTree.jsw_single(node, direction)
                    parent = parent[last]
                elif not RBTree.is_red(node[1 - direction]):
                    sibling = parent[1 - last]
                    if sibling is not None:
                        if (not RBTree.is_red(sibling[1 - last])) and (not RBTree.is_red(sibling[last])):
                            
                            parent.red = False
                            sibling.red = True
                            node.red = True
                        else:
                            direction2 = 1 if grand_parent.right is parent else 0
                            if RBTree.is_red(sibling[last]):
                                grand_parent[direction2] = RBTree.jsw_double(parent, last)
                            elif RBTree.is_red(sibling[1-last]):
                                grand_parent[direction2] = RBTree.jsw_single(parent, last)
                            
                            grand_parent[direction2].red = True
                            node.red = True
                            grand_parent[direction2].left.red = False
                            grand_parent[direction2].right.red = False

        
        if found is not None:
            found.key = node.key
            found.value = node.value
            parent[int(parent.right is node)] = node[int(node.left is None)]
            node.free()
            self._count -= 1

        
        self._root = head.right
        if self._root is not None:
            self._root.red = False
        if not found:
            raise KeyError(str(key))
